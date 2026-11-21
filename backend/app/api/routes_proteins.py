"""Federated discovery across UniProt, RCSB PDB, and AlphaFold DB."""
from __future__ import annotations

import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from fastapi import APIRouter, HTTPException, Query

router = APIRouter(prefix="/api/v1/proteins", tags=["Protein search"])
USER_AGENT = "DockLab/1.0 (public biological data search)"
_CACHE: dict[tuple[str, int, bool], tuple[float, dict]] = {}
_CACHE_LOCK = Lock()


def _json_request(url: str, payload: dict | None = None, timeout: int = 18):
    data = json.dumps(payload).encode("utf-8") if payload is not None else None
    headers = {"User-Agent": USER_AGENT, "Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    request = Request(url, data=data, headers=headers, method="POST" if data is not None else "GET")
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def _uniprot_search(term: str, limit: int) -> list[dict]:
    query = f"accession:{term}" if re.fullmatch(r"[A-Z0-9]{6,10}", term.upper()) else term
    params = urlencode({"query": query, "format": "json", "size": limit,
                        "fields": "accession,protein_name,gene_names,organism_name,length,xref_pdb"})
    data = _json_request(f"https://rest.uniprot.org/uniprotkb/search?{params}")
    result = []
    for row in data.get("results", []):
        accession = row.get("primaryAccession")
        desc = row.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
        genes = [x.get("geneName", {}).get("value") for x in row.get("genes", []) if x.get("geneName", {}).get("value")]
        organism = row.get("organism", {}).get("scientificName", "")
        pdb_ids = sorted({x.get("id", "").upper() for x in row.get("uniProtKBCrossReferences", []) if x.get("database") == "PDB"})
        result.append({"accession": accession, "name": desc or accession, "genes": genes,
                       "organism": organism, "length": row.get("sequence", {}).get("length"),
                       "pdb_ids": pdb_ids})
    return result


def _rcsb_entity_search(term: str, limit: int) -> list[str]:
    payload = {"query": {"type": "terminal", "service": "full_text", "parameters": {"value": term}},
               "return_type": "polymer_entity",
               "request_options": {"paginate": {"start": 0, "rows": limit},
                                   "results_content_type": ["experimental"]}}
    try:
        data = _json_request("https://search.rcsb.org/rcsbsearch/v2/query", payload)
    except HTTPError as exc:
        if exc.code == 204:
            return []
        raise
    return data.get("result_set", []) and [x.get("identifier", "") for x in data["result_set"] if x.get("identifier")]


def _rcsb_entity_detail(identifier: str) -> dict | None:
    try:
        pdb_id, entity_id = identifier.split("_", 1)
        entity = _json_request(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")
        entry = _json_request(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")
    except Exception:
        return None
    ids = entity.get("rcsb_polymer_entity_container_identifiers", {})
    def as_list(value):
        if isinstance(value, list): return value
        if isinstance(value, str): return [x.strip() for x in value.split(",") if x.strip()]
        return []
    info = entry.get("rcsb_entry_info", {})
    organism = entity.get("rcsb_entity_source_organism", [])
    methods = [x.get("method") for x in entry.get("exptl", []) if x.get("method")]
    resolutions = info.get("resolution_combined") or []
    return {"id": f"{pdb_id}_{entity_id}", "pdb_id": pdb_id.upper(), "entity_id": entity_id,
            "name": entity.get("entity_poly", {}).get("pdbx_description") or entry.get("struct", {}).get("title") or pdb_id.upper(),
            "title": entry.get("struct", {}).get("title", ""),
            "chain_ids": as_list(ids.get("auth_asym_ids") or ids.get("asym_ids")),
            "uniprot_ids": as_list(ids.get("uniprot_ids")),
            "organism": (organism[0].get("ncbi_scientific_name", "") if organism else ""),
            "method": ", ".join(methods), "resolution": resolutions[0] if resolutions else None,
            "kind": "experimental", "source": "RCSB PDB",
            "source_url": f"https://www.rcsb.org/structure/{pdb_id.upper()}",
            "download_url": f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"}


def _alphafold_models(accession: str) -> list[dict]:
    try:
        data = _json_request(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
    except Exception:
        return []
    records = data if isinstance(data, list) else [data]
    result = []
    for row in records[:3]:
        model_id = row.get("modelEntityId") or row.get("entryId") or f"AF-{accession}-F1"
        result.append({"id": model_id, "accession": accession,
                       "name": row.get("uniprotDescription") or row.get("proteinName") or f"Predicted model {accession}",
                       "organism": row.get("organismScientificName", ""),
                       "chain_ids": [], "uniprot_ids": [accession], "kind": "predicted",
                       "source": "AlphaFold DB · predicted, not experimental",
                       "source_url": f"https://alphafold.ebi.ac.uk/entry/{model_id}",
                       "download_url": row.get("pdbUrl") or row.get("cifUrl"),
                       "pdb_url": row.get("pdbUrl"), "cif_url": row.get("cifUrl"),
                       "confidence": row.get("globalMetricValue"),
                       "sequence_length": row.get("sequenceEnd") or row.get("uniprotEnd")})
    return result


def _search_cached(term: str, limit: int, include_predicted: bool):
    cache_key = (term.casefold(), limit, include_predicted)
    with _CACHE_LOCK:
        cached = _CACHE.get(cache_key)
        if cached and time.time() - cached[0] < 300:
            return cached[1]
    uniprot: list[dict] = []
    entities: list[str] = []
    errors: list[str] = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        f_uni = pool.submit(_uniprot_search, term, min(limit, 10))
        f_rcsb = pool.submit(_rcsb_entity_search, term, min(limit, 20))
        try: uniprot = f_uni.result()
        except Exception as exc: errors.append(f"UniProt: {exc.__class__.__name__}")
        try: entities = f_rcsb.result()
        except Exception as exc: errors.append(f"RCSB PDB: {exc.__class__.__name__}")
    seen_ids = set(entities)
    for protein in uniprot:
        for pdb_id in protein["pdb_ids"]:
            # The RCSB search handles names; UniProt cross-references fill in known structure links.
            if re.fullmatch(r"[A-Z0-9]{4}", pdb_id):
                seen_ids.add(f"{pdb_id}_1")
    structures = []
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures = [pool.submit(_rcsb_entity_detail, identifier) for identifier in list(seen_ids)[:limit]]
        for future in as_completed(futures):
            value = future.result()
            if value:
                structures.append(value)
    # Put high resolution experimental models first, when resolution is available.
    structures.sort(key=lambda row: (row["resolution"] is None, row["resolution"] or 99, row["pdb_id"]))
    predicted = []
    if include_predicted:
        with ThreadPoolExecutor(max_workers=4) as pool:
            futures = [pool.submit(_alphafold_models, row["accession"]) for row in uniprot[:min(5, limit)]]
            for future in as_completed(futures):
                for model in future.result():
                    protein = next((x for x in uniprot if x["accession"] == model["accession"]), {})
                    model["name"] = protein.get("name") or model["name"]
                    model["organism"] = protein.get("organism") or model["organism"]
                    model["genes"] = protein.get("genes", [])
                    predicted.append(model)
    result = {"query": term, "proteins": uniprot, "structures": structures[:limit],
              "predicted_models": predicted[:limit], "sources": ["UniProtKB", "RCSB PDB"] + (["AlphaFold DB"] if include_predicted else []),
              "warnings": errors, "cached_seconds": 300}
    with _CACHE_LOCK:
        _CACHE[cache_key] = (time.time(), result)
    return result


@router.get("/search")
async def search_proteins(q: str = Query(min_length=2, max_length=120), limit: int = Query(default=10, ge=1, le=20), include_predicted: bool = True):
    term = " ".join(q.split())
    try:
        return _search_cached(term, limit, include_predicted)
    except Exception as exc:
        raise HTTPException(502, detail={"error": {"code": "PROTEIN_SEARCH_UNAVAILABLE", "message": "Protein databases did not respond. Try again or use a PDB ID.", "details": {"reason": exc.__class__.__name__}}}) from exc
