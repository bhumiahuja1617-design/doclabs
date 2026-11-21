"""Self-contained HTML and PDF report rendering."""
from __future__ import annotations
from html import escape
from typing import Any


class ReportService:
    @staticmethod
    def render_html(job: Any) -> str:
        protein = job.settings.get("protein", {})
        target = protein.get("target") or protein
        if not isinstance(target, dict):
            target = {}
        rows = job.results or []
        target_name = target.get("name") or target.get("title") or protein.get("pdb_id") or "Target not set"
        source_urls = target.get("sources", [])
        if isinstance(source_urls, str):
            source_urls = [source_urls]
        if target.get("source_url"):
            source_urls = [target["source_url"], *source_urls]
        sources = " ".join(f'<a href="{escape(str(url), quote=True)}">{escape(str(url))}</a>'
                            for url in dict.fromkeys(source_urls) if str(url).startswith(("http://", "https://")))
        facts = [
            ("Structure", target.get("structure_type") or protein.get("structure_type")),
            ("PDB / model", target.get("pdb_id") or protein.get("pdb_id")),
            ("Organism", target.get("organism")),
            ("UniProt", target.get("uniprot_id") or ", ".join(target.get("uniprot_ids", []))),
            ("Method", target.get("method") or target.get("structure_method")),
            ("Resolution", target.get("resolution") or target.get("resolution_angstroms")),
            ("Chains", ", ".join(target.get("chain_ids") or protein.get("chains", []))),
            ("Source", target.get("source") or protein.get("source")),
        ]
        fact_rows = "".join(f"<tr><th>{escape(label)}</th><td>{escape(str(value))}</td></tr>"
                             for label, value in facts if value not in (None, "", []))
        summary = target.get("function_summary") or "No functional summary is available for this structure."
        crystal_ligand = target.get("co_crystallized_ligand") or {}
        crystal_code = crystal_ligand.get("residue_name") if isinstance(crystal_ligand, dict) else None
        chemical_html = (f'<p>Co-crystallized ligand: {escape(str(crystal_ligand.get("name") or crystal_code))} '
                         f'(<a href="https://pdbj.org/chemie/pdb/{escape(str(crystal_code), quote=True)}">PDBj chemical component {escape(str(crystal_code))}</a>)</p>'
                         if crystal_code else "")
        residues = target.get("key_residues", [])
        residue_html = ("<h3>Key residues</h3><ul>" + "".join(
            f"<li>{escape(str(item.get('residue', '')))} · {escape(str(item.get('role', '')))}</li>"
            for item in residues if isinstance(item, dict)) + "</ul>") if residues else ""
        inhibitors = target.get("known_inhibitors") or target.get("approved_inhibitors") or []
        inhibitor_html = ("<p>Known inhibitors: " + ", ".join(
            escape(str(item.get("name", ""))) for item in inhibitors if isinstance(item, dict)) + "</p>") if inhibitors else ""
        target_html = (f'<h2>Selected target</h2><h3>{escape(str(target_name))}</h3>'
                       f'<p>{escape(str(summary))}</p><table><tbody>{fact_rows}</tbody></table>'
                       f'{chemical_html}{residue_html}{inhibitor_html}<p>{sources}</p>')
        qc_html = "<p>QC not run.</p>"
        if job.qc_report:
            checks = "".join(f"<tr><td>{escape(c['severity'].upper())}</td><td>{escape(c['title'])}</td><td>{escape(c['detail'])}</td></tr>" for c in job.qc_report.get("checks", []))
            qc_html = f"<p>{escape(job.qc_report.get('summary', ''))}</p><table><thead><tr><th>Status</th><th>Check</th><th>Finding</th></tr></thead><tbody>{checks}</tbody></table>"
        result_rows = "".join(f"<tr><td>{escape(str(r.get('rank','—')))}</td><td>{escape(str(r.get('name','')))}</td><td>{escape(str(r.get('label','unknown')))}</td><td>{escape(str(r.get('score','—')))}</td><td>{escape(str(r.get('status','')))}</td><td>{escape(str(r.get('error','')))}</td></tr>" for r in rows)
        results_html = f"<table><thead><tr><th>Rank</th><th>Molecule</th><th>Dataset label</th><th>Vina score</th><th>Status</th><th>Details</th></tr></thead><tbody>{result_rows}</tbody></table>"
        progress = getattr(job, "progress_done", 0) or 0
        total = getattr(job, "progress_total", 0) or 0
        interim_note = (f"<p><strong>Interim report:</strong> {progress} of {total} molecules have finished. Reopen or refresh this report for updated results.</p>"
                        if job.status in {"queued", "running", "reliability"} else "")
        labels = [str(item.get("label", "unknown")).lower() for item in (job.ligands or [])]
        active_count, decoy_count = labels.count("active"), labels.count("decoy")
        analysis = job.settings.get("analysis")
        analysis_html = (f"<p>Screening analysis: {escape(str(analysis))}</p>" if analysis else
                         "<p>Screening metrics require at least three docked actives and three docked decoys.</p>")
        return f"""<!doctype html><html><head><meta charset="utf-8"><style>
        @page{{size:A4;margin:18mm}}body{{font:11pt/1.5 Georgia,serif;color:#161616}}h1,h2{{font-weight:400}}h1{{font-size:25pt}}h2{{border-bottom:1px solid #aaa;padding-bottom:5px;margin-top:24pt}}table{{width:100%;border-collapse:collapse;font:8.5pt monospace}}th,td{{padding:5px;border-bottom:1px solid #ccc;text-align:left;vertical-align:top}}.note{{border-top:1px solid #bbb;margin-top:30pt;padding-top:8pt;font-size:9pt;color:#555}}</style></head><body>
        <h1>DockLab scientific report</h1><p>Job {escape(job.id)} · {escape(job.status)} · created {escape(job.created_at.isoformat())}</p>{interim_note}
        {target_html}<h2>Inputs and search box</h2><p>Protein: {escape(str(protein.get('pdb_id', '')))} · Molecules: {len(job.ligands or [])} · labeled actives: {active_count} · labeled decoys: {decoy_count}</p><pre>{escape(str(job.box))}</pre>
        <h2>Pre-docking QC</h2>{qc_html}<h2>Results</h2>{results_html}<h2>Screening interpretation</h2>{analysis_html}
        <h2>Settings and versions</h2><pre>{escape(str(job.settings))}</pre><p class="note">Docking scores are estimates. They do not guarantee real binding. Reliability measures repeatability of the computation, not biological activity.</p>
        </body></html>"""

    @classmethod
    def render_pdf(cls, job: Any) -> bytes:
        from weasyprint import HTML
        return HTML(string=cls.render_html(job), base_url="").write_pdf()
