"""Small-molecule depiction endpoints."""
from fastapi import APIRouter, HTTPException, Response
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

router = APIRouter(prefix="/api/v1/ligands", tags=["Ligands"])


class DepictRequest(BaseModel):
    smiles: str = Field(min_length=1, max_length=4000)
    name: str = Field(default="", max_length=120)
    width: int = Field(default=440, ge=200, le=1000)
    height: int = Field(default=280, ge=140, le=800)


@router.post("/depict")
async def depict_ligand(payload: DepictRequest):
    mol = Chem.MolFromSmiles(payload.smiles)
    if mol is None:
        raise HTTPException(422, detail={"error": {"code": "INVALID_SMILES", "message": "This SMILES string could not be parsed into a molecular structure.", "details": {}}})
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(payload.width, payload.height)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.bondLineWidth = 2
    options.addStereoAnnotation = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    return Response(svg, media_type="image/svg+xml; charset=utf-8", headers={"Cache-Control": "private, max-age=300"})


@router.post("/summary")
async def summarize_ligand(payload: DepictRequest):
    """Return a compact chemical-component record and 2D drawing for a ligand."""
    mol = Chem.MolFromSmiles(payload.smiles)
    if mol is None:
        raise HTTPException(422, detail={"error": {"code": "INVALID_SMILES", "message": "This SMILES string could not be parsed into a molecular structure.", "details": {}}})
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(payload.width, payload.height)
    drawer.drawOptions().clearBackground = False
    drawer.drawOptions().bondLineWidth = 2
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return {"name": payload.name,
            "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 3),
            "formal_charge": sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
            "atom_count": mol.GetNumAtoms(),
            "svg": drawer.GetDrawingText().replace("svg:", "")}
