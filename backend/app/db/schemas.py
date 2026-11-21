"""Request and response validation models."""
from typing import Literal
from pydantic import BaseModel, Field, ConfigDict


class LigandInput(BaseModel):
    name: str = Field(min_length=1, max_length=120)
    smiles: str = Field(min_length=1, max_length=4000)
    label: Literal["active", "decoy", "unknown"] = "unknown"


class LigandsRequest(BaseModel):
    ligands: list[LigandInput] = Field(min_length=1, max_length=200)


class BoxModel(BaseModel):
    model_config = ConfigDict(extra="forbid")
    center: dict[str, float]
    size: dict[str, float]


class SettingsRequest(BaseModel):
    exhaustiveness: int = Field(default=8, ge=1, le=32)
    seed: int = 42
    reliability: bool = False
    top_k: int = Field(default=10, ge=1, le=20)
    n_seeds: int = Field(default=5, ge=2, le=10)


class RunRequest(BaseModel):
    acknowledge_qc_failures: bool = False
