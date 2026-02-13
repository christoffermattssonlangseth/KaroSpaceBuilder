from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any
import json

SCHEMA_VERSION = "1.0.0"


@dataclass(slots=True)
class AssetFile:
    kind: str
    path: str
    size_bytes: int


@dataclass(slots=True)
class ExpressionChunk:
    path: str
    start: int
    end: int
    size_bytes: int


@dataclass(slots=True)
class ExpressionGene:
    gene: str
    index: int
    chunks: list[ExpressionChunk] = field(default_factory=list)


@dataclass(slots=True)
class ExpressionSpec:
    format: str
    dtype: str
    n_cells: int
    gzip: bool
    genes: list[ExpressionGene] = field(default_factory=list)


@dataclass(slots=True)
class ImageSpec:
    path: str
    width: int
    height: int


@dataclass(slots=True)
class Manifest:
    dataset_name: str
    n_cells: int
    n_genes_exported: int
    coordinate_source: str
    annotation_columns: list[str]
    gene_list: list[str]
    asset_files: list[AssetFile]
    expression: ExpressionSpec
    schema_version: str = SCHEMA_VERSION
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    preview_path: str | None = None
    image: ImageSpec | None = None

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        if payload.get("image") is None:
            payload.pop("image", None)
        if payload.get("preview_path") is None:
            payload.pop("preview_path", None)
        return payload


def write_manifest(path: Path, manifest: Manifest) -> None:
    path.write_text(json.dumps(manifest.to_dict(), indent=2), encoding="utf-8")
