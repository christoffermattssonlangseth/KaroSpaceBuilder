from __future__ import annotations

from pathlib import Path
import gzip
import json

import pandas as pd


def file_size(path: Path) -> int:
    return path.stat().st_size


def write_json(path: Path, payload: object, gzip_enabled: bool = False) -> Path:
    text = json.dumps(payload, separators=(",", ":"))
    if gzip_enabled:
        gz_path = Path(f"{path}.gz")
        with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
            handle.write(text)
        return gz_path
    path.write_text(text, encoding="utf-8")
    return path


def write_binary(path: Path, payload: bytes, gzip_enabled: bool = False) -> Path:
    if gzip_enabled:
        gz_path = Path(f"{path}.gz")
        with gzip.open(gz_path, "wb") as handle:
            handle.write(payload)
        return gz_path
    path.write_bytes(payload)
    return path


def estimate_rows_per_chunk(df: pd.DataFrame, max_asset_bytes: int, fmt: str) -> int:
    if len(df) == 0:
        return 1

    probe = df.head(min(1000, len(df)))
    if fmt == "csv":
        encoded = probe.to_csv(index=False).encode("utf-8")
    elif fmt == "json":
        encoded = probe.to_json(orient="records").encode("utf-8")
    else:
        raise ValueError(f"Unsupported chunking format: {fmt}")

    avg_row_bytes = max(1, int(len(encoded) / max(1, len(probe))))
    rows = max(1, int(max_asset_bytes / avg_row_bytes))
    return rows


def write_dataframe_chunks(
    df: pd.DataFrame,
    outdir: Path,
    stem: str,
    max_asset_bytes: int,
    gzip_enabled: bool,
    fmt: str,
) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    rows_per_chunk = estimate_rows_per_chunk(df, max_asset_bytes=max_asset_bytes, fmt=fmt)

    paths: list[Path] = []
    chunk_id = 0
    for start in range(0, len(df), rows_per_chunk):
        end = min(len(df), start + rows_per_chunk)
        chunk = df.iloc[start:end]

        if fmt == "csv":
            ext = "csv"
            path = outdir / f"{stem}_{chunk_id:03d}.{ext}"
            if gzip_enabled:
                path = Path(f"{path}.gz")
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    chunk.to_csv(handle, index=False)
            else:
                chunk.to_csv(path, index=False)
        elif fmt == "json":
            ext = "json"
            path = outdir / f"{stem}_{chunk_id:03d}.{ext}"
            payload = chunk.to_dict(orient="records")
            if gzip_enabled:
                path = Path(f"{path}.gz")
                with gzip.open(path, "wt", encoding="utf-8") as handle:
                    json.dump(payload, handle, separators=(",", ":"))
            else:
                path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        else:
            raise ValueError(f"Unsupported output format: {fmt}")

        paths.append(path)
        chunk_id += 1

    if not paths:
        empty = outdir / f"{stem}_000.{fmt}"
        if fmt == "csv":
            if gzip_enabled:
                empty = Path(f"{empty}.gz")
                with gzip.open(empty, "wt", encoding="utf-8") as handle:
                    df.to_csv(handle, index=False)
            else:
                df.to_csv(empty, index=False)
        else:
            if gzip_enabled:
                empty = Path(f"{empty}.gz")
                with gzip.open(empty, "wt", encoding="utf-8") as handle:
                    json.dump([], handle)
            else:
                empty.write_text("[]", encoding="utf-8")
        paths.append(empty)

    return paths
