from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(slots=True)
class GeneMode:
    kind: str
    value: int | Path


def parse_genes_mode(raw: str) -> GeneMode:
    if ":" not in raw:
        raise ValueError(
            "Invalid --genes value. Expected one of hvgs:N, top_mean:N, list:PATH"
        )

    kind, value = raw.split(":", 1)
    kind = kind.strip().lower()
    value = value.strip()

    if kind in {"hvgs", "top_mean"}:
        try:
            count = int(value)
        except ValueError as exc:
            raise ValueError(f"Gene mode {kind} expects an integer count") from exc
        if count <= 0:
            raise ValueError(f"Gene count for {kind} must be > 0")
        return GeneMode(kind=kind, value=count)

    if kind == "list":
        path = Path(value)
        if not path.exists():
            raise ValueError(f"Gene list file not found: {path}")
        return GeneMode(kind=kind, value=path)

    raise ValueError("Invalid --genes mode. Use hvgs:N, top_mean:N, or list:PATH")


def _compute_means_by_blocks(
    adata,
    obs_indices: np.ndarray,
    block_size: int = 512,
) -> np.ndarray:
    means = np.zeros(adata.n_vars, dtype=np.float64)
    full = len(obs_indices) == adata.n_obs and np.array_equal(obs_indices, np.arange(adata.n_obs))
    selector = slice(None) if full else obs_indices

    for start in range(0, adata.n_vars, block_size):
        end = min(adata.n_vars, start + block_size)
        x_block = adata[selector, start:end].X

        if sparse.issparse(x_block):
            block_means = np.asarray(x_block.mean(axis=0)).ravel()
        else:
            x_arr = np.asarray(x_block)
            block_means = x_arr.mean(axis=0)

        means[start:end] = block_means

    return means


def _pick_top_mean_genes(adata, obs_indices: np.ndarray, count: int) -> list[str]:
    means = _compute_means_by_blocks(adata, obs_indices)
    ranked_idx = np.argsort(means)[::-1]
    top_idx = ranked_idx[:count]
    return [str(adata.var_names[i]) for i in top_idx]


def _pick_hvgs(adata, count: int, obs_indices: np.ndarray) -> list[str]:
    var: pd.DataFrame = adata.var

    if "highly_variable" in var.columns:
        hv = var[var["highly_variable"].astype(bool)]

        if not hv.empty:
            if "highly_variable_rank" in hv.columns:
                hv = hv.sort_values("highly_variable_rank", ascending=True)
            elif "dispersions_norm" in hv.columns:
                hv = hv.sort_values("dispersions_norm", ascending=False)

            selected = [str(name) for name in hv.index[:count].tolist()]
            if len(selected) >= count:
                return selected

            extra = _pick_top_mean_genes(adata, obs_indices=obs_indices, count=count)
            seen = set(selected)
            for gene in extra:
                if gene not in seen:
                    selected.append(gene)
                    seen.add(gene)
                if len(selected) == count:
                    break
            return selected

    return _pick_top_mean_genes(adata, obs_indices=obs_indices, count=count)


def resolve_gene_names(adata, mode: GeneMode, obs_indices: np.ndarray) -> list[str]:
    var_names = pd.Index([str(x) for x in adata.var_names])

    if mode.kind == "top_mean":
        genes = _pick_top_mean_genes(adata, obs_indices=obs_indices, count=int(mode.value))
    elif mode.kind == "hvgs":
        genes = _pick_hvgs(adata, count=int(mode.value), obs_indices=obs_indices)
    elif mode.kind == "list":
        gene_lines = Path(mode.value).read_text(encoding="utf-8").splitlines()
        wanted = [line.strip() for line in gene_lines if line.strip()]
        missing = [gene for gene in wanted if gene not in var_names]
        if missing:
            preview = ", ".join(missing[:10])
            raise ValueError(f"{len(missing)} genes from list are missing in var_names: {preview}")
        genes = wanted
    else:
        raise ValueError(f"Unsupported gene mode: {mode.kind}")

    deduped: list[str] = []
    seen: set[str] = set()
    for gene in genes:
        if gene not in seen:
            deduped.append(gene)
            seen.add(gene)

    return deduped
