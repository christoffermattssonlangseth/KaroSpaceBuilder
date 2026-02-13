from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def build_synthetic_adata(
    n_cells: int = 600,
    n_genes: int = 120,
    seed: int = 123,
) -> ad.AnnData:
    rng = np.random.default_rng(seed)

    groups = np.array(["tumor", "stroma", "immune"])
    cell_type = rng.choice(groups, size=n_cells, p=[0.5, 0.3, 0.2])
    leiden = rng.integers(0, 6, size=n_cells).astype(str)

    x = rng.normal(loc=0.0, scale=1.0, size=n_cells)
    y = rng.normal(loc=0.0, scale=1.0, size=n_cells)

    base = rng.gamma(shape=2.0, scale=1.0, size=(n_cells, n_genes)).astype(np.float32)

    marker_0 = np.where(cell_type == "tumor", 3.0, 0.1)
    marker_1 = np.where(cell_type == "immune", 2.5, 0.1)
    marker_2 = np.where(cell_type == "stroma", 2.0, 0.1)

    base[:, 0] += marker_0
    base[:, 1] += marker_1
    base[:, 2] += marker_2

    x_sparse = sparse.csr_matrix(base)

    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(cell_type),
            "leiden": pd.Categorical(leiden),
            "sample_id": "toy_sample",
            "centroid_x": x,
            "centroid_y": y,
        },
        index=[f"cell_{i:05d}" for i in range(n_cells)],
    )

    var = pd.DataFrame(index=[f"Gene{i:04d}" for i in range(n_genes)])
    var["highly_variable"] = False
    var.loc[var.index[: min(30, n_genes)], "highly_variable"] = True
    var["highly_variable_rank"] = np.arange(n_genes)

    adata = ad.AnnData(X=x_sparse, obs=obs, var=var)
    adata.obsm["spatial"] = np.column_stack([x, y]).astype(np.float32)
    return adata


def main() -> None:
    out = Path("synthetic_tiny.h5ad")
    adata = build_synthetic_adata()
    adata.write_h5ad(out)
    print(f"Wrote {out.resolve()}")


if __name__ == "__main__":
    main()
