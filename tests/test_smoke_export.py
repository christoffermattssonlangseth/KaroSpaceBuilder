from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from karospace_export.export import ExportConfig, export_h5ad


def _make_test_adata() -> ad.AnnData:
    rng = np.random.default_rng(7)
    n_cells = 80
    n_genes = 24

    x = sparse.csr_matrix(rng.gamma(shape=1.8, scale=1.2, size=(n_cells, n_genes)).astype(np.float32))
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(rng.choice(["A", "B", "C"], size=n_cells)),
            "leiden": pd.Categorical(rng.integers(0, 4, size=n_cells).astype(str)),
        },
        index=[f"cell_{i:04d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"Gene{i:03d}" for i in range(n_genes)])
    var["highly_variable"] = False
    var.loc[var.index[:8], "highly_variable"] = True
    var["highly_variable_rank"] = np.arange(n_genes)

    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obsm["spatial"] = np.column_stack(
        [rng.normal(0, 1, size=n_cells), rng.normal(0, 1, size=n_cells)]
    ).astype(np.float32)
    return adata


def test_smoke_export(tmp_path):
    h5ad_path = tmp_path / "toy.h5ad"
    outdir = tmp_path / "export"

    _make_test_adata().write_h5ad(h5ad_path)

    config = ExportConfig(
        h5ad_path=h5ad_path,
        outdir=outdir,
        annotation_columns=["cell_type", "leiden"],
        genes_mode="hvgs:5",
        gzip=False,
        max_asset_mb=2,
    )

    manifest = export_h5ad(config)

    assert manifest.n_cells == 80
    assert manifest.n_genes_exported == 5

    assert (outdir / "index.html").exists()
    assert (outdir / "manifest.json").exists()

    manifest_json = json.loads((outdir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest_json["schema_version"] == "1.0.0"
    assert manifest_json["expression"]["format"] == "per_gene_float32_v1"

    for asset in manifest_json["asset_files"]:
        assert (outdir / asset["path"]).exists(), asset["path"]

    html = (outdir / "index.html").read_text(encoding="utf-8")
    assert "manifest.json" in html
    assert "assets" in html
