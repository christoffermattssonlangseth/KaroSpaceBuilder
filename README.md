# karospace-export

A small Python app + CLI that exports an AnnData `.h5ad` file into a fully static KaroSpace-compatible HTML viewer bundle.

## Proposed Export Schema

Output layout:

```text
export_dir/
  index.html
  manifest.json
  preview.png                  # optional
  assets/
    obs_000.csv.gz             # metadata table (chunked)
    obs_records_000.json.gz    # same metadata for viewer fast-load
    var_000.csv.gz             # exported genes table
    coords_000.json.gz         # coordinates table (chunked)
    image.png                  # optional tissue image
    expression/
      gene_00000_part_000.bin.gz
      gene_00001_part_000.bin.gz
      ...
```

`manifest.json` contains:

- `dataset_name`
- `n_cells`
- `n_genes_exported`
- `coordinate_source`
- `annotation_columns`
- `gene_list`
- `asset_files` (path + byte size + kind)
- `expression` spec (`per_gene_float32_v1`)
- `schema_version`
- optional `image`, `preview_path`

## Why per-gene vectors (format B)

Expression is exported as per-gene `float32` vectors (`gene_XXXXX_part_YYY.bin[.gz]`) instead of a full matrix blob. This allows lazy loading in the browser: the viewer loads only the selected gene on demand, which is memory-friendly for large datasets.

## Installation

```bash
pip install -e .
```

## GUI App

Launch the desktop app:

```bash
karospace-export-app
```

What the GUI provides:

- file pickers for input/output/image/gene-list paths
- simple controls for coordinates, annotations, genes mode, downsampling, gzip, and chunk size
- background export with runtime logs
- optional one-click local preview server

## CLI Usage

```bash
karospace-export \
  --h5ad /path/to/input.h5ad \
  --outdir /path/to/export \
  --anno cell_type \
  --anno leiden \
  --genes hvgs:500 \
  --gzip true \
  --max-asset-mb 16
```

Arguments:

- `--h5ad PATH` (required)
- `--outdir PATH` (required)
- `--coords {obsm:spatial,obs:centroid_x_y}` (optional, auto-detected otherwise)
- `--anno COLNAME` (repeatable)
- `--genes MODE`
  - `hvgs:N`
  - `top_mean:N`
  - `list:PATH`
- `--image PATH` (optional)
- `--downsample N` (optional)
- `--gzip true/false` (default: `true`)
- `--max-asset-mb NUMBER` (default: `16`)
- `--serve` (optional local static server)
- `--port` (default `8000`)

## Python API

```python
from pathlib import Path
from karospace_export import ExportConfig, export_h5ad

manifest = export_h5ad(
    ExportConfig(
        h5ad_path=Path("input.h5ad"),
        outdir=Path("export"),
        genes_mode="top_mean:300",
        annotation_columns=["cell_type", "leiden"],
    )
)
```

## Static Viewer Notes

- `index.html` loads all assets with relative paths.
- On static hosting (Cloudflare Pages/R2), it works directly.
- On `file://` browsers that block local fetch, the viewer prompts once for the export folder and reads assets via the browser File API (still no server required).

## Example Dataset

Generate a tiny synthetic `.h5ad`:

```bash
python examples/generate_synthetic_h5ad.py
```

## Test

```bash
pytest
```
