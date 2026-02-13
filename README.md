# KaroSpaceBuilder

![KaroSpaceBuilder logo](assets/logo_KSB.png)

Build fully static KaroSpace viewer exports from AnnData (`.h5ad`) with either:

- a modern tabbed GUI (`KaroSpaceBuilder`)
- the CLI (`karospace-export`)

The export output is static files (`index.html`, `manifest.json`, `assets/`) and can be hosted without a Python backend.

## Install

```bash
cd /Users/christoffer/work/karolinska/development/karospace-builder
pip install -e .
```

If your Python is missing Tkinter, create a GUI-capable environment first:

```bash
conda create -n ks-gui -c conda-forge python=3.12 pip tk python.app anndata numpy pandas scipy pillow
conda activate ks-gui
pip install -e /Users/christoffer/work/karolinska/development/karospace-builder
```

## Run KaroSpaceBuilder (GUI)

Preferred command:

```bash
karospace-builder
```

Also supported:

```bash
karospace-export-app
```

Fallback:

```bash
python -m karospace_export.app
```

## GUI Overview

`KaroSpaceBuilder` has four tabs plus runtime logs:

1. `Basic`
2. `Colors & Genes`
3. `Advanced` (expand/collapse)
4. `Help`

### Key GUI features

- Scrollable app layout (works in smaller/non-fullscreen windows)
- `Inspect H5AD` reads the input file and populates searchable dropdowns from:
  - `adata.obs.columns` for `additional_colors` and `groupby`
  - `adata.var_names` for manual genes
- Searchable list builders with:
  - `+ Add`
  - `Remove`
  - `Clear`
- Genes modes:
  - `hvgs` + count
  - `top_mean` + count
  - `list_file` (one gene per line in a text file)
  - `manual_list` (pick genes directly in GUI)
- Preset profiles:
  - `Default`
  - `Pancreas`
  - `Lightweight`

### Inputs by tab

`Basic`:

- `Input .h5ad`: path to AnnData file
- `Output directory`: folder to write static export
- `Coordinates`:
  - `auto`
  - `obsm:spatial`
  - `obs:centroid_x_y`
- `Optional tissue image`: optional override image path
- `Downsample`: optional integer number of cells

`Colors & Genes`:

- `additional_colors (obs columns)`: categorical color fields in viewer
- `groupby list (obs columns)`: extra grouping fields (merged into export annotations)
- `genes mode`:
  - `hvgs` / `top_mean`: requires `N genes`
  - `list_file`: choose gene file path
  - `manual_list`: build genes directly from searchable var names list

`Advanced`:

- `Asset split limit` (`max_asset_mb`)
- `Gzip assets`
- `Write preview.png`
- `Serve after export`
- `Serve port`

## Quick CLI Example

```bash
karospace-export \
  --h5ad /absolute/path/to/input.h5ad \
  --outdir /absolute/path/to/export \
  --anno cell_type \
  --anno leiden \
  --genes hvgs:500 \
  --gzip true \
  --max-asset-mb 16
```

Local preview server:

```bash
karospace-export \
  --h5ad /absolute/path/to/input.h5ad \
  --outdir /absolute/path/to/export \
  --genes top_mean:300 \
  --serve \
  --port 8000
```

## Export Output

Typical output:

```text
export_dir/
  index.html
  manifest.json
  preview.png                  # optional
  assets/
    obs_000.csv.gz
    obs_records_000.json.gz
    var_000.csv.gz
    coords_000.json.gz
    image.png                  # optional
    expression/
      gene_00000_part_000.bin.gz
      gene_00001_part_000.bin.gz
      ...
```

## Input Assumptions

Expected in `.h5ad`:

- `adata.obs` with metadata columns
- coordinates from `adata.obsm["spatial"]` or `adata.obs[["centroid_x", "centroid_y"]]`
- optional image from `adata.uns["spatial"]` or explicit image path

## Troubleshooting

`No module named '_tkinter'`:

- Use a conda env with `tk` (see install section).

`karospace-builder` not found:

```bash
pip install -e /Users/christoffer/work/karolinska/development/karospace-builder
```

Then run:

```bash
karospace-builder
```

## Development

Run tests:

```bash
python -m pytest -q
```
