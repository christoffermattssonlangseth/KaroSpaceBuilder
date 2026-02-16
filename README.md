# KaroSpaceBuilder

![KaroSpaceBuilder logo](assets/logo_KSB.png)

Export an AnnData `.h5ad` into a fully static KaroSpace-compatible viewer bundle.

The end product is static HTML viewer output (`index.html`).
No backend service is generated or required.

## Installation

```bash
git clone https://github.com/christoffermattssonlangseth/KaroSpaceBuilder
cd karospace-builder
python -m pip install -e .
```

If your Python does not include Tkinter (`_tkinter` error), create a GUI-ready environment first:

```bash
conda create -n ks-gui -c conda-forge python=3.12 pip tk python.app anndata numpy pandas scipy pillow
conda activate ks-gui
python -m pip install -e .
```

## Desktop Distribution (Recommended for Non-technical Users)

End users do not need Python or terminal access.

1. Go to the GitHub Releases page.
2. Download the file for your OS:
   - `KaroSpaceBuilder-macos.zip`
   - `KaroSpaceBuilder-windows.zip`
3. Unzip and open `KaroSpaceBuilder`.
4. Export your dataset; the output is `index.html`.

## Build Desktop Releases (Maintainers)

This repository includes automated packaging with PyInstaller:

- Workflow: `.github/workflows/desktop-release.yml`
- Trigger:
  - `workflow_dispatch` (manual)
  - `release.published` (auto-attach binaries to the release)

Local build command:

```bash
python -m pip install -e .[build]
python scripts/build_desktop.py
```

Build output:

- macOS: `dist/KaroSpaceBuilder.app`
- Windows: `dist/KaroSpaceBuilder/KaroSpaceBuilder.exe`

## What You Get

Default output (single-file mode):

```text
export_dir/
  index.html
```

That `index.html` is the final viewer artifact.

If you explicitly want split files, use CLI flag `--multi-file` to also keep `manifest.json` and `assets/`.

## Quick Start (GUI App)

### 1) Create a GUI-capable environment

If your current Python lacks `_tkinter`, use this:

```bash
conda create -n ks-gui -c conda-forge python=3.12 pip tk python.app anndata numpy pandas scipy pillow
conda activate ks-gui
```

### 2) Install this project

```bash
python -m pip install -e .
```

### 3) Launch app

```bash
KaroSpaceBuilder
```

If the command is not found, run:

```bash
karospace-builder
```

or:

```bash
karospace-export-app
```

or:

```bash
python -m karospace_export.app
```

After export finishes, the builder writes `index.html` in your output folder.

GUI features:

- File pickers for input/output/image/gene-list paths.
- Controls for coordinate source, annotations, gene mode, downsampling, gzip, and asset chunk size.
- Background export with logs.
- Optional local preview server.

## Quick Start (CLI)

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

Optional local serving:

```bash
karospace-export \
  --h5ad /absolute/path/to/input.h5ad \
  --outdir /absolute/path/to/export \
  --anno cell_type \
  --genes top_mean:300 \
  --serve \
  --port 8000
```

## CLI Arguments

- `--h5ad PATH` required.
- `--outdir PATH` required.
- `--coords {obsm:spatial,obs:centroid_x_y}` optional (auto-detect if omitted).
- `--anno COLNAME` repeatable.
- `--genes MODE` where mode is `hvgs:N`, `top_mean:N`, or `list:PATH`.
- `--image PATH` optional.
- `--downsample N` optional.
- `--gzip true/false` default `true`.
- `--max-asset-mb NUMBER` default `16`.
- `--serve` optional preview server.
- `--port` optional server port (default `8000`).
- `--no-preview` disable `preview.png`.
- `--multi-file` keep `manifest.json` + `assets/` files instead of embedding everything into `index.html`.

## Input Assumptions

Expected in `.h5ad`:

- `adata.obs` with annotation columns.
- Coordinates from `adata.obsm["spatial"]` or `obs` columns `centroid_x` and `centroid_y`.
- Optional tissue image from `adata.uns["spatial"]` or `--image PATH`.

## Export Format (Multi-file Mode)

When you use `--multi-file`, `manifest.json` includes:

- `dataset_name`
- `n_cells`
- `n_genes_exported`
- `coordinate_source`
- `annotation_columns`
- `gene_list`
- `asset_files`
- `expression`
- `schema_version`
- optional `image`
- optional `preview_path`

Expression is stored as per-gene `float32` vectors in chunked binary files (`per_gene_float32_v1`) so the browser can lazy-load one gene at a time.
In default single-file mode, these files are embedded into `index.html`.

## Viewer Behavior

`index.html` is dependency-light and supports:

- Fast point rendering in canvas.
- Color by annotation or selected gene expression.
- Hover tooltip with `obs` fields.
- Category filter.
- Point size and opacity controls.

In default single-file mode, no folder prompt is needed because data is embedded in the HTML.

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

## Example Dataset

Generate a tiny synthetic `.h5ad`:

```bash
python examples/generate_synthetic_h5ad.py
```

## Troubleshooting

### `No module named '_tkinter'`

Your Python build does not include Tkinter. Use a dedicated Conda environment with `python=3.12` and `tk`:

```bash
conda create -n ks-gui -c conda-forge python=3.12 pip tk python.app anndata numpy pandas scipy pillow
conda activate ks-gui
python -m pip install -e .
KaroSpaceBuilder
```

### `KaroSpaceBuilder` command not found

Reinstall in the active environment:

```bash
python -m pip install -e .
```

Then run either:

```bash
KaroSpaceBuilder
```

or:

```bash
karospace-builder
```

or:

```bash
karospace-export-app
```

or:

```bash
python -m karospace_export.app
```

## Development

Run tests:

```bash
python -m pytest -q
```
