from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any
import shutil

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .schema import (
    AssetFile,
    ExpressionChunk,
    ExpressionGene,
    ExpressionSpec,
    ImageSpec,
    Manifest,
    write_manifest,
)
from .utils import file_size, parse_genes_mode, resolve_gene_names, write_binary, write_dataframe_chunks


@dataclass(slots=True)
class ExportConfig:
    h5ad_path: Path
    outdir: Path
    coords: str | None = None
    annotation_columns: list[str] | None = None
    genes_mode: str = "hvgs:500"
    image_path: Path | None = None
    downsample: int | None = None
    gzip: bool = True
    max_asset_mb: float = 16.0
    dataset_name: str | None = None
    preview: bool = True


def _max_asset_bytes(max_asset_mb: float) -> int:
    if max_asset_mb <= 0:
        raise ValueError("--max-asset-mb must be > 0")
    return int(max_asset_mb * 1024 * 1024)


def _read_adata(path: Path):
    try:
        return ad.read_h5ad(path, backed="r")
    except Exception:
        return ad.read_h5ad(path)


def _close_adata(adata: Any) -> None:
    if getattr(adata, "isbacked", False):
        file_obj = getattr(adata, "file", None)
        if file_obj is not None:
            file_obj.close()


def _resolve_obs_indices(n_obs: int, downsample: int | None) -> np.ndarray:
    indices = np.arange(n_obs)
    if downsample is None or downsample >= n_obs:
        return indices

    if downsample <= 0:
        raise ValueError("--downsample must be > 0")

    rng = np.random.default_rng(42)
    chosen = np.sort(rng.choice(indices, size=downsample, replace=False))
    return chosen


def _resolve_coordinates(adata, obs_idx: np.ndarray, requested: str | None) -> tuple[np.ndarray, str]:
    obs = adata.obs

    def use_obsm_spatial() -> tuple[np.ndarray, str]:
        spatial = adata.obsm["spatial"]
        arr = np.asarray(spatial)
        if arr.ndim != 2 or arr.shape[1] < 2:
            raise ValueError("obsm['spatial'] must have shape (n_cells, 2+)")
        return arr[obs_idx, :2].astype(np.float32), "obsm:spatial"

    def use_obs_centroid() -> tuple[np.ndarray, str]:
        for col in ("centroid_x", "centroid_y"):
            if col not in obs.columns:
                raise ValueError(
                    "--coords obs:centroid_x_y requested but centroid_x/centroid_y are missing in obs"
                )
        coords = obs.iloc[obs_idx][["centroid_x", "centroid_y"]].to_numpy(dtype=np.float32)
        return coords, "obs:centroid_x_y"

    if requested is not None:
        if requested == "obsm:spatial":
            if "spatial" not in adata.obsm:
                raise ValueError("--coords obsm:spatial requested but obsm['spatial'] is missing")
            return use_obsm_spatial()
        if requested == "obs:centroid_x_y":
            return use_obs_centroid()
        raise ValueError("--coords must be one of obsm:spatial, obs:centroid_x_y")

    if "spatial" in adata.obsm:
        return use_obsm_spatial()

    if {"centroid_x", "centroid_y"}.issubset(set(obs.columns)):
        return use_obs_centroid()

    raise ValueError(
        "Could not auto-detect coordinates. Provide --coords and ensure either obsm['spatial'] or obs centroid_x/centroid_y exists."
    )


def _resolve_annotation_columns(obs: pd.DataFrame, requested: list[str] | None) -> list[str]:
    if requested:
        missing = [name for name in requested if name not in obs.columns]
        if missing:
            raise ValueError(f"Missing annotation columns in obs: {', '.join(missing)}")
        return requested

    auto_cols: list[str] = []
    for col in obs.columns:
        series = obs[col]
        if str(series.dtype) in {"object", "category", "bool"}:
            auto_cols.append(col)
        if len(auto_cols) >= 4:
            break

    if not auto_cols:
        auto_cols = [str(c) for c in obs.columns[: min(4, len(obs.columns))].tolist()]

    if not auto_cols:
        raise ValueError("No annotation columns available in obs. Pass --anno with valid columns.")

    return auto_cols


def _extract_gene_vector(adata, obs_idx: np.ndarray, gene_idx: int) -> np.ndarray:
    full = len(obs_idx) == adata.n_obs and np.array_equal(obs_idx, np.arange(adata.n_obs))
    selector = slice(None) if full else obs_idx
    x = adata[selector, gene_idx].X

    if sparse.issparse(x):
        arr = np.asarray(x.toarray()).ravel()
    else:
        arr = np.asarray(x).reshape(-1)

    arr = arr.astype(np.float32, copy=False)
    if np.isnan(arr).any():
        arr = np.nan_to_num(arr, copy=False)
    return arr


def _build_obs_export(obs: pd.DataFrame, obs_idx: np.ndarray, anno_cols: list[str]) -> pd.DataFrame:
    subset = obs.iloc[obs_idx].copy()
    out = pd.DataFrame({"cell_id": subset.index.astype(str)})

    for col in anno_cols:
        series = subset[col]
        if str(series.dtype) == "category":
            out[col] = series.astype(str)
        elif pd.api.types.is_numeric_dtype(series):
            out[col] = series.to_numpy()
        else:
            out[col] = series.astype(str)

    return out


def _build_var_export(adata, selected_genes: list[str]) -> pd.DataFrame:
    positions = {str(name): i for i, name in enumerate(adata.var_names)}
    rows = []
    for gene in selected_genes:
        pos = positions[gene]
        rows.append({"gene": gene, "source_index": int(pos)})
    return pd.DataFrame(rows)


def _write_expression_assets(
    adata,
    obs_idx: np.ndarray,
    selected_genes: list[str],
    outdir: Path,
    expr_dir: Path,
    gzip_assets: bool,
    max_asset_bytes: int,
) -> tuple[ExpressionSpec, list[AssetFile]]:
    expr_dir.mkdir(parents=True, exist_ok=True)
    mapping = {str(name): i for i, name in enumerate(adata.var_names)}

    values_per_chunk = max(1, max_asset_bytes // np.dtype("<f4").itemsize)

    asset_files: list[AssetFile] = []
    genes_out: list[ExpressionGene] = []

    for gene_index, gene in enumerate(selected_genes):
        source_idx = mapping[gene]
        vec = _extract_gene_vector(adata, obs_idx=obs_idx, gene_idx=source_idx)

        chunks: list[ExpressionChunk] = []
        part = 0
        for start in range(0, len(vec), values_per_chunk):
            end = min(len(vec), start + values_per_chunk)
            payload = vec[start:end].astype("<f4", copy=False).tobytes(order="C")

            base_path = expr_dir / f"gene_{gene_index:05d}_part_{part:03d}.bin"
            written = write_binary(base_path, payload, gzip_enabled=gzip_assets)
            rel_path = written.relative_to(outdir).as_posix()
            size = file_size(written)

            chunks.append(ExpressionChunk(path=rel_path, start=start, end=end, size_bytes=size))
            asset_files.append(AssetFile(kind="expression_chunk", path=rel_path, size_bytes=size))
            part += 1

        genes_out.append(ExpressionGene(gene=gene, index=gene_index, chunks=chunks))

    spec = ExpressionSpec(
        format="per_gene_float32_v1",
        dtype="float32",
        n_cells=len(obs_idx),
        gzip=gzip_assets,
        genes=genes_out,
    )
    return spec, asset_files


def _export_image(
    adata,
    outdir: Path,
    assets_dir: Path,
    image_path: Path | None,
) -> tuple[ImageSpec | None, list[AssetFile]]:
    asset_files: list[AssetFile] = []

    def _image_spec(path: Path) -> ImageSpec | None:
        try:
            from PIL import Image

            with Image.open(path) as img:
                width, height = img.size
        except Exception:
            return None
        return ImageSpec(path=path.relative_to(outdir).as_posix(), width=width, height=height)

    if image_path is not None:
        if not image_path.exists():
            raise ValueError(f"--image not found: {image_path}")
        target = assets_dir / f"image{image_path.suffix.lower() or '.png'}"
        shutil.copy2(image_path, target)
        spec = _image_spec(target)
        asset_files.append(
            AssetFile(kind="image", path=target.relative_to(outdir).as_posix(), size_bytes=file_size(target))
        )
        return spec, asset_files

    spatial = adata.uns.get("spatial") if hasattr(adata, "uns") else None
    if not isinstance(spatial, dict):
        return None, asset_files

    try:
        from PIL import Image
    except Exception:
        return None, asset_files

    for _, payload in spatial.items():
        if not isinstance(payload, dict):
            continue
        images = payload.get("images")
        if not isinstance(images, dict):
            continue

        for key in ("hires", "lowres"):
            img_array = images.get(key)
            if img_array is None:
                continue

            arr = np.asarray(img_array)
            if arr.ndim not in {2, 3}:
                continue
            if arr.dtype != np.uint8:
                arr = np.clip(arr, 0.0, 1.0) * 255.0
                arr = arr.astype(np.uint8)

            pil = Image.fromarray(arr)
            target = assets_dir / "image.png"
            pil.save(target)
            spec = ImageSpec(path=target.relative_to(outdir).as_posix(), width=pil.width, height=pil.height)
            asset_files.append(
                AssetFile(kind="image", path=target.relative_to(outdir).as_posix(), size_bytes=file_size(target))
            )
            return spec, asset_files

    return None, asset_files


def _write_preview(outdir: Path, coords: np.ndarray) -> str | None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return None

    if coords.size == 0:
        return None

    x = coords[:, 0].astype(np.float32)
    y = coords[:, 1].astype(np.float32)

    x_min, x_max = float(np.min(x)), float(np.max(x))
    y_min, y_max = float(np.min(y)), float(np.max(y))

    w, h = 512, 512
    pad = 20
    span_x = max(1e-6, x_max - x_min)
    span_y = max(1e-6, y_max - y_min)

    px = pad + ((x - x_min) / span_x) * (w - 2 * pad)
    py = pad + ((y - y_min) / span_y) * (h - 2 * pad)

    img = Image.new("RGB", (w, h), "white")
    draw = ImageDraw.Draw(img)
    for xi, yi in zip(px, py):
        draw.ellipse((xi - 1, yi - 1, xi + 1, yi + 1), fill=(38, 70, 83))

    preview = outdir / "preview.png"
    img.save(preview)
    return preview.name


def _build_viewer_html() -> str:
    return """<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\" />
  <title>KaroSpace Static Viewer</title>
  <style>
    :root {
      --bg: #f5f1e9;
      --panel: #fffdf8;
      --ink: #182024;
      --muted: #5f6b73;
      --line: #d6d1c5;
      --accent: #2a9d8f;
      --danger: #b22222;
    }
    * { box-sizing: border-box; }
    html, body { margin: 0; height: 100%; font-family: ui-sans-serif, system-ui, -apple-system, Segoe UI, sans-serif; color: var(--ink); background: radial-gradient(circle at 10% 10%, #fffef8 0%, var(--bg) 70%); }
    #app { display: grid; grid-template-columns: 320px 1fr; height: 100vh; gap: 14px; padding: 14px; }
    .panel {
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 14px;
      padding: 12px;
      box-shadow: 0 6px 24px rgba(0,0,0,0.06);
      overflow: auto;
    }
    .control { margin-bottom: 12px; }
    label { display: block; font-size: 12px; text-transform: uppercase; letter-spacing: .04em; margin-bottom: 4px; color: var(--muted); }
    select, input[type=\"range\"], button {
      width: 100%;
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 6px 8px;
      background: #fff;
      color: var(--ink);
    }
    button { cursor: pointer; }
    .status { font-size: 13px; color: var(--muted); min-height: 24px; }
    .status.error { color: var(--danger); }
    .stage {
      position: relative;
      background: var(--panel);
      border: 1px solid var(--line);
      border-radius: 14px;
      overflow: hidden;
      min-height: 420px;
    }
    canvas { width: 100%; height: 100%; display: block; }
    #tooltip {
      position: absolute;
      pointer-events: none;
      background: rgba(24,32,36,.92);
      color: #fff;
      border-radius: 8px;
      padding: 8px;
      font-size: 12px;
      max-width: 260px;
      display: none;
      z-index: 30;
    }
    #legend {
      position: absolute;
      right: 10px;
      bottom: 10px;
      background: rgba(255,255,255,.92);
      border: 1px solid var(--line);
      border-radius: 8px;
      padding: 8px;
      max-width: 280px;
      font-size: 12px;
      color: var(--ink);
      z-index: 20;
    }
    .legend-item { display: flex; align-items: center; gap: 6px; margin: 2px 0; }
    .swatch { width: 12px; height: 12px; border-radius: 3px; border: 1px solid rgba(0,0,0,.15); }
    #localPrompt {
      position: absolute;
      inset: 0;
      background: rgba(255,255,255,.94);
      display: none;
      align-items: center;
      justify-content: center;
      text-align: center;
      z-index: 40;
      padding: 24px;
    }
    #localPromptInner { max-width: 420px; }
    @media (max-width: 960px) {
      #app { grid-template-columns: 1fr; grid-template-rows: auto 1fr; }
      .panel { max-height: 42vh; }
    }
  </style>
</head>
<body>
  <div id=\"app\">
    <section class=\"panel\">
      <h2 style=\"margin:0 0 8px 0; font-size:18px;\">KaroSpace Viewer</h2>
      <div id=\"datasetMeta\" style=\"font-size:13px; color:var(--muted); margin-bottom: 14px;\"></div>

      <div class=\"control\">
        <label for=\"annotationSelect\">Annotation Column</label>
        <select id=\"annotationSelect\"></select>
      </div>

      <div class=\"control\">
        <label for=\"filterSelect\">Category Filter</label>
        <select id=\"filterSelect\"></select>
      </div>

      <div class=\"control\">
        <label for=\"geneSelect\">Gene Expression Overlay</label>
        <select id=\"geneSelect\"></select>
      </div>

      <div class=\"control\">
        <label for=\"sizeSlider\">Point Size</label>
        <input id=\"sizeSlider\" type=\"range\" min=\"1\" max=\"8\" step=\"1\" value=\"2\" />
      </div>

      <div class=\"control\">
        <label for=\"opacitySlider\">Point Opacity</label>
        <input id=\"opacitySlider\" type=\"range\" min=\"0.05\" max=\"1\" step=\"0.05\" value=\"0.8\" />
      </div>

      <div id=\"status\" class=\"status\"></div>
    </section>

    <section class=\"stage\" id=\"stage\">
      <canvas id=\"plot\"></canvas>
      <div id=\"tooltip\"></div>
      <div id=\"legend\"></div>
      <div id=\"localPrompt\">
        <div id=\"localPromptInner\">
          <p style=\"font-size:14px;\">Local file mode detected. If browser blocks direct asset reads, pick this export folder once so the viewer can read the files.</p>
          <button id=\"pickFolder\">Pick Export Folder</button>
          <input id=\"folderInput\" type=\"file\" webkitdirectory multiple style=\"display:none\" />
        </div>
      </div>
    </section>
  </div>

  <script>
    const MANIFEST_PATH = "manifest.json";
    const palette = [
      "#264653", "#2a9d8f", "#e9c46a", "#f4a261", "#e76f51", "#457b9d", "#8d99ae", "#ef476f", "#06d6a0", "#118ab2"
    ];

    const state = {
      manifest: null,
      obs: [],
      coords: [],
      pointSize: 2,
      opacity: 0.8,
      annotationCol: null,
      filterValue: "__ALL__",
      activeGene: "__NONE__",
      geneCache: new Map(),
      width: 0,
      height: 0,
      px: null,
      py: null,
      image: null,
      hoverGrid: null,
      currentColorMap: new Map(),
      currentMin: 0,
      currentMax: 1,
    };

    const stage = document.getElementById("stage");
    const canvas = document.getElementById("plot");
    const tooltip = document.getElementById("tooltip");
    const legend = document.getElementById("legend");
    const statusEl = document.getElementById("status");
    const datasetMeta = document.getElementById("datasetMeta");
    const annotationSelect = document.getElementById("annotationSelect");
    const filterSelect = document.getElementById("filterSelect");
    const geneSelect = document.getElementById("geneSelect");
    const sizeSlider = document.getElementById("sizeSlider");
    const opacitySlider = document.getElementById("opacitySlider");
    const localPrompt = document.getElementById("localPrompt");
    const pickFolder = document.getElementById("pickFolder");
    const folderInput = document.getElementById("folderInput");
    const ctx = canvas.getContext("2d", { alpha: false });

    function setStatus(msg, isError = false) {
      statusEl.textContent = msg;
      statusEl.classList.toggle("error", isError);
    }

    const localFiles = {
      map: new Map(),
      ready: false,
      waiting: null,

      async ensureReady() {
        if (this.ready) return;
        if (this.waiting) return this.waiting;

        localPrompt.style.display = "flex";
        this.waiting = new Promise((resolve, reject) => {
          const onPick = () => folderInput.click();
          pickFolder.addEventListener("click", onPick, { once: true });

          folderInput.addEventListener("change", () => {
            const files = Array.from(folderInput.files || []);
            if (!files.length) {
              reject(new Error("No files selected."));
              return;
            }

            this.map.clear();
            for (const file of files) {
              const rel = (file.webkitRelativePath || file.name).replace(/^\\.\\//, "");
              this.map.set(rel, file);

              const parts = rel.split("/");
              if (parts.length > 1) {
                const stripped = parts.slice(1).join("/");
                this.map.set(stripped, file);
              }

              this.map.set(file.name, file);
            }

            this.ready = true;
            localPrompt.style.display = "none";
            resolve();
          }, { once: true });
        });

        return this.waiting;
      },

      get(path) {
        const clean = path.replace(/^\\.\\//, "");
        return this.map.get(clean) || this.map.get(path);
      },

      async readArrayBuffer(path) {
        const file = this.get(path);
        if (!file) {
          throw new Error(`Local file not found: ${path}. Pick the folder containing manifest.json and assets/.`);
        }
        return file.arrayBuffer();
      }
    };

    async function maybeGunzip(buffer, path) {
      if (!path.endsWith(".gz")) {
        return buffer;
      }

      const bytes = new Uint8Array(buffer);
      const isGzipPayload = bytes.length >= 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
      if (!isGzipPayload) {
        return buffer;
      }

      if (!("DecompressionStream" in window)) {
        throw new Error("Gzip assets require browser DecompressionStream support. Re-export with --gzip false if needed.");
      }

      const stream = new Response(buffer).body.pipeThrough(new DecompressionStream("gzip"));
      return new Response(stream).arrayBuffer();
    }

    async function readBinary(path) {
      const normalized = path.replace(/^\\.\\//, "");
      try {
        const resp = await fetch(normalized);
        if (!resp.ok) {
          throw new Error(`HTTP ${resp.status}`);
        }
        const raw = await resp.arrayBuffer();
        return maybeGunzip(raw, normalized);
      } catch (err) {
        if (location.protocol !== "file:") {
          throw err;
        }

        await localFiles.ensureReady();
        const raw = await localFiles.readArrayBuffer(normalized);
        return maybeGunzip(raw, normalized);
      }
    }

    async function readText(path) {
      const buffer = await readBinary(path);
      return new TextDecoder().decode(buffer);
    }

    async function readJSON(path) {
      const text = await readText(path);
      return JSON.parse(text);
    }

    function relPaths(kind) {
      return state.manifest.asset_files.filter((x) => x.kind === kind).map((x) => x.path);
    }

    async function loadObs() {
      const paths = relPaths("obs_json");
      const rows = [];
      for (const path of paths) {
        const chunk = await readJSON(path);
        rows.push(...chunk);
      }
      return rows;
    }

    async function loadCoords() {
      const paths = relPaths("coords_json");
      const rows = [];
      for (const path of paths) {
        const chunk = await readJSON(path);
        rows.push(...chunk);
      }
      return rows;
    }

    async function loadImage() {
      if (!state.manifest.image || !state.manifest.image.path) {
        return null;
      }
      const raw = await readBinary(state.manifest.image.path);
      const blob = new Blob([raw]);
      const url = URL.createObjectURL(blob);
      const img = new Image();
      await new Promise((resolve, reject) => {
        img.onload = () => resolve();
        img.onerror = () => reject(new Error("Unable to load image asset."));
        img.src = url;
      });
      return img;
    }

    function resizeCanvas() {
      const dpr = window.devicePixelRatio || 1;
      const w = stage.clientWidth;
      const h = stage.clientHeight;
      canvas.width = Math.floor(w * dpr);
      canvas.height = Math.floor(h * dpr);
      canvas.style.width = `${w}px`;
      canvas.style.height = `${h}px`;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      state.width = w;
      state.height = h;
      projectCoords();
      rebuildHoverGrid();
      draw();
    }

    function projectCoords() {
      const n = state.coords.length;
      if (!n) return;

      let xMin = Number.POSITIVE_INFINITY;
      let xMax = Number.NEGATIVE_INFINITY;
      let yMin = Number.POSITIVE_INFINITY;
      let yMax = Number.NEGATIVE_INFINITY;
      const xs = new Float32Array(n);
      const ys = new Float32Array(n);

      for (let i = 0; i < n; i += 1) {
        const x = Number(state.coords[i].x);
        const y = Number(state.coords[i].y);
        xs[i] = x;
        ys[i] = y;
        if (x < xMin) xMin = x;
        if (x > xMax) xMax = x;
        if (y < yMin) yMin = y;
        if (y > yMax) yMax = y;
      }

      const pad = 22;
      const spanX = Math.max(1e-6, xMax - xMin);
      const spanY = Math.max(1e-6, yMax - yMin);

      state.px = new Float32Array(n);
      state.py = new Float32Array(n);
      for (let i = 0; i < n; i += 1) {
        state.px[i] = pad + ((xs[i] - xMin) / spanX) * (state.width - 2 * pad);
        state.py[i] = pad + ((ys[i] - yMin) / spanY) * (state.height - 2 * pad);
      }
    }

    function rebuildHoverGrid() {
      const n = state.coords.length;
      if (!n || !state.px || !state.py) {
        state.hoverGrid = null;
        return;
      }

      const cell = Math.max(6, state.pointSize * 4);
      const map = new Map();

      for (let i = 0; i < n; i += 1) {
        const gx = Math.floor(state.px[i] / cell);
        const gy = Math.floor(state.py[i] / cell);
        const key = `${gx}:${gy}`;
        const bucket = map.get(key);
        if (bucket) {
          bucket.push(i);
        } else {
          map.set(key, [i]);
        }
      }

      state.hoverGrid = { cell, map };
    }

    function annotationValue(i) {
      const key = state.annotationCol;
      return key ? String(state.obs[i][key]) : "";
    }

    function isVisible(i) {
      if (state.filterValue === "__ALL__") return true;
      return annotationValue(i) === state.filterValue;
    }

    function lerp(a, b, t) {
      return Math.round(a + (b - a) * t);
    }

    function geneColor(v, lo, hi) {
      const t = Math.max(0, Math.min(1, (v - lo) / Math.max(1e-9, hi - lo)));
      const r = lerp(20, 239, t);
      const g = lerp(44, 138, t);
      const b = lerp(89, 98, t);
      return `rgb(${r},${g},${b})`;
    }

    function draw() {
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(0, 0, state.width, state.height);

      if (state.image) {
        ctx.globalAlpha = 0.32;
        ctx.drawImage(state.image, 0, 0, state.width, state.height);
      }

      ctx.globalAlpha = state.opacity;
      const n = state.coords.length;
      const size = state.pointSize;
      const hasGene = state.activeGene !== "__NONE__";
      const geneVec = hasGene ? state.geneCache.get(state.activeGene) : null;

      for (let i = 0; i < n; i += 1) {
        if (!isVisible(i)) continue;

        let color = "#666";
        if (hasGene && geneVec) {
          color = geneColor(geneVec[i], state.currentMin, state.currentMax);
        } else {
          const category = annotationValue(i);
          color = state.currentColorMap.get(category) || "#888";
        }

        ctx.fillStyle = color;
        const x = state.px[i];
        const y = state.py[i];
        if (size <= 2) {
          ctx.fillRect(x, y, size, size);
        } else {
          ctx.beginPath();
          ctx.arc(x, y, size, 0, Math.PI * 2);
          ctx.fill();
        }
      }

      ctx.globalAlpha = 1;
      renderLegend();
    }

    function renderLegend() {
      if (state.activeGene !== "__NONE__") {
        legend.innerHTML = `
          <div style=\"font-weight:600; margin-bottom:6px;\">Gene: ${state.activeGene}</div>
          <div style=\"height:10px; background:linear-gradient(90deg, rgb(20,44,89), rgb(239,138,98)); border-radius:5px; margin-bottom:4px;\"></div>
          <div style=\"display:flex; justify-content:space-between; gap:8px;\"><span>${state.currentMin.toFixed(2)}</span><span>${state.currentMax.toFixed(2)}</span></div>
        `;
        return;
      }

      const values = Array.from(state.currentColorMap.entries());
      const top = values.slice(0, 20)
        .map(([label, color]) => `<div class=\"legend-item\"><span class=\"swatch\" style=\"background:${color}\"></span><span>${label || "(empty)"}</span></div>`)
        .join("");

      legend.innerHTML = `<div style=\"font-weight:600; margin-bottom:6px;\">${state.annotationCol}</div>${top}`;
    }

    function updateAnnotationColors() {
      const categories = new Map();
      for (let i = 0; i < state.obs.length; i += 1) {
        const value = annotationValue(i);
        if (!categories.has(value)) {
          categories.set(value, palette[categories.size % palette.length]);
        }
      }
      state.currentColorMap = categories;
      refreshFilterOptions();
    }

    function refreshFilterOptions() {
      const current = state.filterValue;
      filterSelect.innerHTML = "";
      filterSelect.add(new Option("All", "__ALL__"));
      for (const key of state.currentColorMap.keys()) {
        filterSelect.add(new Option(key || "(empty)", key));
      }
      filterSelect.value = state.currentColorMap.has(current) || current === "__ALL__" ? current : "__ALL__";
      state.filterValue = filterSelect.value;
    }

    async function ensureGeneLoaded(gene) {
      if (state.geneCache.has(gene)) {
        return state.geneCache.get(gene);
      }

      const spec = state.manifest.expression.genes.find((x) => x.gene === gene);
      if (!spec) {
        throw new Error(`Gene asset not found: ${gene}`);
      }

      setStatus(`Loading gene ${gene}...`);
      const out = new Float32Array(state.manifest.expression.n_cells);
      for (const chunk of spec.chunks) {
        const buffer = await readBinary(chunk.path);
        const arr = new Float32Array(buffer);
        out.set(arr, chunk.start);
      }

      state.geneCache.set(gene, out);
      return out;
    }

    async function onGeneChange() {
      const gene = geneSelect.value;
      state.activeGene = gene;

      if (gene === "__NONE__") {
        setStatus("Showing annotation colors.");
        draw();
        return;
      }

      try {
        const vec = await ensureGeneLoaded(gene);
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        for (let i = 0; i < vec.length; i += 1) {
          const v = vec[i];
          if (v < min) min = v;
          if (v > max) max = v;
        }
        state.currentMin = Number.isFinite(min) ? min : 0;
        state.currentMax = Number.isFinite(max) ? max : 1;
        setStatus(`Showing gene expression for ${gene}.`);
        draw();
      } catch (err) {
        setStatus(err.message || String(err), true);
      }
    }

    function setupControls() {
      annotationSelect.innerHTML = "";
      for (const col of state.manifest.annotation_columns) {
        annotationSelect.add(new Option(col, col));
      }
      state.annotationCol = state.manifest.annotation_columns[0] || null;
      annotationSelect.value = state.annotationCol || "";
      updateAnnotationColors();

      geneSelect.innerHTML = "";
      geneSelect.add(new Option("(None)", "__NONE__"));
      for (const gene of state.manifest.gene_list) {
        geneSelect.add(new Option(gene, gene));
      }
      geneSelect.value = "__NONE__";

      annotationSelect.addEventListener("change", () => {
        state.annotationCol = annotationSelect.value;
        state.activeGene = "__NONE__";
        geneSelect.value = "__NONE__";
        updateAnnotationColors();
        draw();
      });

      filterSelect.addEventListener("change", () => {
        state.filterValue = filterSelect.value;
        draw();
      });

      geneSelect.addEventListener("change", onGeneChange);

      sizeSlider.addEventListener("input", () => {
        state.pointSize = Number(sizeSlider.value);
        rebuildHoverGrid();
        draw();
      });

      opacitySlider.addEventListener("input", () => {
        state.opacity = Number(opacitySlider.value);
        draw();
      });
    }

    function nearestPoint(mouseX, mouseY) {
      const grid = state.hoverGrid;
      if (!grid) return null;

      const gx = Math.floor(mouseX / grid.cell);
      const gy = Math.floor(mouseY / grid.cell);

      let best = null;
      let bestDist = Math.max(16, state.pointSize * 6);

      for (let ox = -1; ox <= 1; ox += 1) {
        for (let oy = -1; oy <= 1; oy += 1) {
          const key = `${gx + ox}:${gy + oy}`;
          const bucket = grid.map.get(key);
          if (!bucket) continue;

          for (const i of bucket) {
            if (!isVisible(i)) continue;
            const dx = state.px[i] - mouseX;
            const dy = state.py[i] - mouseY;
            const dist = Math.hypot(dx, dy);
            if (dist < bestDist) {
              bestDist = dist;
              best = i;
            }
          }
        }
      }

      return best;
    }

    function bindHover() {
      stage.addEventListener("mousemove", (ev) => {
        const rect = canvas.getBoundingClientRect();
        const mx = ev.clientX - rect.left;
        const my = ev.clientY - rect.top;
        const idx = nearestPoint(mx, my);

        if (idx == null) {
          tooltip.style.display = "none";
          return;
        }

        const row = state.obs[idx];
        const gene = state.activeGene !== "__NONE__" ? state.activeGene : null;
        const geneVal = gene ? state.geneCache.get(gene)?.[idx] : null;

        let html = `<div><strong>${row.cell_id}</strong></div>`;
        for (const col of state.manifest.annotation_columns) {
          html += `<div>${col}: ${row[col]}</div>`;
        }
        if (gene) {
          html += `<div>${gene}: ${Number(geneVal || 0).toFixed(3)}</div>`;
        }

        tooltip.innerHTML = html;
        tooltip.style.display = "block";
        tooltip.style.left = `${mx + 14}px`;
        tooltip.style.top = `${my + 14}px`;
      });

      stage.addEventListener("mouseleave", () => {
        tooltip.style.display = "none";
      });
    }

    async function boot() {
      try {
        setStatus("Loading manifest...");
        state.manifest = await readJSON(MANIFEST_PATH);

        datasetMeta.textContent = `${state.manifest.dataset_name} | cells: ${state.manifest.n_cells.toLocaleString()} | genes: ${state.manifest.n_genes_exported.toLocaleString()}`;

        setStatus("Loading metadata...");
        const [obs, coords, image] = await Promise.all([
          loadObs(),
          loadCoords(),
          loadImage()
        ]);

        state.obs = obs;
        state.coords = coords;
        state.image = image;

        if (state.obs.length !== state.coords.length) {
          throw new Error(`obs rows (${state.obs.length}) and coords rows (${state.coords.length}) do not match.`);
        }

        setupControls();
        bindHover();
        resizeCanvas();
        window.addEventListener("resize", resizeCanvas);
        setStatus("Ready.");
      } catch (err) {
        setStatus(err.message || String(err), true);
      }
    }

    boot();
  </script>
</body>
</html>
"""


def export_h5ad(config: ExportConfig) -> Manifest:
    h5ad_path = Path(config.h5ad_path)
    outdir = Path(config.outdir)
    assets_dir = outdir / "assets"
    expr_dir = assets_dir / "expression"

    outdir.mkdir(parents=True, exist_ok=True)
    assets_dir.mkdir(parents=True, exist_ok=True)
    expr_dir.mkdir(parents=True, exist_ok=True)

    adata = _read_adata(h5ad_path)

    try:
        obs_idx = _resolve_obs_indices(adata.n_obs, config.downsample)
        coord_arr, coord_source = _resolve_coordinates(adata, obs_idx=obs_idx, requested=config.coords)
        anno_cols = _resolve_annotation_columns(adata.obs, config.annotation_columns)

        obs_export = _build_obs_export(adata.obs, obs_idx=obs_idx, anno_cols=anno_cols)
        genes_mode = parse_genes_mode(config.genes_mode)
        selected_genes = resolve_gene_names(adata, genes_mode, obs_indices=obs_idx)

        if not selected_genes:
            raise ValueError("No genes selected for export. Adjust --genes mode.")

        var_export = _build_var_export(adata, selected_genes=selected_genes)

        max_asset_bytes = _max_asset_bytes(config.max_asset_mb)
        asset_files: list[AssetFile] = []

        obs_csv_paths = write_dataframe_chunks(
            obs_export,
            outdir=assets_dir,
            stem="obs",
            max_asset_bytes=max_asset_bytes,
            gzip_enabled=config.gzip,
            fmt="csv",
        )
        for path in obs_csv_paths:
            asset_files.append(
                AssetFile(kind="obs_csv", path=path.relative_to(outdir).as_posix(), size_bytes=file_size(path))
            )

        obs_json_paths = write_dataframe_chunks(
            obs_export,
            outdir=assets_dir,
            stem="obs_records",
            max_asset_bytes=max_asset_bytes,
            gzip_enabled=config.gzip,
            fmt="json",
        )
        for path in obs_json_paths:
            asset_files.append(
                AssetFile(kind="obs_json", path=path.relative_to(outdir).as_posix(), size_bytes=file_size(path))
            )

        coords_df = pd.DataFrame(
            {
                "cell_id": obs_export["cell_id"],
                "x": coord_arr[:, 0],
                "y": coord_arr[:, 1],
            }
        )
        coords_paths = write_dataframe_chunks(
            coords_df,
            outdir=assets_dir,
            stem="coords",
            max_asset_bytes=max_asset_bytes,
            gzip_enabled=config.gzip,
            fmt="json",
        )
        for path in coords_paths:
            asset_files.append(
                AssetFile(kind="coords_json", path=path.relative_to(outdir).as_posix(), size_bytes=file_size(path))
            )

        var_paths = write_dataframe_chunks(
            var_export,
            outdir=assets_dir,
            stem="var",
            max_asset_bytes=max_asset_bytes,
            gzip_enabled=config.gzip,
            fmt="csv",
        )
        for path in var_paths:
            asset_files.append(
                AssetFile(kind="var_csv", path=path.relative_to(outdir).as_posix(), size_bytes=file_size(path))
            )

        expr_spec, expr_assets = _write_expression_assets(
            adata,
            obs_idx=obs_idx,
            selected_genes=selected_genes,
            outdir=outdir,
            expr_dir=expr_dir,
            gzip_assets=config.gzip,
            max_asset_bytes=max_asset_bytes,
        )
        asset_files.extend(expr_assets)

        image_spec, image_assets = _export_image(
            adata,
            outdir=outdir,
            assets_dir=assets_dir,
            image_path=Path(config.image_path) if config.image_path else None,
        )
        asset_files.extend(image_assets)

        preview_path = _write_preview(outdir, coord_arr) if config.preview else None
        if preview_path is not None:
            preview_abs = outdir / preview_path
            asset_files.append(
                AssetFile(kind="preview", path=preview_path, size_bytes=file_size(preview_abs))
            )

        index_path = outdir / "index.html"
        index_path.write_text(_build_viewer_html(), encoding="utf-8")

        manifest = Manifest(
            dataset_name=config.dataset_name or h5ad_path.stem,
            n_cells=len(obs_idx),
            n_genes_exported=len(selected_genes),
            coordinate_source=coord_source,
            annotation_columns=anno_cols,
            gene_list=selected_genes,
            asset_files=asset_files,
            expression=expr_spec,
            preview_path=preview_path,
            image=image_spec,
        )

        write_manifest(outdir / "manifest.json", manifest)

        return manifest
    finally:
        _close_adata(adata)
