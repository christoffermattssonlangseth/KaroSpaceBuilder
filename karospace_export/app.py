from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from functools import partial
from pathlib import Path
import http.server
import os
import queue
import socketserver
import subprocess
import sys
import threading
import traceback
import webbrowser

import anndata as ad

from .export import ExportConfig, export_h5ad

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, ttk
except Exception as exc:  # pragma: no cover - platform/runtime guard
    tk = None
    filedialog = None
    messagebox = None
    ttk = None
    TK_IMPORT_ERROR = exc
else:
    TK_IMPORT_ERROR = None


@dataclass(slots=True)
class AppResult:
    outdir: Path
    n_cells: int
    n_genes_exported: int


class _ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


class ExportApp(tk.Tk if tk is not None else object):
    def __init__(self) -> None:
        super().__init__()
        self.title("KaroSpace Exporter")
        self.geometry("1180x760")
        self.minsize(1024, 700)

        self._queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self._export_thread: threading.Thread | None = None
        self._server: _ThreadingHTTPServer | None = None
        self._server_thread: threading.Thread | None = None
        self._last_outdir: Path | None = None

        self._build_style()
        self._build_variables()
        self._build_layout()
        self.after(120, self._poll_events)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_style(self) -> None:
        self.configure(bg="#efe9dd")

        style = ttk.Style(self)
        style.theme_use("clam")

        style.configure("Root.TFrame", background="#efe9dd")
        style.configure("Card.TFrame", background="#fffdf7")
        style.configure("Header.TLabel", background="#fffdf7", foreground="#1f2933", font=("Avenir Next", 20, "bold"))
        style.configure("Subheader.TLabel", background="#fffdf7", foreground="#52606d", font=("Avenir Next", 10))
        style.configure("FieldLabel.TLabel", background="#fffdf7", foreground="#334e68", font=("Avenir Next", 10, "bold"))
        style.configure("Body.TLabel", background="#fffdf7", foreground="#243b53", font=("Avenir Next", 10))

        style.configure("Primary.TButton", padding=(10, 8), font=("Avenir Next", 10, "bold"))
        style.map(
            "Primary.TButton",
            background=[("!disabled", "#2f855a"), ("active", "#276749")],
            foreground=[("!disabled", "#ffffff")],
        )

        style.configure("Secondary.TButton", padding=(10, 8), font=("Avenir Next", 10))
        style.map(
            "Secondary.TButton",
            background=[("!disabled", "#edf2f7"), ("active", "#e2e8f0")],
            foreground=[("!disabled", "#1a202c")],
        )

        style.configure("Good.Horizontal.TProgressbar", troughcolor="#d9e2ec", background="#2f855a")

    def _build_variables(self) -> None:
        self.h5ad_var = tk.StringVar()
        self.outdir_var = tk.StringVar()
        self.coords_var = tk.StringVar(value="auto")
        self.anno_var = tk.StringVar(value="cell_type,leiden")

        self.genes_mode_var = tk.StringVar(value="hvgs")
        self.genes_count_var = tk.StringVar(value="500")
        self.gene_list_path_var = tk.StringVar()

        self.image_var = tk.StringVar()
        self.downsample_var = tk.StringVar()
        self.gzip_var = tk.BooleanVar(value=True)
        self.max_asset_mb_var = tk.StringVar(value="16")
        self.preview_var = tk.BooleanVar(value=True)

        self.serve_var = tk.BooleanVar(value=False)
        self.port_var = tk.StringVar(value="8000")

        self.status_var = tk.StringVar(value="Ready")

    def _build_layout(self) -> None:
        root = ttk.Frame(self, style="Root.TFrame", padding=16)
        root.pack(fill="both", expand=True)

        root.columnconfigure(0, weight=3)
        root.columnconfigure(1, weight=2)
        root.rowconfigure(0, weight=1)

        controls = ttk.Frame(root, style="Card.TFrame", padding=18)
        controls.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        side = ttk.Frame(root, style="Card.TFrame", padding=18)
        side.grid(row=0, column=1, sticky="nsew")

        controls.columnconfigure(1, weight=1)
        side.columnconfigure(0, weight=1)
        side.rowconfigure(4, weight=1)

        ttk.Label(controls, text="KaroSpace Export App", style="Header.TLabel").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Label(
            controls,
            text="Create a fully static viewer export from AnnData without memorizing CLI flags.",
            style="Subheader.TLabel",
        ).grid(row=1, column=0, columnspan=3, sticky="w", pady=(2, 18))

        row = 2
        row = self._path_field(controls, row, "Input .h5ad", self.h5ad_var, choose_file=True)
        row = self._path_field(controls, row, "Output directory", self.outdir_var, choose_file=False)

        row = self._option_row(
            controls,
            row,
            "Coordinates",
            widget=self._coords_dropdown(controls),
            hint="Auto detects obsm['spatial'] or obs centroid_x/centroid_y.",
        )
        row = self._option_row(
            controls,
            row,
            "Annotation columns",
            widget=self._entry(controls, self.anno_var),
            hint="Comma separated list, e.g. cell_type,leiden. Leave blank for auto-select.",
        )

        row = self._option_row(
            controls,
            row,
            "Genes mode",
            widget=self._genes_mode_row(controls),
            hint="Choose HVGs/top mean with a count, or provide a gene list file.",
        )

        row = self._path_field(controls, row, "Optional tissue image", self.image_var, choose_file=True, optional=True)

        downsample_container = ttk.Frame(controls, style="Card.TFrame")
        ttk.Entry(downsample_container, textvariable=self.downsample_var, width=12).pack(side="left")
        ttk.Label(downsample_container, text="cells (blank = all)", style="Body.TLabel").pack(side="left", padx=(8, 0))
        row = self._option_row(
            controls,
            row,
            "Downsample",
            widget=downsample_container,
            hint="Useful for lightweight demos.",
        )

        max_asset_container = ttk.Frame(controls, style="Card.TFrame")
        ttk.Entry(max_asset_container, textvariable=self.max_asset_mb_var, width=12).pack(side="left")
        ttk.Label(max_asset_container, text="MB per asset file", style="Body.TLabel").pack(side="left", padx=(8, 0))
        row = self._option_row(
            controls,
            row,
            "Asset split limit",
            widget=max_asset_container,
            hint="Large datasets are chunked below this per-file target.",
        )

        toggles = ttk.Frame(controls, style="Card.TFrame")
        ttk.Checkbutton(toggles, text="Gzip assets", variable=self.gzip_var).pack(side="left")
        ttk.Checkbutton(toggles, text="Write preview.png", variable=self.preview_var).pack(side="left", padx=(14, 0))
        ttk.Checkbutton(toggles, text="Serve after export", variable=self.serve_var).pack(side="left", padx=(14, 0))
        row = self._option_row(controls, row, "Options", widget=toggles)

        serve_row = ttk.Frame(controls, style="Card.TFrame")
        ttk.Entry(serve_row, textvariable=self.port_var, width=10).pack(side="left")
        ttk.Label(serve_row, text="port", style="Body.TLabel").pack(side="left", padx=(8, 0))
        row = self._option_row(controls, row, "Serve port", widget=serve_row)

        button_row = ttk.Frame(controls, style="Card.TFrame")
        self.inspect_btn = ttk.Button(button_row, text="Inspect H5AD", style="Secondary.TButton", command=self._inspect_h5ad)
        self.inspect_btn.pack(side="left")

        self.export_btn = ttk.Button(button_row, text="Export", style="Primary.TButton", command=self._on_export)
        self.export_btn.pack(side="left", padx=(10, 0))

        self.stop_server_btn = ttk.Button(button_row, text="Stop Server", style="Secondary.TButton", command=self._stop_server)
        self.stop_server_btn.pack(side="left", padx=(10, 0))

        row = self._option_row(controls, row, "Actions", widget=button_row)

        ttk.Label(side, text="Runtime", style="Header.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(side, textvariable=self.status_var, style="Subheader.TLabel").grid(row=1, column=0, sticky="w", pady=(2, 12))

        self.progress = ttk.Progressbar(side, mode="indeterminate", style="Good.Horizontal.TProgressbar")
        self.progress.grid(row=2, column=0, sticky="ew")

        launch_row = ttk.Frame(side, style="Card.TFrame")
        launch_row.grid(row=3, column=0, sticky="ew", pady=(10, 10))
        ttk.Button(launch_row, text="Open Output Folder", style="Secondary.TButton", command=self._open_output_folder).pack(
            side="left"
        )
        ttk.Button(launch_row, text="Open Viewer", style="Secondary.TButton", command=self._open_viewer).pack(side="left", padx=(10, 0))

        log_wrap = ttk.Frame(side, style="Card.TFrame")
        log_wrap.grid(row=4, column=0, sticky="nsew")
        log_wrap.columnconfigure(0, weight=1)
        log_wrap.rowconfigure(0, weight=1)

        self.log_text = tk.Text(
            log_wrap,
            height=14,
            wrap="word",
            background="#0f172a",
            foreground="#cbd5e1",
            insertbackground="#cbd5e1",
            relief="flat",
            padx=10,
            pady=10,
        )
        self.log_text.grid(row=0, column=0, sticky="nsew")
        log_scroll = ttk.Scrollbar(log_wrap, orient="vertical", command=self.log_text.yview)
        log_scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=log_scroll.set)
        self.log_text.configure(state="disabled")

        self.genes_mode_var.trace_add("write", lambda *_: self._update_genes_mode_visibility())
        self._update_genes_mode_visibility()

    def _path_field(
        self,
        parent: ttk.Frame,
        row: int,
        label: str,
        variable: tk.StringVar,
        choose_file: bool,
        optional: bool = False,
    ) -> int:
        ttk.Label(parent, text=label, style="FieldLabel.TLabel").grid(row=row, column=0, sticky="nw", pady=(0, 6), padx=(0, 12))
        entry = self._entry(parent, variable)
        entry.grid(row=row, column=1, sticky="ew", pady=(0, 6))

        if choose_file:
            button = ttk.Button(
                parent,
                text="Browse",
                style="Secondary.TButton",
                command=lambda: self._choose_file(variable, optional=optional),
            )
        else:
            button = ttk.Button(parent, text="Browse", style="Secondary.TButton", command=lambda: self._choose_dir(variable))
        button.grid(row=row, column=2, sticky="e", pady=(0, 6), padx=(10, 0))
        return row + 1

    def _option_row(
        self,
        parent: ttk.Frame,
        row: int,
        label: str,
        widget: tk.Widget,
        hint: str | None = None,
    ) -> int:
        ttk.Label(parent, text=label, style="FieldLabel.TLabel").grid(row=row, column=0, sticky="nw", pady=(0, 2), padx=(0, 12))
        widget.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(0, 2))
        row += 1
        if hint:
            ttk.Label(parent, text=hint, style="Subheader.TLabel").grid(row=row, column=1, columnspan=2, sticky="w", pady=(0, 8))
            row += 1
        return row

    def _entry(self, parent: ttk.Frame, variable: tk.StringVar) -> ttk.Entry:
        return ttk.Entry(parent, textvariable=variable)

    def _coords_dropdown(self, parent: ttk.Frame) -> ttk.Combobox:
        combo = ttk.Combobox(
            parent,
            textvariable=self.coords_var,
            values=["auto", "obsm:spatial", "obs:centroid_x_y"],
            state="readonly",
        )
        return combo

    def _genes_mode_row(self, parent: ttk.Frame) -> ttk.Frame:
        wrap = ttk.Frame(parent, style="Card.TFrame")

        self.genes_mode_combo = ttk.Combobox(
            wrap,
            textvariable=self.genes_mode_var,
            values=["hvgs", "top_mean", "list"],
            state="readonly",
            width=12,
        )
        self.genes_mode_combo.pack(side="left")

        self.genes_count_entry = ttk.Entry(wrap, textvariable=self.genes_count_var, width=10)
        self.genes_count_entry.pack(side="left", padx=(10, 0))

        self.genes_count_label = ttk.Label(wrap, text="N genes", style="Body.TLabel")
        self.genes_count_label.pack(side="left", padx=(8, 0))

        self.gene_list_entry = ttk.Entry(wrap, textvariable=self.gene_list_path_var, width=30)
        self.gene_list_button = ttk.Button(
            wrap,
            text="Gene List",
            style="Secondary.TButton",
            command=lambda: self._choose_file(self.gene_list_path_var, optional=False),
        )

        return wrap

    def _update_genes_mode_visibility(self) -> None:
        mode = self.genes_mode_var.get().strip().lower()
        if mode == "list":
            self.genes_count_entry.pack_forget()
            self.genes_count_label.pack_forget()
            if not self.gene_list_entry.winfo_manager():
                self.gene_list_entry.pack(side="left", padx=(10, 0))
                self.gene_list_button.pack(side="left", padx=(10, 0))
            return

        self.gene_list_entry.pack_forget()
        self.gene_list_button.pack_forget()
        if not self.genes_count_entry.winfo_manager():
            self.genes_count_entry.pack(side="left", padx=(10, 0))
            self.genes_count_label.pack(side="left", padx=(8, 0))

    def _choose_file(self, variable: tk.StringVar, optional: bool = False) -> None:
        initial_dir = str(Path(variable.get()).expanduser().parent) if variable.get() else str(Path.cwd())
        path = filedialog.askopenfilename(initialdir=initial_dir)
        if path:
            variable.set(path)
        elif not optional and not variable.get():
            self._log("File selection canceled.")

    def _choose_dir(self, variable: tk.StringVar) -> None:
        initial_dir = variable.get() or str(Path.cwd())
        path = filedialog.askdirectory(initialdir=initial_dir)
        if path:
            variable.set(path)

    def _log(self, message: str) -> None:
        stamp = datetime.now().strftime("%H:%M:%S")
        line = f"[{stamp}] {message}\n"
        self.log_text.configure(state="normal")
        self.log_text.insert("end", line)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def _inspect_h5ad(self) -> None:
        path_text = self.h5ad_var.get().strip()
        if not path_text:
            messagebox.showerror("Missing input", "Pick an input .h5ad first.")
            return

        path = Path(path_text).expanduser()
        if not path.exists():
            messagebox.showerror("Missing file", f"Input file not found:\n{path}")
            return

        self._log(f"Inspecting {path}")
        adata = None
        try:
            try:
                adata = ad.read_h5ad(path, backed="r")
            except Exception:
                adata = ad.read_h5ad(path)

            obs_cols = [str(c) for c in adata.obs.columns]
            if obs_cols:
                defaults = [c for c in obs_cols if c in {"cell_type", "leiden", "sample", "sample_id"}]
                if not defaults:
                    defaults = obs_cols[: min(4, len(obs_cols))]
                self.anno_var.set(",".join(defaults))

            has_spatial = "spatial" in adata.obsm
            has_centroid = {"centroid_x", "centroid_y"}.issubset(set(obs_cols))
            if has_spatial:
                self.coords_var.set("obsm:spatial")
            elif has_centroid:
                self.coords_var.set("obs:centroid_x_y")
            else:
                self.coords_var.set("auto")

            self._log(
                f"obs columns: {len(obs_cols)} | cells: {adata.n_obs} | genes: {adata.n_vars} | "
                f"coords: {'obsm:spatial' if has_spatial else 'obs centroids' if has_centroid else 'not detected'}"
            )
            self.status_var.set("Inspection complete")
        except Exception as exc:
            messagebox.showerror("Inspect failed", str(exc))
            self._log(f"Inspect failed: {exc}")
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    def _parse_config(self) -> ExportConfig:
        h5ad_text = self.h5ad_var.get().strip()
        outdir_text = self.outdir_var.get().strip()

        if not h5ad_text:
            raise ValueError("Input .h5ad is required.")
        if not outdir_text:
            raise ValueError("Output directory is required.")

        h5ad_path = Path(h5ad_text).expanduser()
        outdir = Path(outdir_text).expanduser()
        if not h5ad_path.exists():
            raise ValueError(f"Input .h5ad not found: {h5ad_path}")

        coords = None if self.coords_var.get() == "auto" else self.coords_var.get()
        annos = [x.strip() for x in self.anno_var.get().split(",") if x.strip()]
        if not annos:
            annos = None

        mode = self.genes_mode_var.get().strip().lower()
        if mode in {"hvgs", "top_mean"}:
            count_text = self.genes_count_var.get().strip()
            if not count_text:
                raise ValueError("Set a gene count for hvgs/top_mean mode.")
            try:
                count = int(count_text)
            except ValueError as exc:
                raise ValueError("Gene count must be an integer.") from exc
            if count <= 0:
                raise ValueError("Gene count must be > 0.")
            genes_mode = f"{mode}:{count}"
        elif mode == "list":
            list_path_text = self.gene_list_path_var.get().strip()
            if not list_path_text:
                raise ValueError("Choose a gene list file for list mode.")
            list_path = Path(list_path_text).expanduser()
            if not list_path.exists():
                raise ValueError(f"Gene list file not found: {list_path}")
            genes_mode = f"list:{list_path}"
        else:
            raise ValueError("Genes mode must be hvgs, top_mean, or list.")

        image_text = self.image_var.get().strip()
        image_path = Path(image_text).expanduser() if image_text else None

        downsample_text = self.downsample_var.get().strip()
        downsample = None
        if downsample_text:
            try:
                downsample = int(downsample_text)
            except ValueError as exc:
                raise ValueError("Downsample must be an integer.") from exc
            if downsample <= 0:
                raise ValueError("Downsample must be > 0.")

        max_asset_text = self.max_asset_mb_var.get().strip()
        try:
            max_asset_mb = float(max_asset_text)
        except ValueError as exc:
            raise ValueError("Asset split limit must be a number.") from exc
        if max_asset_mb <= 0:
            raise ValueError("Asset split limit must be > 0.")

        return ExportConfig(
            h5ad_path=h5ad_path,
            outdir=outdir,
            coords=coords,
            annotation_columns=annos,
            genes_mode=genes_mode,
            image_path=image_path,
            downsample=downsample,
            gzip=bool(self.gzip_var.get()),
            max_asset_mb=max_asset_mb,
            preview=bool(self.preview_var.get()),
        )

    def _set_busy(self, busy: bool) -> None:
        widgets = [self.export_btn, self.inspect_btn]
        for widget in widgets:
            if busy:
                widget.state(["disabled"])
            else:
                widget.state(["!disabled"])

        if busy:
            self.progress.start(12)
            self.status_var.set("Export running...")
        else:
            self.progress.stop()
            self.status_var.set("Ready")

    def _on_export(self) -> None:
        if self._export_thread and self._export_thread.is_alive():
            messagebox.showinfo("Export running", "An export is already running.")
            return

        try:
            config = self._parse_config()
        except Exception as exc:
            messagebox.showerror("Invalid options", str(exc))
            return

        self._set_busy(True)
        self._log(f"Starting export: {config.h5ad_path} -> {config.outdir}")

        thread = threading.Thread(target=self._run_export, args=(config,), daemon=True)
        self._export_thread = thread
        thread.start()

    def _run_export(self, config: ExportConfig) -> None:
        try:
            manifest = export_h5ad(config)
            result = AppResult(outdir=config.outdir, n_cells=manifest.n_cells, n_genes_exported=manifest.n_genes_exported)
            self._queue.put(("done", result))
            if self.serve_var.get():
                self._queue.put(("start_server", config.outdir))
        except Exception:
            self._queue.put(("error", traceback.format_exc()))

    def _poll_events(self) -> None:
        while True:
            try:
                kind, payload = self._queue.get_nowait()
            except queue.Empty:
                break

            if kind == "done":
                self._set_busy(False)
                result = payload
                assert isinstance(result, AppResult)
                self._last_outdir = result.outdir
                self._log(
                    f"Export complete. cells={result.n_cells}, genes={result.n_genes_exported}, outdir={result.outdir}"
                )
                self.status_var.set("Export complete")
            elif kind == "start_server":
                outdir = payload
                assert isinstance(outdir, Path)
                self._start_server(outdir)
            elif kind == "error":
                self._set_busy(False)
                details = str(payload)
                self._log("Export failed. See traceback in popup.")
                self.status_var.set("Export failed")
                messagebox.showerror("Export failed", details)

        self.after(120, self._poll_events)

    def _start_server(self, outdir: Path | None = None) -> None:
        target = outdir or self._last_outdir
        if target is None:
            messagebox.showinfo("No export", "Run an export first.")
            return

        try:
            port = int(self.port_var.get().strip())
        except ValueError:
            messagebox.showerror("Invalid port", "Port must be an integer.")
            return

        if self._server is not None:
            self._stop_server()

        handler = partial(http.server.SimpleHTTPRequestHandler, directory=str(target))

        try:
            server = _ThreadingHTTPServer(("127.0.0.1", port), handler)
        except OSError as exc:
            messagebox.showerror("Server start failed", str(exc))
            self._log(f"Could not start server: {exc}")
            return

        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()

        self._server = server
        self._server_thread = thread
        self._log(f"Serving {target} at http://127.0.0.1:{port}")
        self.status_var.set(f"Serving on :{port}")
        webbrowser.open_new_tab(f"http://127.0.0.1:{port}")

    def _stop_server(self) -> None:
        if self._server is None:
            return

        self._server.shutdown()
        self._server.server_close()
        self._server = None
        self._server_thread = None
        self._log("Preview server stopped")
        if not (self._export_thread and self._export_thread.is_alive()):
            self.status_var.set("Ready")

    def _open_output_folder(self) -> None:
        path = self._last_outdir or (Path(self.outdir_var.get().strip()).expanduser() if self.outdir_var.get().strip() else None)
        if path is None:
            messagebox.showinfo("No output", "Pick an output directory first.")
            return
        self._open_path(path)

    def _open_viewer(self) -> None:
        if self._server is not None:
            try:
                port = int(self.port_var.get().strip())
            except ValueError:
                port = 8000
            webbrowser.open_new_tab(f"http://127.0.0.1:{port}")
            return

        outdir = self._last_outdir or (Path(self.outdir_var.get().strip()).expanduser() if self.outdir_var.get().strip() else None)
        if outdir is None:
            messagebox.showinfo("No output", "Run an export first.")
            return

        viewer = outdir / "index.html"
        if not viewer.exists():
            messagebox.showerror("Viewer missing", f"Expected file not found:\n{viewer}")
            return

        webbrowser.open_new_tab(viewer.resolve().as_uri())

    def _open_path(self, path: Path) -> None:
        try:
            if sys.platform == "darwin":
                subprocess.Popen(["open", str(path)])
            elif os.name == "nt":  # pragma: no cover - platform branch
                os.startfile(str(path))
            else:
                subprocess.Popen(["xdg-open", str(path)])
        except Exception as exc:
            messagebox.showerror("Open failed", str(exc))
            self._log(f"Open path failed: {exc}")

    def _on_close(self) -> None:
        self._stop_server()
        self.destroy()


def main() -> int:
    if tk is None:
        raise RuntimeError(
            "Tkinter is not available in this Python environment. "
            f"Original import error: {TK_IMPORT_ERROR}"
        )

    app = ExportApp()
    app.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
