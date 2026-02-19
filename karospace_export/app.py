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
import tempfile
import threading
import time
import traceback
import webbrowser

from .utils import parse_genes_mode, resolve_gene_names

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

_TTK_FRAME_BASE = ttk.Frame if ttk is not None else object

_KI_COLORS = {
    "plum_dark": "#4F0433",
    "orange": "#FF876F",
    "light_orange": "#FEEEEB",
    "light_blue": "#EDF4F4",
    "plum": "#870052",
}

_KAROSPACE_LIGHT_PALETTE = {
    "plum_dark": _KI_COLORS["plum_dark"],
    "orange": _KI_COLORS["orange"],
    "light_orange": _KI_COLORS["light_orange"],
    "light_blue": _KI_COLORS["light_blue"],
    "plum": _KI_COLORS["plum"],
    "background": "#f5f5f5",
    "text": "#1a1a1a",
    "panel_bg": "#ffffff",
    "border": "#e0e0e0",
    "input_bg": "#ffffff",
    "muted": "#666666",
    "hover_bg": "#f0f0f0",
    "accent": _KI_COLORS["plum"],
    "accent_strong": _KI_COLORS["plum_dark"],
    "on_accent": "#ffffff",
}

_KAROSPACE_DARK_PALETTE = {
    "plum_dark": _KI_COLORS["plum_dark"],
    "orange": _KI_COLORS["orange"],
    "light_orange": _KI_COLORS["light_orange"],
    "light_blue": _KI_COLORS["light_blue"],
    "plum": _KI_COLORS["plum"],
    "background": "#171717",
    "text": "#e9e9e9",
    "panel_bg": "#222222",
    "border": "#3a3a3a",
    "input_bg": "#2a2a2a",
    "muted": "#a0a0a0",
    "hover_bg": "#303030",
    "accent": _KI_COLORS["orange"],
    "accent_strong": _KI_COLORS["plum"],
    "on_accent": "#1a1a1a",
}


def _palette_for_mode(mode: str) -> dict[str, str]:
    return (
        dict(_KAROSPACE_DARK_PALETTE)
        if str(mode).strip().lower() == "dark"
        else dict(_KAROSPACE_LIGHT_PALETTE)
    )


def _get_anndata():
    import anndata as ad

    return ad


def _get_numpy():
    import numpy as np

    return np


class SearchableListEditor(_TTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        *,
        label: str,
        height: int = 8,
        help_text: str | None = None,
    ) -> None:
        super().__init__(parent, style="Card.TFrame")
        self._choices: list[str] = []
        self._input_var = tk.StringVar(value="")

        self.columnconfigure(0, weight=1)
        ttk.Label(self, text=label, style="FieldLabel.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 6))

        controls = ttk.Frame(self, style="Card.TFrame")
        controls.grid(row=1, column=0, sticky="ew")
        controls.columnconfigure(0, weight=1)

        self.entry = ttk.Combobox(controls, textvariable=self._input_var, state="normal")
        self.entry.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self.entry.bind("<KeyRelease>", self._on_search)
        self.entry.bind("<Return>", lambda _event: self.add_current())

        self.add_btn = ttk.Button(controls, text="+ Add", style="Secondary.TButton", command=self.add_current)
        self.add_btn.grid(row=0, column=1, padx=(0, 6))
        self.remove_btn = ttk.Button(controls, text="Remove", style="Secondary.TButton", command=self.remove_selected)
        self.remove_btn.grid(row=0, column=2, padx=(0, 6))
        self.clear_btn = ttk.Button(controls, text="Clear", style="Secondary.TButton", command=self.clear)
        self.clear_btn.grid(row=0, column=3)

        list_wrap = ttk.Frame(self, style="Card.TFrame")
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=height,
            selectmode="extended",
            activestyle="none",
            background=_KAROSPACE_LIGHT_PALETTE["input_bg"],
            foreground=_KAROSPACE_LIGHT_PALETTE["text"],
            selectbackground=_KAROSPACE_LIGHT_PALETTE["accent"],
            selectforeground=_KAROSPACE_LIGHT_PALETTE["on_accent"],
            relief="solid",
            bd=1,
            highlightthickness=1,
            highlightbackground=_KAROSPACE_LIGHT_PALETTE["border"],
            highlightcolor=_KAROSPACE_LIGHT_PALETTE["accent"],
            font=("Segoe UI", 10),
        )
        self.listbox.grid(row=0, column=0, sticky="ew")
        scroll = ttk.Scrollbar(list_wrap, orient="vertical", command=self.listbox.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.listbox.configure(yscrollcommand=scroll.set)

        if help_text:
            ttk.Label(self, text=help_text, style="Subheader.TLabel").grid(row=3, column=0, sticky="w", pady=(6, 0))

    def _on_search(self, _event) -> None:
        self._update_choices(self._input_var.get())

    def _update_choices(self, query: str = "") -> None:
        needle = query.strip().lower()
        if not needle:
            values = self._choices
        else:
            values = [item for item in self._choices if needle in item.lower()]
        self.entry.configure(values=values[:300])

    def set_choices(self, values: list[str] | tuple[str, ...]) -> None:
        seen: set[str] = set()
        ordered: list[str] = []
        for raw in values:
            value = str(raw).strip()
            if not value or value in seen:
                continue
            seen.add(value)
            ordered.append(value)
        self._choices = ordered
        self._update_choices(self._input_var.get())

    def add_current(self) -> None:
        value = self._input_var.get().strip()
        if not value:
            return
        existing = self.get_items()
        if value in existing:
            idx = existing.index(value)
            self.listbox.selection_clear(0, "end")
            self.listbox.selection_set(idx)
            self.listbox.see(idx)
            self._input_var.set("")
            return
        self.listbox.insert("end", value)
        self._input_var.set("")
        self._update_choices("")

    def remove_selected(self) -> None:
        for idx in reversed(self.listbox.curselection()):
            self.listbox.delete(idx)

    def clear(self) -> None:
        self.listbox.delete(0, "end")

    def set_items(self, values: list[str] | tuple[str, ...]) -> None:
        self.clear()
        seen: set[str] = set()
        for raw in values:
            value = str(raw).strip()
            if not value or value in seen:
                continue
            seen.add(value)
            self.listbox.insert("end", value)

    def get_items(self) -> list[str]:
        return [str(v) for v in self.listbox.get(0, "end")]

    def set_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        self.entry.configure(state=state)
        self.add_btn.configure(state=state)
        self.remove_btn.configure(state=state)
        self.clear_btn.configure(state=state)
        self.listbox.configure(state=state)

    def apply_palette(self, palette: dict[str, str]) -> None:
        self.listbox.configure(
            background=palette["input_bg"],
            foreground=palette["text"],
            selectbackground=palette["accent"],
            selectforeground=palette.get("on_accent", "#ffffff"),
            disabledforeground=palette["muted"],
            highlightbackground=palette["border"],
            highlightcolor=palette["accent"],
        )


class SearchableMultiSelectEditor(_TTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        *,
        label: str,
        height: int = 8,
        help_text: str | None = None,
        max_visible: int = 500,
    ) -> None:
        super().__init__(parent, style="Card.TFrame")
        self._choices: list[str] = []
        self._selected_values: set[str] = set()
        self._search_var = tk.StringVar(value="")
        self._info_var = tk.StringVar(value="")
        self._max_visible = max(50, int(max_visible))

        self.columnconfigure(0, weight=1)
        ttk.Label(self, text=label, style="FieldLabel.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 6))

        controls = ttk.Frame(self, style="Card.TFrame")
        controls.grid(row=1, column=0, sticky="ew")
        controls.columnconfigure(0, weight=1)

        self.search_entry = ttk.Entry(controls, textvariable=self._search_var)
        self.search_entry.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self.search_entry.bind("<KeyRelease>", lambda _event: self._on_search_changed())

        self.select_all_btn = ttk.Button(
            controls,
            text="Select all",
            style="Secondary.TButton",
            command=self.select_all_matches,
        )
        self.select_all_btn.grid(row=0, column=1, padx=(0, 6))
        self.clear_btn = ttk.Button(
            controls,
            text="Clear",
            style="Secondary.TButton",
            command=self.clear_selection,
        )
        self.clear_btn.grid(row=0, column=2)

        list_wrap = ttk.Frame(self, style="Card.TFrame")
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=height,
            selectmode="extended",
            activestyle="none",
            background=_KAROSPACE_LIGHT_PALETTE["input_bg"],
            foreground=_KAROSPACE_LIGHT_PALETTE["text"],
            selectbackground=_KAROSPACE_LIGHT_PALETTE["accent"],
            selectforeground=_KAROSPACE_LIGHT_PALETTE["on_accent"],
            relief="solid",
            bd=1,
            highlightthickness=1,
            highlightbackground=_KAROSPACE_LIGHT_PALETTE["border"],
            highlightcolor=_KAROSPACE_LIGHT_PALETTE["accent"],
            font=("Segoe UI", 10),
        )
        self.listbox.grid(row=0, column=0, sticky="ew")
        self.listbox.bind("<<ListboxSelect>>", lambda _event: self._capture_visible_selection())
        scroll = ttk.Scrollbar(list_wrap, orient="vertical", command=self.listbox.yview)
        scroll.grid(row=0, column=1, sticky="ns")
        self.listbox.configure(yscrollcommand=scroll.set)

        ttk.Label(self, textvariable=self._info_var, style="Subheader.TLabel").grid(row=3, column=0, sticky="w", pady=(6, 0))

        if help_text:
            ttk.Label(self, text=help_text, style="Subheader.TLabel").grid(row=4, column=0, sticky="w", pady=(4, 0))

    def _all_matches(self, needle: str) -> list[str]:
        if not needle:
            return list(self._choices)
        return [name for name in self._choices if needle in name.lower()]

    def _filtered_choices(self) -> tuple[list[str], int]:
        needle = self._search_var.get().strip().lower()
        matches = self._all_matches(needle)
        return matches[: self._max_visible], len(matches)

    def _capture_visible_selection(self) -> None:
        visible = [str(v) for v in self.listbox.get(0, "end")]
        if not visible:
            return
        visible_set = set(visible)
        self._selected_values = {name for name in self._selected_values if name not in visible_set}
        for index in self.listbox.curselection():
            if 0 <= int(index) < len(visible):
                self._selected_values.add(visible[int(index)])

    def _render(self) -> None:
        visible, total_matches = self._filtered_choices()
        self.listbox.delete(0, "end")
        for idx, name in enumerate(visible):
            self.listbox.insert("end", name)
            if name in self._selected_values:
                self.listbox.selection_set(idx)
        if total_matches > len(visible):
            self._info_var.set(
                f"Showing first {len(visible)} of {total_matches} matches. Type more characters to narrow."
            )
        else:
            self._info_var.set(f"{total_matches} matches.")

    def _on_search_changed(self) -> None:
        self._capture_visible_selection()
        self._render()

    def set_choices(self, values: list[str] | tuple[str, ...]) -> None:
        self._capture_visible_selection()
        seen: set[str] = set()
        ordered: list[str] = []
        for raw in values:
            value = str(raw).strip()
            if not value or value in seen:
                continue
            seen.add(value)
            ordered.append(value)
        self._choices = ordered
        allowed = set(self._choices)
        self._selected_values = {name for name in self._selected_values if name in allowed}
        self._render()

    def set_selected(self, values: list[str] | tuple[str, ...]) -> None:
        allowed = set(self._choices)
        self._selected_values = {str(v).strip() for v in values if str(v).strip() in allowed}
        self._render()

    def get_selected(self) -> list[str]:
        self._capture_visible_selection()
        selected = self._selected_values
        return [name for name in self._choices if name in selected]

    def select_all_matches(self) -> None:
        self._capture_visible_selection()
        needle = self._search_var.get().strip().lower()
        for name in self._all_matches(needle):
            self._selected_values.add(name)
        self._render()

    def clear_selection(self) -> None:
        self._selected_values.clear()
        self._render()

    def set_enabled(self, enabled: bool) -> None:
        state = "normal" if enabled else "disabled"
        self.search_entry.configure(state=state)
        self.select_all_btn.configure(state=state)
        self.clear_btn.configure(state=state)
        self.listbox.configure(state=state)

    def apply_palette(self, palette: dict[str, str]) -> None:
        self.listbox.configure(
            background=palette["input_bg"],
            foreground=palette["text"],
            selectbackground=palette["accent"],
            selectforeground=palette.get("on_accent", "#ffffff"),
            disabledforeground=palette["muted"],
            highlightbackground=palette["border"],
            highlightcolor=palette["accent"],
        )


@dataclass(slots=True)
class AppResult:
    outdir: Path
    n_cells: int
    n_sections: int
    output_html: Path


@dataclass(slots=True)
class BuilderConfig:
    h5ad_path: Path
    outdir: Path
    coords_mode: str | None
    section_groupby: str
    initial_color: str
    title: str
    theme: str
    outline_by: str | None
    min_panel_size: int
    spot_size: float | str | None
    downsample: int | None
    additional_colors: list[str] | None
    genes: list[str] | None
    use_hvgs: bool
    hvg_limit: int
    marker_genes_groupby: list[str] | None
    genes_mode: str
    marker_genes_top_n: int
    neighbor_stats_groupby: list[str] | None
    neighbor_stats_permutations: int | None
    neighbor_stats_seed: int
    interaction_markers_enabled: bool
    interaction_markers_groupby: list[str] | None
    interaction_markers_top_targets: int
    interaction_markers_top_genes: int
    interaction_markers_min_cells: int
    interaction_markers_min_neighbors: int


class _ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


class ExportApp(tk.Tk if tk is not None else object):
    _OUTPUT_HTML_BASENAME = "KaroSpace"
    _NEIGHBOR_MATRIX_BUDGET_BYTES = 512 * 1024 * 1024
    _MARKER_GROUPBY_MAX_UNIQUE = 2048
    _INTERACTION_GROUPBY_MAX_UNIQUE = 512

    def __init__(self) -> None:
        super().__init__()
        self.title("KaroSpaceBuilder")
        self.geometry("1180x760")
        self.minsize(1024, 700)

        self._queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self._export_thread: threading.Thread | None = None
        self._server: _ThreadingHTTPServer | None = None
        self._server_thread: threading.Thread | None = None
        self._last_outdir: Path | None = None
        self._last_output_html: Path | None = None
        self._inspected_h5ad_path: Path | None = None
        self._inspected_var_name_set: set[str] | None = None
        self._inspected_coords_mode: str | None = None

        self._build_style()
        self._build_variables()
        self._build_layout()
        self.after(120, self._poll_events)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_style(self) -> None:
        self._style = ttk.Style(self)
        self._style.theme_use("clam")
        self._app_palette = _palette_for_mode("dark")
        self._apply_app_theme("dark")

    def _apply_app_theme(self, mode: str) -> None:
        palette = _palette_for_mode(mode)
        self._app_palette = palette
        self.configure(bg=palette["background"])

        style = self._style
        style.configure("Root.TFrame", background=palette["background"])
        style.configure("Card.TFrame", background=palette["panel_bg"])
        style.configure("Header.TLabel", background=palette["panel_bg"], foreground=palette["text"], font=("Segoe UI", 20, "bold"))
        style.configure("Subheader.TLabel", background=palette["panel_bg"], foreground=palette["muted"], font=("Segoe UI", 10))
        style.configure("FieldLabel.TLabel", background=palette["panel_bg"], foreground=palette["text"], font=("Segoe UI", 10, "bold"))
        style.configure("Body.TLabel", background=palette["panel_bg"], foreground=palette["text"], font=("Segoe UI", 10))

        style.configure(
            "Primary.TButton",
            padding=(10, 8),
            font=("Segoe UI", 10, "bold"),
            borderwidth=1,
            relief="solid",
            background=palette["accent"],
            foreground=palette.get("on_accent", "#ffffff"),
        )
        style.map(
            "Primary.TButton",
            background=[("disabled", palette["hover_bg"]), ("active", palette["accent_strong"]), ("!disabled", palette["accent"])],
            foreground=[("disabled", palette["muted"]), ("!disabled", palette.get("on_accent", "#ffffff"))],
        )

        style.configure(
            "Secondary.TButton",
            padding=(10, 8),
            font=("Segoe UI", 10),
            borderwidth=1,
            relief="solid",
            background=palette["input_bg"],
            foreground=palette["text"],
        )
        style.map(
            "Secondary.TButton",
            background=[("active", palette["hover_bg"]), ("!disabled", palette["input_bg"])],
            foreground=[("disabled", palette["muted"]), ("!disabled", palette["text"])],
        )

        style.configure("TEntry", fieldbackground=palette["input_bg"], foreground=palette["text"], bordercolor=palette["border"])
        style.configure(
            "TCombobox",
            fieldbackground=palette["input_bg"],
            background=palette["input_bg"],
            foreground=palette["text"],
            bordercolor=palette["border"],
            arrowcolor=palette["text"],
        )
        style.map(
            "TCombobox",
            fieldbackground=[("readonly", palette["input_bg"])],
            background=[("readonly", palette["input_bg"])],
            foreground=[("readonly", palette["text"])],
        )
        style.configure("TCheckbutton", background=palette["panel_bg"], foreground=palette["text"], font=("Segoe UI", 10))
        style.map("TCheckbutton", foreground=[("disabled", palette["muted"]), ("!disabled", palette["text"])])

        style.configure("TNotebook", background=palette["background"], borderwidth=0)
        style.configure(
            "TNotebook.Tab",
            background=palette["hover_bg"],
            foreground=palette["muted"],
            font=("Segoe UI", 10, "bold"),
            padding=(12, 8),
            borderwidth=1,
            relief="flat",
        )
        style.map(
            "TNotebook.Tab",
            background=[("selected", palette["panel_bg"]), ("active", palette["hover_bg"])],
            foreground=[("selected", palette["text"]), ("active", palette["text"])],
        )
        style.configure("Good.Horizontal.TProgressbar", troughcolor=palette["hover_bg"], background=palette["accent"])

        if hasattr(self, "scroll_canvas"):
            self.scroll_canvas.configure(background=palette["background"])
        if hasattr(self, "help_text"):
            self.help_text.configure(
                background=palette["input_bg"],
                foreground=palette["text"],
                insertbackground=palette["text"],
                highlightbackground=palette["border"],
                highlightcolor=palette["accent"],
            )
        if hasattr(self, "log_text"):
            self.log_text.configure(
                background=palette["input_bg"],
                foreground=palette["text"],
                insertbackground=palette["text"],
                highlightbackground=palette["border"],
                highlightcolor=palette["accent"],
            )

        for attr in (
            "additional_colors_editor",
            "groupby_editor",
            "manual_genes_editor",
            "selection_additional_picker",
            "selection_groupby_picker",
            "selection_genes_picker",
        ):
            widget = getattr(self, attr, None)
            if widget is not None and hasattr(widget, "apply_palette"):
                widget.apply_palette(palette)

    def _sync_app_theme_to_viewer_setting(self) -> None:
        self._apply_app_theme(self.theme_var.get())

    def _build_variables(self) -> None:
        self.h5ad_var = tk.StringVar()
        self.outdir_var = tk.StringVar()
        self.coords_var = tk.StringVar(value="auto")
        self.section_groupby_var = tk.StringVar(value="sample_id")
        self.initial_color_var = tk.StringVar(value="leiden")
        self.title_var = tk.StringVar(value="KaroSpace")
        self.theme_var = tk.StringVar(value="dark")
        self.outline_by_var = tk.StringVar(value="condition")

        self.genes_mode_var = tk.StringVar(value="hvgs")
        self.genes_count_var = tk.StringVar(value="500")
        self.gene_list_path_var = tk.StringVar()
        self.advanced_open_var = tk.BooleanVar(value=False)
        self.min_panel_size_var = tk.StringVar(value="120")
        self.spot_size_var = tk.StringVar(value="auto")
        self.marker_genes_top_n_var = tk.StringVar(value="50")
        self.neighbor_permutations_var = tk.StringVar(value="25")
        self.neighbor_stats_seed_var = tk.StringVar(value="42")
        self.neighbor_auto_var = tk.BooleanVar(value=True)
        self.interaction_markers_enabled_var = tk.BooleanVar(value=False)
        self.interaction_markers_top_targets_var = tk.StringVar(value="6")
        self.interaction_markers_top_genes_var = tk.StringVar(value="15")
        self.interaction_markers_min_cells_var = tk.StringVar(value="30")
        self.interaction_markers_min_neighbors_var = tk.StringVar(value="1")
        self.selection_mode_var = tk.BooleanVar(value=False)

        self.downsample_var = tk.StringVar()

        self.serve_var = tk.BooleanVar(value=False)
        self.port_var = tk.StringVar(value="8000")

        self.status_var = tk.StringVar(value="Ready")

    def _build_layout(self) -> None:
        palette = self._app_palette
        shell = ttk.Frame(self, style="Root.TFrame")
        shell.pack(fill="both", expand=True)
        shell.columnconfigure(0, weight=1)
        shell.rowconfigure(0, weight=1)

        self.scroll_canvas = tk.Canvas(shell, background=palette["background"], highlightthickness=0, bd=0)
        self.scroll_canvas.grid(row=0, column=0, sticky="nsew")
        self.scrollbar = ttk.Scrollbar(shell, orient="vertical", command=self.scroll_canvas.yview)
        self.scrollbar.grid(row=0, column=1, sticky="ns")
        self.scroll_canvas.configure(yscrollcommand=self.scrollbar.set)

        root = ttk.Frame(self.scroll_canvas, style="Root.TFrame", padding=16)
        self._scroll_window = self.scroll_canvas.create_window((0, 0), window=root, anchor="nw")
        root.bind("<Configure>", lambda _event: self.scroll_canvas.configure(scrollregion=self.scroll_canvas.bbox("all")))
        self.scroll_canvas.bind("<Configure>", self._on_canvas_configure)
        self._bind_mousewheel_scroll()

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

        ttk.Label(controls, text="KaroSpaceBuilder", style="Header.TLabel").grid(row=0, column=0, columnspan=3, sticky="w")
        ttk.Label(
            controls,
            text="Export AnnData into a static KaroSpace viewer bundle with guided inputs and inspected field pickers.",
            style="Subheader.TLabel",
        ).grid(row=1, column=0, columnspan=3, sticky="w", pady=(2, 18))

        preset_row = ttk.Frame(controls, style="Card.TFrame")
        preset_row.grid(row=2, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        ttk.Button(preset_row, text="Default", style="Secondary.TButton", command=lambda: self._apply_preset("default")).pack(
            side="left"
        )
        ttk.Label(preset_row, text="Default profile for fast setup.", style="Subheader.TLabel").pack(side="left", padx=(12, 0))

        notebook = ttk.Notebook(controls)
        notebook.grid(row=3, column=0, columnspan=3, sticky="nsew")
        controls.rowconfigure(3, weight=1)

        basic_tab = ttk.Frame(notebook, style="Card.TFrame", padding=10)
        colors_tab = ttk.Frame(notebook, style="Card.TFrame", padding=10)
        advanced_tab = ttk.Frame(notebook, style="Card.TFrame", padding=10)
        help_tab = ttk.Frame(notebook, style="Card.TFrame", padding=10)
        notebook.add(basic_tab, text="Basic")
        notebook.add(colors_tab, text="Colors & Genes")
        notebook.add(advanced_tab, text="Advanced")
        notebook.add(help_tab, text="Help")

        basic_tab.columnconfigure(1, weight=1)
        row = 0
        row = self._path_field(basic_tab, row, "Input .h5ad", self.h5ad_var, choose_file=True)
        row = self._path_field(basic_tab, row, "Output directory", self.outdir_var, choose_file=False)
        row = self._option_row(
            basic_tab,
            row,
            "Coordinates",
            widget=self._coords_dropdown(basic_tab),
            hint="auto | obsm:spatial | obs:centroid_x_y. obs centroid mode is converted to temporary obsm['spatial'] before KaroSpace export.",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Section groupby",
            widget=self._groupby_dropdown(basic_tab),
            hint="obs column used to split sections (same as load_spatial_data(groupby=...)).",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Initial color",
            widget=self._color_dropdown(basic_tab),
            hint="Initial viewer color (obs column or gene).",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Theme",
            widget=self._theme_dropdown(basic_tab),
            hint="KaroSpace viewer theme.",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Outline by",
            widget=self._outline_dropdown(basic_tab),
            hint="Optional metadata column for panel outlines.",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Viewer title",
            widget=self._entry(basic_tab, self.title_var),
            hint="Title shown in the exported HTML.",
        )
        downsample_container = ttk.Frame(basic_tab, style="Card.TFrame")
        ttk.Entry(downsample_container, textvariable=self.downsample_var, width=12).pack(side="left")
        ttk.Label(downsample_container, text="cells per section (blank = all)", style="Body.TLabel").pack(
            side="left", padx=(8, 0)
        )
        self._option_row(
            basic_tab,
            row,
            "Downsample",
            widget=downsample_container,
            hint="Maps to export_to_html(downsample=...).",
        )

        colors_tab.columnconfigure(0, weight=1)
        self.additional_colors_editor = SearchableListEditor(
            colors_tab,
            label="additional_colors (obs columns)",
            height=6,
            help_text="These become color options in KaroSpace. Use Inspect to load obs columns.",
        )
        self.additional_colors_editor.grid(row=0, column=0, sticky="ew", pady=(0, 10))

        self.groupby_editor = SearchableListEditor(
            colors_tab,
            label="groupby lists (marker/neighbor/interaction)",
            height=6,
            help_text="Used for marker_genes_groupby and optionally neighbor/interaction groupby.",
        )
        self.groupby_editor.grid(row=1, column=0, sticky="ew", pady=(0, 12))

        genes_card = ttk.Frame(colors_tab, style="Card.TFrame")
        genes_card.grid(row=2, column=0, sticky="ew")
        genes_card.columnconfigure(0, weight=1)
        ttk.Label(genes_card, text="Gene Selection", style="FieldLabel.TLabel").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._genes_mode_row(genes_card).grid(row=1, column=0, sticky="ew")
        ttk.Label(
            genes_card,
            text=(
                "hvgs:N maps to use_hvgs=True/hvg_limit=N. top_mean/list_file/manual_list map to explicit genes list "
                "with use_hvgs=False."
            ),
            style="Subheader.TLabel",
        ).grid(row=2, column=0, sticky="w", pady=(4, 8))
        self.manual_genes_editor = SearchableListEditor(
            genes_card,
            label="genes list (manual_list)",
            height=8,
            help_text="Search var_names and build the genes list with + Add / Remove.",
        )
        self.manual_genes_editor.grid(row=3, column=0, sticky="ew")

        selection_card = ttk.Frame(colors_tab, style="Card.TFrame")
        selection_card.grid(row=3, column=0, sticky="ew", pady=(12, 0))
        selection_card.columnconfigure(0, weight=1)
        self.selection_mode_check = ttk.Checkbutton(
            selection_card,
            text="Enable tick selection mode (from Inspect H5AD)",
            variable=self.selection_mode_var,
        )
        self.selection_mode_check.grid(row=0, column=0, sticky="w")
        ttk.Label(
            selection_card,
            text=(
                "When enabled, selected items below are used during export for additional colors, "
                "groupby lists, and manual genes."
            ),
            style="Subheader.TLabel",
        ).grid(row=1, column=0, sticky="w", pady=(4, 8))

        self.selection_mode_content = ttk.Frame(selection_card, style="Card.TFrame")
        self.selection_mode_content.grid(row=2, column=0, sticky="ew")
        self.selection_mode_content.columnconfigure(0, weight=1)
        self.selection_mode_content.columnconfigure(1, weight=1)

        self.selection_additional_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick additional_colors (obs)",
            height=7,
            help_text="Multi-select obs fields to include as additional viewer colors.",
        )
        self.selection_additional_picker.grid(row=0, column=0, sticky="ew", padx=(0, 8))

        self.selection_groupby_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick groupby lists (obs)",
            height=7,
            help_text="Multi-select obs columns for marker/neighbor/interaction groupby lists.",
        )
        self.selection_groupby_picker.grid(row=0, column=1, sticky="ew")

        self.selection_genes_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick genes (manual_list mode)",
            height=8,
            help_text="Multi-select genes from var_names. Used when genes mode is manual_list.",
        )
        self.selection_genes_picker.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(10, 0))

        selection_actions = ttk.Frame(self.selection_mode_content, style="Card.TFrame")
        selection_actions.grid(row=2, column=0, columnspan=2, sticky="w", pady=(10, 0))
        self.selection_apply_btn = ttk.Button(
            selection_actions,
            text="Apply selections to inputs",
            style="Secondary.TButton",
            command=self._apply_tick_selections_to_inputs,
        )
        self.selection_apply_btn.pack(side="left")
        self.selection_sync_btn = ttk.Button(
            selection_actions,
            text="Sync from current inputs",
            style="Secondary.TButton",
            command=self._load_inputs_into_tick_selection,
        )
        self.selection_sync_btn.pack(side="left", padx=(8, 0))

        advanced_tab.columnconfigure(0, weight=1)
        adv_header = ttk.Frame(advanced_tab, style="Card.TFrame")
        adv_header.grid(row=0, column=0, sticky="ew")
        self.advanced_toggle_btn = ttk.Button(
            adv_header,
            text="Show Advanced Options",
            style="Secondary.TButton",
            command=self._toggle_advanced,
        )
        self.advanced_toggle_btn.pack(anchor="w")
        ttk.Label(
            adv_header,
            text="Analytics and rendering parameters from KaroSpace export_to_html.",
            style="Subheader.TLabel",
        ).pack(anchor="w", pady=(6, 0))

        self.advanced_content = ttk.Frame(advanced_tab, style="Card.TFrame")
        self.advanced_content.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        self.advanced_content.columnconfigure(1, weight=1)

        min_panel_row = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Entry(min_panel_row, textvariable=self.min_panel_size_var, width=10).pack(side="left")
        ttk.Label(min_panel_row, text="px", style="Body.TLabel").pack(side="left", padx=(8, 0))
        self._option_row(
            self.advanced_content,
            0,
            "Min panel size",
            widget=min_panel_row,
            hint="Minimum section panel width in exported KaroSpace HTML.",
        )

        spot_row = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Combobox(spot_row, textvariable=self.spot_size_var, values=["auto", "adaptive", "density"], width=12).pack(
            side="left"
        )
        ttk.Label(spot_row, text="or numeric value", style="Body.TLabel").pack(side="left", padx=(8, 0))
        self._option_row(
            self.advanced_content,
            2,
            "Spot size",
            widget=spot_row,
            hint="Use auto/adaptive/density or a positive number.",
        )

        marker_row = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Entry(marker_row, textvariable=self.marker_genes_top_n_var, width=10).pack(side="left")
        ttk.Label(marker_row, text="top genes per group", style="Body.TLabel").pack(side="left", padx=(8, 0))
        self._option_row(
            self.advanced_content,
            4,
            "Marker genes top N",
            widget=marker_row,
            hint="Maps to marker_genes_top_n.",
        )

        neighbor_row = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Checkbutton(neighbor_row, text="Auto neighbor groupby = [initial color]", variable=self.neighbor_auto_var).pack(
            side="left"
        )
        ttk.Label(neighbor_row, text="Permutations", style="Body.TLabel").pack(side="left", padx=(12, 4))
        ttk.Entry(neighbor_row, textvariable=self.neighbor_permutations_var, width=8).pack(side="left")
        ttk.Label(neighbor_row, text="Seed", style="Body.TLabel").pack(side="left", padx=(12, 4))
        ttk.Entry(neighbor_row, textvariable=self.neighbor_stats_seed_var, width=8).pack(side="left")
        self._option_row(
            self.advanced_content,
            6,
            "Neighbor stats",
            widget=neighbor_row,
            hint="Permutations accepts integer or 'auto'.",
        )

        interaction_row_1 = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Checkbutton(
            interaction_row_1,
            text="Enable interaction markers (uses groupby list)",
            variable=self.interaction_markers_enabled_var,
        ).pack(side="left")
        self._option_row(self.advanced_content, 8, "Interaction markers", widget=interaction_row_1)

        interaction_row_2 = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Label(interaction_row_2, text="Top targets", style="Body.TLabel").pack(side="left")
        ttk.Entry(interaction_row_2, textvariable=self.interaction_markers_top_targets_var, width=8).pack(
            side="left", padx=(4, 10)
        )
        ttk.Label(interaction_row_2, text="Top genes", style="Body.TLabel").pack(side="left")
        ttk.Entry(interaction_row_2, textvariable=self.interaction_markers_top_genes_var, width=8).pack(
            side="left", padx=(4, 10)
        )
        ttk.Label(interaction_row_2, text="Min cells", style="Body.TLabel").pack(side="left")
        ttk.Entry(interaction_row_2, textvariable=self.interaction_markers_min_cells_var, width=8).pack(
            side="left", padx=(4, 10)
        )
        ttk.Label(interaction_row_2, text="Min neighbors", style="Body.TLabel").pack(side="left")
        ttk.Entry(interaction_row_2, textvariable=self.interaction_markers_min_neighbors_var, width=8).pack(
            side="left", padx=(4, 0)
        )
        self._option_row(
            self.advanced_content,
            9,
            "Interaction limits",
            widget=interaction_row_2,
            hint="Maps to interaction_markers_top_targets/top_genes/min_cells/min_neighbors.",
        )

        serve_row = ttk.Frame(self.advanced_content, style="Card.TFrame")
        ttk.Checkbutton(serve_row, text="Serve after export", variable=self.serve_var).pack(side="left")
        ttk.Label(serve_row, text="Port", style="Body.TLabel").pack(side="left", padx=(12, 4))
        ttk.Entry(serve_row, textvariable=self.port_var, width=10).pack(side="left")
        self._option_row(
            self.advanced_content,
            11,
            "Preview server",
            widget=serve_row,
            hint="Optional local server to open the latest generated KaroSpace_*.html file.",
        )
        self._set_advanced_visible(False)

        help_tab.columnconfigure(0, weight=1)
        help_tab.rowconfigure(0, weight=1)
        self.help_text = tk.Text(
            help_tab,
            wrap="word",
            background=palette["input_bg"],
            foreground=palette["text"],
            relief="solid",
            bd=1,
            highlightthickness=1,
            highlightbackground=palette["border"],
            highlightcolor=palette["accent"],
            insertbackground=palette["text"],
            padx=10,
            pady=10,
            height=18,
        )
        self.help_text.grid(row=0, column=0, sticky="nsew")
        help_scroll = ttk.Scrollbar(help_tab, orient="vertical", command=self.help_text.yview)
        help_scroll.grid(row=0, column=1, sticky="ns")
        self.help_text.configure(yscrollcommand=help_scroll.set)
        self.help_text.insert(
            "1.0",
            "Basic tab\n"
            "- Input .h5ad: absolute path to your AnnData file.\n"
            "- Output directory: folder where KaroSpace_YYYYMMDD_HHMMSS.html is written.\n"
            "- Coordinates: auto/obsm/obs-centroid modes are converted to KaroSpace spatial input.\n"
            "- Section groupby: section split column used by load_spatial_data.\n"
            "- Initial color/theme/outline/title map directly to export_to_html.\n"
            "- Downsample: integer cells per section (blank keeps all).\n\n"
            "Colors & Genes tab\n"
            "- additional_colors: obs columns offered as categorical coloring fields.\n"
            "- groupby lists: columns used for marker/neighbor/interaction analytics.\n"
            "- Tick selection mode: optional inspected multi-select mode. Tick obs/genes and export directly without manually adding rows.\n"
            "- genes mode:\n"
            "  hvgs -> use_hvgs=True with hvg_limit from count\n"
            "  top_mean -> compute top_mean genes list from adata\n"
            "  list_file -> choose a text file with one gene name per line\n"
            "  manual_list -> build list from var_names picker\n"
            "  (non-hvgs modes set use_hvgs=False)\n\n"
            "Advanced tab\n"
            "- Min panel size, spot size, marker top N.\n"
            "- Neighbor stats permutations/seed and auto groupby mode.\n"
            "- Interaction markers limits and enable/disable toggle.\n"
            "- Serve after export starts a local preview server for the latest KaroSpace_*.html export.\n\n"
            "Presets\n"
            "- Default: balanced defaults.\n"
            "\n"
            "Tip: click Inspect H5AD to load searchable dropdown choices from adata.obs and adata.var_names."
        )
        self.help_text.configure(state="disabled")

        button_row = ttk.Frame(controls, style="Card.TFrame")
        self.inspect_btn = ttk.Button(button_row, text="Inspect H5AD", style="Secondary.TButton", command=self._inspect_h5ad)
        self.inspect_btn.pack(side="left")

        self.export_btn = ttk.Button(button_row, text="Export", style="Primary.TButton", command=self._on_export)
        self.export_btn.pack(side="left", padx=(10, 0))

        self.stop_server_btn = ttk.Button(button_row, text="Stop Server", style="Secondary.TButton", command=self._stop_server)
        self.stop_server_btn.pack(side="left", padx=(10, 0))

        button_row.grid(row=4, column=0, columnspan=3, sticky="w", pady=(12, 0))

        ttk.Label(side, text="Runtime", style="Header.TLabel").grid(row=0, column=0, sticky="w")
        ttk.Label(side, textvariable=self.status_var, style="Subheader.TLabel").grid(row=1, column=0, sticky="w", pady=(2, 12))

        self.progress = ttk.Progressbar(side, mode="determinate", maximum=100, value=0, style="Good.Horizontal.TProgressbar")
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
            background=palette["input_bg"],
            foreground=palette["text"],
            insertbackground=palette["text"],
            relief="solid",
            bd=1,
            highlightthickness=1,
            highlightbackground=palette["border"],
            highlightcolor=palette["accent"],
            padx=10,
            pady=10,
            font=("Consolas", 10),
        )
        self.log_text.grid(row=0, column=0, sticky="nsew")
        log_scroll = ttk.Scrollbar(log_wrap, orient="vertical", command=self.log_text.yview)
        log_scroll.grid(row=0, column=1, sticky="ns")
        self.log_text.configure(yscrollcommand=log_scroll.set)
        self.log_text.configure(state="disabled")

        self.genes_mode_var.trace_add("write", lambda *_: self._update_genes_mode_visibility())
        self.neighbor_auto_var.trace_add("write", lambda *_: self._update_neighbor_groupby_state())
        self.selection_mode_var.trace_add("write", lambda *_: self._update_selection_mode_visibility())
        self.theme_var.trace_add("write", lambda *_: self._sync_app_theme_to_viewer_setting())
        self._apply_preset("default", log=False)
        self._sync_app_theme_to_viewer_setting()
        self._update_genes_mode_visibility()
        self._update_neighbor_groupby_state()
        self._set_selection_mode_visible(False)

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

    def _groupby_dropdown(self, parent: ttk.Frame) -> ttk.Combobox:
        self.groupby_combo = ttk.Combobox(parent, textvariable=self.section_groupby_var, state="normal")
        return self.groupby_combo

    def _color_dropdown(self, parent: ttk.Frame) -> ttk.Combobox:
        self.color_combo = ttk.Combobox(parent, textvariable=self.initial_color_var, state="normal")
        return self.color_combo

    def _outline_dropdown(self, parent: ttk.Frame) -> ttk.Combobox:
        self.outline_combo = ttk.Combobox(parent, textvariable=self.outline_by_var, state="normal")
        return self.outline_combo

    def _theme_dropdown(self, parent: ttk.Frame) -> ttk.Combobox:
        combo = ttk.Combobox(parent, textvariable=self.theme_var, values=["light", "dark"], state="readonly")
        return combo

    def _genes_mode_row(self, parent: ttk.Frame) -> ttk.Frame:
        wrap = ttk.Frame(parent, style="Card.TFrame")

        self.genes_mode_combo = ttk.Combobox(
            wrap,
            textvariable=self.genes_mode_var,
            values=["hvgs", "top_mean", "list_file", "manual_list"],
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
        if mode == "manual_list":
            self.genes_count_entry.pack_forget()
            self.genes_count_label.pack_forget()
            self.gene_list_entry.pack_forget()
            self.gene_list_button.pack_forget()
            if hasattr(self, "manual_genes_editor"):
                self.manual_genes_editor.grid()
            return

        if hasattr(self, "manual_genes_editor"):
            self.manual_genes_editor.grid_remove()

        if mode == "list_file":
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

    def _update_neighbor_groupby_state(self) -> None:
        # Groupby list is shared by marker/neighbor/interaction settings.
        # Keep it editable even when neighbor auto mode is enabled.
        self.groupby_editor.set_enabled(True)

    def _update_selection_mode_visibility(self) -> None:
        self._set_selection_mode_visible(bool(self.selection_mode_var.get()))

    def _set_selection_mode_visible(self, visible: bool) -> None:
        if visible:
            self.selection_mode_content.grid()
        else:
            self.selection_mode_content.grid_remove()

    def _load_inputs_into_tick_selection(self, *, log: bool = True) -> None:
        self.selection_additional_picker.set_selected(self.additional_colors_editor.get_items())
        self.selection_groupby_picker.set_selected(self.groupby_editor.get_items())
        self.selection_genes_picker.set_selected(self.manual_genes_editor.get_items())
        if log:
            self._log("Selection mode synced from current input lists.")

    def _apply_tick_selections_to_inputs(self) -> None:
        additional = self.selection_additional_picker.get_selected()
        groupby = self.selection_groupby_picker.get_selected()
        genes = self.selection_genes_picker.get_selected()

        self.additional_colors_editor.set_items(additional)
        self.groupby_editor.set_items(groupby)
        if genes:
            self.manual_genes_editor.set_items(genes)

        self.status_var.set("Selection picks applied")
        self._log(
            f"Applied selection mode picks: additional_colors={len(additional)}, groupby={len(groupby)}, genes={len(genes)}"
        )

    def _toggle_advanced(self) -> None:
        self._set_advanced_visible(not bool(self.advanced_open_var.get()))

    def _set_advanced_visible(self, visible: bool) -> None:
        self.advanced_open_var.set(bool(visible))
        if visible:
            self.advanced_content.grid()
            self.advanced_toggle_btn.configure(text="Hide Advanced Options")
        else:
            self.advanced_content.grid_remove()
            self.advanced_toggle_btn.configure(text="Show Advanced Options")

    def _on_canvas_configure(self, event) -> None:
        if hasattr(self, "_scroll_window"):
            self.scroll_canvas.itemconfigure(self._scroll_window, width=event.width)

    def _bind_mousewheel_scroll(self) -> None:
        self.bind_all("<MouseWheel>", self._on_mousewheel, add="+")
        self.bind_all("<Button-4>", self._on_mousewheel_linux_up, add="+")
        self.bind_all("<Button-5>", self._on_mousewheel_linux_down, add="+")

    def _on_mousewheel(self, event) -> None:
        if not hasattr(self, "scroll_canvas"):
            return
        if sys.platform == "darwin":
            delta = -1 * int(event.delta)
        else:
            delta = -1 * int(event.delta / 120) if event.delta else 0
        if delta == 0:
            return
        self.scroll_canvas.yview_scroll(delta, "units")

    def _on_mousewheel_linux_up(self, _event) -> None:
        if hasattr(self, "scroll_canvas"):
            self.scroll_canvas.yview_scroll(-1, "units")

    def _on_mousewheel_linux_down(self, _event) -> None:
        if hasattr(self, "scroll_canvas"):
            self.scroll_canvas.yview_scroll(1, "units")

    @staticmethod
    def _merge_unique(*groups: list[str]) -> list[str]:
        seen: set[str] = set()
        merged: list[str] = []
        for group in groups:
            for raw in group:
                value = str(raw).strip()
                if not value or value in seen:
                    continue
                seen.add(value)
                merged.append(value)
        return merged

    def _matches_inspected_h5ad(self, h5ad_path: Path) -> bool:
        inspected = self._inspected_h5ad_path
        if inspected is None:
            return False
        return inspected == h5ad_path.expanduser().resolve()

    def _apply_preset(self, name: str, *, log: bool = True) -> None:
        # Shared baseline.
        if not self.h5ad_var.get().strip():
            self.h5ad_var.set("/absolute/path/to/input.h5ad")
        if not self.outdir_var.get().strip():
            self.outdir_var.set(str((Path.cwd() / "karospace_export").resolve()))
        self.coords_var.set("auto")
        self.serve_var.set(False)
        self.port_var.set("8000")
        self.downsample_var.set("")
        self.genes_count_var.set("100")
        self.gene_list_path_var.set("")
        self.section_groupby_var.set("sample_id")
        self.initial_color_var.set("leiden")
        self.title_var.set("KaroSpace")
        if not self.theme_var.get().strip():
            self.theme_var.set("dark")
        self.outline_by_var.set("condition")
        self.min_panel_size_var.set("120")
        self.spot_size_var.set("auto")
        self.marker_genes_top_n_var.set("50")
        self.neighbor_auto_var.set(True)
        self.neighbor_permutations_var.set("auto")
        self.neighbor_stats_seed_var.set("0")
        self.interaction_markers_enabled_var.set(False)
        self.interaction_markers_top_targets_var.set("8")
        self.interaction_markers_top_genes_var.set("20")
        self.interaction_markers_min_cells_var.set("30")
        self.interaction_markers_min_neighbors_var.set("1")
        self.section_groupby_var.set("sample_id")
        self.initial_color_var.set("leiden")
        self.outline_by_var.set("condition")
        self.min_panel_size_var.set("150")
        self.additional_colors_editor.set_items(["leiden_1", "leiden_2", "gmm_mana_10"])
        self.groupby_editor.set_items([])
        self.genes_mode_var.set("hvgs")
        self.genes_count_var.set("20")
        self.manual_genes_editor.set_items(["Cd4", "Cd8a", "Gfap", "Mki67"])
        self.marker_genes_top_n_var.set("30")
        self.neighbor_auto_var.set(True)
        self.neighbor_permutations_var.set("auto")
        self.interaction_markers_enabled_var.set(False)
        self.status_var.set("Preset loaded: Default")
        label = "Default"

        self._update_genes_mode_visibility()
        self._update_neighbor_groupby_state()
        self._load_inputs_into_tick_selection(log=False)
        if log:
            self._log(f"Applied preset: {label}")

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

        path = Path(path_text).expanduser().resolve()
        if not path.exists():
            messagebox.showerror("Missing file", f"Input file not found:\n{path}")
            return

        self._log(f"Inspecting {path}")
        adata = None
        try:
            ad_mod = _get_anndata()
            try:
                adata = ad_mod.read_h5ad(path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(path)

            obs_cols = [str(c) for c in adata.obs.columns]
            obs_col_set = set(obs_cols)
            total_var_count = int(adata.n_vars)
            max_gene_choices = 120000
            var_name_set: set[str] | None = None
            if total_var_count > max_gene_choices:
                var_names = [str(g) for g in adata.var_names[:max_gene_choices]]
                self._log(
                    f"Large gene table detected ({total_var_count}). "
                    f"Loaded first {max_gene_choices} genes into pickers for responsiveness."
                )
            else:
                var_names = [str(g) for g in adata.var_names]
                var_name_set = set(var_names)

            self.additional_colors_editor.set_choices(obs_cols)
            self.groupby_editor.set_choices(obs_cols)
            self.manual_genes_editor.set_choices(var_names)
            self.selection_additional_picker.set_choices(obs_cols)
            self.selection_groupby_picker.set_choices(obs_cols)
            self.selection_genes_picker.set_choices(var_names)
            if hasattr(self, "groupby_combo"):
                self.groupby_combo.configure(values=obs_cols)
            if hasattr(self, "color_combo"):
                self.color_combo.configure(values=obs_cols + var_names[:200])
            if hasattr(self, "outline_combo"):
                self.outline_combo.configure(values=[""] + obs_cols)

            if self.section_groupby_var.get().strip() not in obs_col_set:
                if "sample_id" in obs_col_set:
                    self.section_groupby_var.set("sample_id")
                elif obs_cols:
                    self.section_groupby_var.set(obs_cols[0])
            initial_color = self.initial_color_var.get().strip()
            initial_in_obs = initial_color in obs_col_set
            initial_in_var = True if var_name_set is None else initial_color in var_name_set
            if not initial_in_obs and not initial_in_var:
                if "leiden" in obs_col_set:
                    self.initial_color_var.set("leiden")
                elif obs_cols:
                    self.initial_color_var.set(obs_cols[0])
                elif var_names:
                    self.initial_color_var.set(var_names[0])
            if self.outline_by_var.get().strip() and self.outline_by_var.get().strip() not in obs_col_set:
                if "condition" in obs_col_set:
                    self.outline_by_var.set("condition")
                else:
                    self.outline_by_var.set("")

            existing_additional = [name for name in self.additional_colors_editor.get_items() if name in obs_col_set]
            if not existing_additional:
                existing_additional = [c for c in obs_cols if c in {"cell_type", "leiden", "sample", "sample_id", "condition"}]
                if not existing_additional:
                    existing_additional = obs_cols[: min(4, len(obs_cols))]
            self.additional_colors_editor.set_items(existing_additional)

            existing_groupby = [name for name in self.groupby_editor.get_items() if name in obs_col_set]
            if not existing_groupby:
                existing_groupby = [c for c in obs_cols if c in {"sample_id", "sample", "condition", "batch", "donor"}]
                if not existing_groupby and obs_cols:
                    existing_groupby = [obs_cols[0]]
            self.groupby_editor.set_items(existing_groupby)

            if var_name_set is None:
                existing_genes = self.manual_genes_editor.get_items()
            else:
                existing_genes = [name for name in self.manual_genes_editor.get_items() if name in var_name_set]
            if not existing_genes:
                if var_name_set is None:
                    existing_genes = [g for g in ["Mki67", "Cd4", "Cd8a", "Gfap"] if g in var_names]
                else:
                    existing_genes = [g for g in ["Mki67", "Cd4", "Cd8a", "Gfap"] if g in var_name_set]
                if not existing_genes:
                    existing_genes = var_names[: min(10, len(var_names))]
            self.manual_genes_editor.set_items(existing_genes)
            self._load_inputs_into_tick_selection(log=False)

            if obs_cols:
                self._log(f"Loaded {len(obs_cols)} obs columns into additional_colors/groupby pickers.")
            if var_names:
                self._log(f"Loaded {len(var_names)} genes into manual genes picker.")

            has_spatial = "spatial" in adata.obsm
            has_centroid = {"centroid_x", "centroid_y"}.issubset(set(obs_cols))
            if has_spatial:
                self.coords_var.set("obsm:spatial")
                inspected_coords_mode = "obsm:spatial"
            elif has_centroid:
                self.coords_var.set("obs:centroid_x_y")
                inspected_coords_mode = "obs:centroid_x_y"
            else:
                self.coords_var.set("auto")
                inspected_coords_mode = None

            self._inspected_h5ad_path = path
            self._inspected_var_name_set = var_name_set
            self._inspected_coords_mode = inspected_coords_mode

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

    @staticmethod
    def _parse_positive_int(label: str, raw: str) -> int:
        text = str(raw).strip()
        try:
            value = int(text)
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer.") from exc
        if value <= 0:
            raise ValueError(f"{label} must be > 0.")
        return value

    @staticmethod
    def _parse_non_negative_int(label: str, raw: str) -> int:
        text = str(raw).strip()
        try:
            value = int(text)
        except ValueError as exc:
            raise ValueError(f"{label} must be an integer.") from exc
        if value < 0:
            raise ValueError(f"{label} must be >= 0.")
        return value

    @staticmethod
    def _parse_spot_size(raw: str) -> float | str | None:
        text = str(raw).strip()
        if not text:
            return "auto"
        if text.lower() in {"auto", "adaptive", "density"}:
            return "auto"
        try:
            value = float(text)
        except ValueError as exc:
            raise ValueError("Spot size must be auto/adaptive/density or a positive number.") from exc
        if value <= 0:
            raise ValueError("Spot size must be > 0.")
        return value

    @staticmethod
    def _parse_neighbor_permutations(raw: str) -> int | None:
        text = str(raw).strip().lower()
        if not text or text == "auto":
            return None
        try:
            value = int(text)
        except ValueError as exc:
            raise ValueError("Neighbor permutations must be an integer or 'auto'.") from exc
        if value < 0:
            raise ValueError("Neighbor permutations must be >= 0.")
        return value

    def _load_var_names(self, h5ad_path: Path) -> set[str]:
        if self._matches_inspected_h5ad(h5ad_path) and self._inspected_var_name_set is not None:
            return set(self._inspected_var_name_set)

        ad_mod = _get_anndata()
        adata = None
        try:
            try:
                adata = ad_mod.read_h5ad(h5ad_path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(h5ad_path)
            return {str(v) for v in adata.var_names}
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    def _compute_top_mean_genes(self, h5ad_path: Path, count: int) -> list[str]:
        ad_mod = _get_anndata()
        np_mod = _get_numpy()
        adata = None
        try:
            try:
                adata = ad_mod.read_h5ad(h5ad_path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(h5ad_path)
            mode = parse_genes_mode(f"top_mean:{count}")
            obs_idx = np_mod.arange(adata.n_obs)
            return resolve_gene_names(adata, mode, obs_indices=obs_idx)
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    def _resolve_gene_settings(
        self,
        h5ad_path: Path,
        mode: str,
        *,
        manual_genes_override: list[str] | None = None,
    ) -> tuple[list[str] | None, bool, int]:
        normalized = mode.strip().lower()
        if normalized == "hvgs":
            hvg_limit = self._parse_positive_int("HVG limit", self.genes_count_var.get())
            return None, True, hvg_limit

        use_hvgs = False
        hvg_limit = self._parse_positive_int("HVG limit", self.genes_count_var.get() or "20")
        if normalized == "top_mean":
            count = self._parse_positive_int("Top-mean genes", self.genes_count_var.get())
            genes = self._compute_top_mean_genes(h5ad_path, count)
            return genes, use_hvgs, hvg_limit

        if normalized == "list_file":
            list_path_text = self.gene_list_path_var.get().strip()
            if not list_path_text:
                raise ValueError("Choose a gene list file for list_file mode.")
            list_path = Path(list_path_text).expanduser()
            if not list_path.exists():
                raise ValueError(f"Gene list file not found: {list_path}")
            genes = [line.strip() for line in list_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        elif normalized == "manual_list":
            genes = manual_genes_override if manual_genes_override is not None else self.manual_genes_editor.get_items()
        else:
            raise ValueError("Genes mode must be hvgs, top_mean, list_file, or manual_list.")

        genes = self._merge_unique(genes)
        if not genes:
            raise ValueError("Add at least one gene for list/manual genes mode.")

        var_names = self._load_var_names(h5ad_path)
        missing = [gene for gene in genes if gene not in var_names]
        if missing:
            preview = ", ".join(missing[:10])
            raise ValueError(f"{len(missing)} genes are missing in var_names: {preview}")
        return genes, use_hvgs, hvg_limit

    def _parse_config(self) -> BuilderConfig:
        h5ad_text = self.h5ad_var.get().strip()
        outdir_text = self.outdir_var.get().strip()

        if not h5ad_text:
            raise ValueError("Input .h5ad is required.")
        if not outdir_text:
            raise ValueError("Output directory is required.")

        h5ad_path = Path(h5ad_text).expanduser().resolve()
        outdir = Path(outdir_text).expanduser().resolve()
        if not h5ad_path.exists():
            raise ValueError(f"Input .h5ad not found: {h5ad_path}")

        outdir.mkdir(parents=True, exist_ok=True)

        coords_raw = self.coords_var.get().strip().lower() or "auto"
        if coords_raw not in {"auto", "obsm:spatial", "obs:centroid_x_y"}:
            raise ValueError("Coordinates must be auto, obsm:spatial, or obs:centroid_x_y.")
        coords_mode = None if coords_raw == "auto" else coords_raw

        section_groupby = self.section_groupby_var.get().strip()
        if not section_groupby:
            raise ValueError("Section groupby is required.")
        initial_color = self.initial_color_var.get().strip()
        if not initial_color:
            raise ValueError("Initial color is required.")
        theme = self.theme_var.get().strip().lower() or "dark"
        if theme not in {"light", "dark"}:
            raise ValueError("Theme must be 'light' or 'dark'.")
        title = self.title_var.get().strip() or "KaroSpace"
        outline_by = self.outline_by_var.get().strip() or None
        min_panel_size = self._parse_positive_int("Min panel size", self.min_panel_size_var.get())
        spot_size = self._parse_spot_size(self.spot_size_var.get())

        downsample_text = self.downsample_var.get().strip()
        downsample = None
        if downsample_text:
            downsample = self._parse_positive_int("Downsample", downsample_text)

        additional_colors = self._merge_unique(self.additional_colors_editor.get_items())
        groupby_lists = self._merge_unique(self.groupby_editor.get_items())
        selection_mode = bool(self.selection_mode_var.get())
        manual_genes_override: list[str] | None = None

        if selection_mode:
            selected_additional = self.selection_additional_picker.get_selected()
            selected_groupby = self.selection_groupby_picker.get_selected()
            if selected_additional:
                additional_colors = self._merge_unique(selected_additional)
            if selected_groupby:
                groupby_lists = self._merge_unique(selected_groupby)
            if self.genes_mode_var.get().strip().lower() == "manual_list":
                selected_genes = self.selection_genes_picker.get_selected()
                if selected_genes:
                    manual_genes_override = self._merge_unique(selected_genes)

        mode = self.genes_mode_var.get().strip().lower()
        genes, use_hvgs, hvg_limit = self._resolve_gene_settings(
            h5ad_path,
            mode,
            manual_genes_override=manual_genes_override,
        )

        marker_genes_top_n = self._parse_positive_int("Marker genes top N", self.marker_genes_top_n_var.get())
        neighbor_permutations = self._parse_neighbor_permutations(self.neighbor_permutations_var.get())
        neighbor_seed = self._parse_non_negative_int("Neighbor stats seed", self.neighbor_stats_seed_var.get() or "0")
        marker_genes_groupby = groupby_lists or None
        if bool(self.neighbor_auto_var.get()):
            neighbor_stats_groupby = [initial_color]
        else:
            neighbor_stats_groupby = groupby_lists or None
        interaction_enabled = bool(self.interaction_markers_enabled_var.get())
        interaction_groupby = (groupby_lists or None) if interaction_enabled else None
        interaction_top_targets = self._parse_positive_int(
            "Interaction top targets", self.interaction_markers_top_targets_var.get()
        )
        interaction_top_genes = self._parse_positive_int(
            "Interaction top genes", self.interaction_markers_top_genes_var.get()
        )
        interaction_min_cells = self._parse_positive_int(
            "Interaction min cells", self.interaction_markers_min_cells_var.get()
        )
        interaction_min_neighbors = self._parse_positive_int(
            "Interaction min neighbors", self.interaction_markers_min_neighbors_var.get()
        )

        return BuilderConfig(
            h5ad_path=h5ad_path,
            outdir=outdir,
            coords_mode=coords_mode,
            section_groupby=section_groupby,
            initial_color=initial_color,
            title=title,
            theme=theme,
            outline_by=outline_by,
            min_panel_size=min_panel_size,
            spot_size=spot_size,
            downsample=downsample,
            additional_colors=additional_colors or None,
            genes=genes,
            use_hvgs=use_hvgs,
            hvg_limit=hvg_limit,
            marker_genes_groupby=marker_genes_groupby,
            genes_mode=mode,
            marker_genes_top_n=marker_genes_top_n,
            neighbor_stats_groupby=neighbor_stats_groupby,
            neighbor_stats_permutations=neighbor_permutations,
            neighbor_stats_seed=neighbor_seed,
            interaction_markers_enabled=interaction_enabled,
            interaction_markers_groupby=interaction_groupby,
            interaction_markers_top_targets=interaction_top_targets,
            interaction_markers_top_genes=interaction_top_genes,
            interaction_markers_min_cells=interaction_min_cells,
            interaction_markers_min_neighbors=interaction_min_neighbors,
        )

    @staticmethod
    def _import_karospace_api():
        import importlib
        import importlib.metadata as importlib_metadata
        import inspect

        # scanpy imports numba with cache=True in some code paths; in certain
        # environments this crashes during import. Disable JIT by default for
        # robust GUI startup/export unless the user explicitly overrides it.
        os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

        def _install_scanpy_inspect_fallback() -> None:
            original_getsource = getattr(inspect, "getsource", None)
            if original_getsource is None:
                return
            if getattr(original_getsource, "_ksb_scanpy_fallback", False):
                return

            def _obj_module_name(obj) -> str:
                module_name = str(getattr(obj, "__module__", "") or "")
                if module_name:
                    return module_name
                if inspect.ismodule(obj):
                    return str(getattr(obj, "__name__", "") or "")
                obj_cls = getattr(obj, "__class__", None)
                if obj_cls is None:
                    return ""
                return str(getattr(obj_cls, "__module__", "") or "")

            def _ksb_safe_getsource(obj) -> str:
                try:
                    return original_getsource(obj)
                except OSError:
                    module_name = _obj_module_name(obj)
                    if module_name.startswith("scanpy"):
                        return ""
                    raise

            setattr(_ksb_safe_getsource, "_ksb_scanpy_fallback", True)
            inspect.getsource = _ksb_safe_getsource

        def _install_metadata_version_fallback() -> None:
            original_version = getattr(importlib_metadata, "version", None)
            if original_version is None:
                return
            if getattr(original_version, "_ksb_scikit_fallback", False):
                return

            def _ksb_version_with_fallback(distribution_name: str) -> str:
                try:
                    return original_version(distribution_name)
                except importlib_metadata.PackageNotFoundError as missing_error:
                    normalized = str(distribution_name).replace("_", "-").lower()
                    if normalized not in {"scikit-learn", "sklearn"}:
                        raise
                    for alias in ("scikit-learn", "scikit_learn", "sklearn"):
                        if alias == distribution_name:
                            continue
                        try:
                            return original_version(alias)
                        except importlib_metadata.PackageNotFoundError:
                            continue
                    try:
                        import sklearn  # type: ignore

                        version = getattr(sklearn, "__version__", "")
                        if version:
                            return str(version)
                    except Exception:
                        pass
                    raise missing_error

            setattr(_ksb_version_with_fallback, "_ksb_scikit_fallback", True)
            importlib_metadata.version = _ksb_version_with_fallback

        def _missing_module_name(error: BaseException) -> str | None:
            current: BaseException | None = error
            while current is not None:
                if isinstance(current, ModuleNotFoundError):
                    return getattr(current, "name", None)
                current = current.__cause__
            return None

        def _missing_metadata_name(error: BaseException) -> str | None:
            current: BaseException | None = error
            while current is not None:
                if isinstance(current, importlib_metadata.PackageNotFoundError):
                    name = getattr(current, "name", None)
                    if name:
                        return str(name)
                current = current.__cause__
            return None

        def _has_scanpy_source_error(error: BaseException) -> bool:
            current: BaseException | None = error
            while current is not None:
                if isinstance(current, OSError) and "could not get source code" in str(current).lower():
                    return True
                current = current.__cause__
            return False

        def _raise_dependency_error(error: BaseException) -> None:
            missing = _missing_module_name(error)
            if missing:
                raise RuntimeError(
                    f"Missing dependency '{missing}' required by KaroSpace. "
                    "Install dependencies in this environment (for example: "
                    "pip install scanpy or pip install -e /path/to/spatial-viewer)."
                ) from error
            missing_metadata = _missing_metadata_name(error)
            if missing_metadata:
                raise RuntimeError(
                    f"Missing package metadata for '{missing_metadata}' required by KaroSpace/scanpy. "
                    "If you are using a desktop binary, rebuild it with updated PyInstaller metadata bundling. "
                    "For source installs, reinstall dependencies in the active environment."
                ) from error
            message = str(error)
            if "cannot cache function" in message and "numba" in message.lower():
                raise RuntimeError(
                    "scanpy/numba failed during import cache initialization. "
                    "Set NUMBA_DISABLE_JIT=1 and restart KaroSpaceBuilder."
                ) from error
            if _has_scanpy_source_error(error):
                raise RuntimeError(
                    "scanpy import failed while inspecting plotting source code in a frozen app. "
                    "Rebuild KaroSpaceBuilder with updated packaging, then retry export."
                ) from error

        _install_scanpy_inspect_fallback()
        _install_metadata_version_fallback()

        try:
            module = importlib.import_module("karospace")
            return module.load_spatial_data, module.export_to_html
        except Exception as exc:
            root_exc = exc
            candidates = [
                (Path(__file__).resolve().parents[2] / "spatial-viewer").resolve(),
                (Path.cwd().parent / "spatial-viewer").resolve(),
            ]
            for candidate in candidates:
                package_dir = candidate / "karospace"
                if not package_dir.exists():
                    continue
                candidate_str = str(candidate)
                if candidate_str not in sys.path:
                    sys.path.insert(0, candidate_str)
                try:
                    module = importlib.import_module("karospace")
                    return module.load_spatial_data, module.export_to_html
                except Exception as inner_exc:
                    _raise_dependency_error(inner_exc)
                    continue
            _raise_dependency_error(root_exc)
            raise RuntimeError(
                "Could not import 'karospace'. Install it in this environment before exporting "
                "(for example: pip install -e /path/to/spatial-viewer)."
            ) from root_exc

    @staticmethod
    def _detect_coords_mode(h5ad_path: Path) -> str:
        ad_mod = _get_anndata()
        adata = None
        try:
            try:
                adata = ad_mod.read_h5ad(h5ad_path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(h5ad_path)
            if "spatial" in adata.obsm:
                return "obsm:spatial"
            obs_cols = set(str(c) for c in adata.obs.columns)
            if {"centroid_x", "centroid_y"}.issubset(obs_cols):
                return "obs:centroid_x_y"
            raise ValueError(
                "Could not detect coordinates. Add adata.obsm['spatial'] or obs columns centroid_x/centroid_y."
            )
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    @staticmethod
    def _build_centroid_spatial_h5ad(h5ad_path: Path) -> Path:
        ad_mod = _get_anndata()
        np_mod = _get_numpy()
        adata = ad_mod.read_h5ad(h5ad_path)
        if "centroid_x" not in adata.obs.columns or "centroid_y" not in adata.obs.columns:
            raise ValueError("coords=obs:centroid_x_y requires obs columns centroid_x and centroid_y.")
        coords = adata.obs[["centroid_x", "centroid_y"]].to_numpy(dtype=np_mod.float32)
        adata.obsm["spatial"] = coords
        with tempfile.NamedTemporaryFile(suffix=".h5ad", prefix="karospace_builder_coords_", delete=False) as handle:
            temp_path = Path(handle.name)
        adata.write_h5ad(temp_path)
        return temp_path

    def _resolve_export_input(self, config: BuilderConfig) -> tuple[Path, str, Path | None]:
        mode = config.coords_mode
        if mode is None and self._matches_inspected_h5ad(config.h5ad_path) and self._inspected_coords_mode:
            mode = self._inspected_coords_mode
        if mode is None:
            mode = self._detect_coords_mode(config.h5ad_path)
        if mode == "obsm:spatial":
            return config.h5ad_path, "spatial", None
        if mode == "obs:centroid_x_y":
            temp_h5ad = self._build_centroid_spatial_h5ad(config.h5ad_path)
            return temp_h5ad, "spatial", temp_h5ad
        raise ValueError(f"Unsupported coordinates mode: {mode}")

    @staticmethod
    def _format_bytes(num_bytes: int) -> str:
        if num_bytes < 1024:
            return f"{num_bytes} B"
        size = float(num_bytes)
        for unit in ("KiB", "MiB", "GiB", "TiB"):
            size /= 1024.0
            if size < 1024.0 or unit == "TiB":
                return f"{size:.1f} {unit}"
        return f"{num_bytes} B"

    @staticmethod
    def _dedupe_names(values: list[str] | None) -> list[str]:
        seen: set[str] = set()
        out: list[str] = []
        for raw in values or []:
            name = str(raw).strip()
            if not name or name in seen:
                continue
            seen.add(name)
            out.append(name)
        return out

    def _sanitize_analytics_groupbys(
        self,
        dataset,
        *,
        marker_genes_groupby: list[str] | None,
        neighbor_stats_groupby: list[str] | None,
        interaction_markers_groupby: list[str] | None,
        neighbor_stats_permutations: int | None,
    ) -> tuple[list[str] | None, list[str], list[str] | None, list[str]]:
        warnings: list[str] = []
        marker = self._dedupe_names(marker_genes_groupby)
        neighbor = self._dedupe_names(neighbor_stats_groupby)
        interaction = self._dedupe_names(interaction_markers_groupby)

        adata = getattr(dataset, "adata", None)
        obs = getattr(adata, "obs", None)
        if obs is None:
            if neighbor_stats_groupby is None:
                warnings.append(
                    "Neighbor stats source columns could not be inspected; auto-neighbor fallback is disabled for safety."
                )
            return marker or None, neighbor, interaction or None, warnings

        obs_cols = {str(c) for c in getattr(obs, "columns", [])}
        unique_cache: dict[str, int | None] = {}

        def unique_count(column: str) -> int | None:
            if column in unique_cache:
                return unique_cache[column]
            try:
                count = int(obs[column].nunique(dropna=True))
            except Exception:
                count = None
            unique_cache[column] = count
            return count

        def drop_missing(columns: list[str], label: str) -> list[str]:
            kept: list[str] = []
            for column in columns:
                if column not in obs_cols:
                    warnings.append(f"Skipping {label} '{column}': not found in adata.obs.")
                    continue
                kept.append(column)
            return kept

        marker = drop_missing(marker, "marker groupby")
        neighbor = drop_missing(neighbor, "neighbor groupby")
        interaction = drop_missing(interaction, "interaction groupby")

        marker_kept: list[str] = []
        for column in marker:
            count = unique_count(column)
            if count is not None and count < 2:
                warnings.append(f"Skipping marker groupby '{column}': only {count} unique value.")
                continue
            if count is not None and count > self._MARKER_GROUPBY_MAX_UNIQUE:
                warnings.append(
                    f"Skipping marker groupby '{column}': {count:,} unique values exceeds "
                    f"{self._MARKER_GROUPBY_MAX_UNIQUE:,}."
                )
                continue
            marker_kept.append(column)

        interaction_kept: list[str] = []
        for column in interaction:
            count = unique_count(column)
            if count is not None and count < 2:
                warnings.append(f"Skipping interaction groupby '{column}': only {count} unique value.")
                continue
            if count is not None and count > self._INTERACTION_GROUPBY_MAX_UNIQUE:
                warnings.append(
                    f"Skipping interaction groupby '{column}': {count:,} unique values exceeds "
                    f"{self._INTERACTION_GROUPBY_MAX_UNIQUE:,}."
                )
                continue
            interaction_kept.append(column)

        permutation_factor = 3 if int(neighbor_stats_permutations or 0) > 0 else 1
        neighbor_kept: list[str] = []
        for column in neighbor:
            count = unique_count(column)
            if count is not None and count < 2:
                warnings.append(f"Skipping neighbor groupby '{column}': only {count} unique value.")
                continue
            if count is not None:
                dense_bytes = 8 * count * count * permutation_factor
                if dense_bytes > self._NEIGHBOR_MATRIX_BUDGET_BYTES:
                    warnings.append(
                        f"Skipping neighbor groupby '{column}': {count:,} unique values would require ~"
                        f"{self._format_bytes(dense_bytes)} dense memory."
                    )
                    continue
            neighbor_kept.append(column)

        if neighbor_stats_groupby is None:
            warnings.append("Neighbor stats auto-fallback disabled in builder; add explicit groupby columns if needed.")

        return marker_kept or None, neighbor_kept, interaction_kept or None, warnings

    def _set_busy(self, busy: bool) -> None:
        widgets = [
            self.export_btn,
            self.inspect_btn,
            self.genes_mode_combo,
            self.genes_count_entry,
            self.gene_list_entry,
            self.gene_list_button,
            self.advanced_toggle_btn,
            self.selection_mode_check,
            self.selection_apply_btn,
            self.selection_sync_btn,
        ]
        for widget in widgets:
            if busy:
                widget.state(["disabled"])
            else:
                widget.state(["!disabled"])
        self.additional_colors_editor.set_enabled(not busy)
        self.groupby_editor.set_enabled(not busy)
        self.manual_genes_editor.set_enabled(not busy)
        self.selection_additional_picker.set_enabled(not busy)
        self.selection_groupby_picker.set_enabled(not busy)
        self.selection_genes_picker.set_enabled(not busy)

        if busy:
            self._set_progress(0, "Queued")
        else:
            if self.status_var.get().startswith("Export running..."):
                self.status_var.set("Ready")

    @staticmethod
    def _coerce_progress_value(value: object) -> int:
        try:
            percent = int(round(float(value)))
        except (TypeError, ValueError):
            percent = 0
        return max(0, min(100, percent))

    def _set_progress(self, value: object, stage: str | None = None) -> None:
        percent = self._coerce_progress_value(value)
        self.progress.configure(value=percent)
        if stage:
            self.status_var.set(f"Export running... {percent}% | {stage}")

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
        self._log(
            "Export options: "
            f"coords={config.coords_mode or 'auto'}, "
            f"groupby={config.section_groupby}, "
            f"initial_color={config.initial_color}, "
            f"theme={config.theme}, "
            f"downsample={config.downsample if config.downsample is not None else 'all'}."
        )
        self._log(
            "Gene settings: "
            f"mode={config.genes_mode}, "
            f"use_hvgs={config.use_hvgs}, "
            f"hvg_limit={config.hvg_limit}, "
            f"manual_genes={len(config.genes or [])}."
        )
        self._log(
            "Analytics settings: "
            f"marker_groupby={len(config.marker_genes_groupby or [])}, "
            f"neighbor_groupby={len(config.neighbor_stats_groupby or [])}, "
            f"neighbor_permutations={config.neighbor_stats_permutations if config.neighbor_stats_permutations is not None else 'auto'}, "
            f"interaction_enabled={config.interaction_markers_enabled}."
        )
        serve_after_export = bool(self.serve_var.get())
        thread = threading.Thread(target=self._run_export, args=(config, serve_after_export), daemon=True)
        self._export_thread = thread
        thread.start()

    def _run_export(self, config: BuilderConfig, serve_after_export: bool) -> None:
        temp_h5ad: Path | None = None
        total_started = time.perf_counter()

        def emit_progress(percent: int, stage: str, detail: str | None = None) -> None:
            self._queue.put(("progress", (percent, stage, detail)))

        try:
            emit_progress(5, "Importing API", "Resolving karospace export functions.")
            load_spatial_data, export_to_html = self._import_karospace_api()

            emit_progress(12, "Preparing input", "Resolving coordinates mode and source data.")
            input_path, spatial_key, temp_h5ad = self._resolve_export_input(config)
            if temp_h5ad is not None:
                self._queue.put(("log", "Converted obs centroid_x/centroid_y to temporary obsm['spatial']."))

            emit_progress(
                25,
                "Loading spatial data",
                f"Reading {input_path} with groupby='{config.section_groupby}' and spatial_key='{spatial_key}'.",
            )
            load_started = time.perf_counter()
            dataset = load_spatial_data(
                str(input_path),
                groupby=config.section_groupby,
                spatial_key=spatial_key,
            )
            load_elapsed = time.perf_counter() - load_started
            emit_progress(
                55,
                "Validating analytics",
                f"Dataset loaded: sections={int(dataset.n_sections)}, cells={int(dataset.n_cells)} in {load_elapsed:.1f}s.",
            )

            marker_groupby, neighbor_groupby, interaction_groupby, guard_warnings = self._sanitize_analytics_groupbys(
                dataset,
                marker_genes_groupby=config.marker_genes_groupby,
                neighbor_stats_groupby=config.neighbor_stats_groupby,
                interaction_markers_groupby=(
                    config.interaction_markers_groupby if config.interaction_markers_enabled else None
                ),
                neighbor_stats_permutations=config.neighbor_stats_permutations,
            )
            for warning in guard_warnings:
                self._queue.put(("log", warning))

            neighbor_permutations = config.neighbor_stats_permutations if neighbor_groupby else 0
            emit_progress(
                68,
                "Preparing output",
                "Resolved analytics groupby columns: "
                f"marker={len(marker_groupby or [])}, "
                f"neighbor={len(neighbor_groupby)}, "
                f"interaction={len(interaction_groupby or []) if config.interaction_markers_enabled else 0}, "
                f"neighbor_permutations={neighbor_permutations if neighbor_permutations else 'off'}.",
            )

            output_html = self._build_output_html_path(config.outdir)
            emit_progress(76, "Writing viewer", f"Exporting HTML viewer to {output_html}.")
            export_started = time.perf_counter()
            output_html_path = Path(
                export_to_html(
                    dataset,
                    output_path=str(output_html),
                    color=config.initial_color,
                    title=config.title,
                    min_panel_size=config.min_panel_size,
                    spot_size=config.spot_size,
                    downsample=config.downsample,
                    theme=config.theme,
                    outline_by=config.outline_by,
                    additional_colors=config.additional_colors,
                    genes=config.genes,
                    use_hvgs=config.use_hvgs,
                    hvg_limit=config.hvg_limit,
                    marker_genes_groupby=marker_groupby,
                    marker_genes_top_n=config.marker_genes_top_n,
                    neighbor_stats_groupby=neighbor_groupby,
                    neighbor_stats_permutations=neighbor_permutations,
                    neighbor_stats_seed=config.neighbor_stats_seed,
                    interaction_markers_groupby=interaction_groupby if config.interaction_markers_enabled else None,
                    interaction_markers_top_targets=config.interaction_markers_top_targets,
                    interaction_markers_top_genes=config.interaction_markers_top_genes,
                    interaction_markers_min_cells=config.interaction_markers_min_cells,
                    interaction_markers_min_neighbors=config.interaction_markers_min_neighbors,
                )
            ).expanduser()
            export_elapsed = time.perf_counter() - export_started
            emit_progress(95, "Finalizing", f"Viewer bundle created in {export_elapsed:.1f}s: {output_html_path}")

            result = AppResult(
                outdir=output_html_path.parent,
                n_cells=int(dataset.n_cells),
                n_sections=int(dataset.n_sections),
                output_html=output_html_path,
            )
            total_elapsed = time.perf_counter() - total_started
            emit_progress(100, "Complete", f"Total export time: {total_elapsed:.1f}s.")
            self._queue.put(("done", result))
            if serve_after_export:
                self._queue.put(("log", "Starting preview server for the exported viewer."))
                self._queue.put(("start_server", output_html_path))
        except Exception:
            self._queue.put(("error", traceback.format_exc()))
        finally:
            if temp_h5ad is not None:
                try:
                    temp_h5ad.unlink(missing_ok=True)
                except Exception:
                    pass

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
                self._last_output_html = result.output_html
                self._log(
                    f"Export complete. sections={result.n_sections}, cells={result.n_cells}, html={result.output_html}"
                )
                self.status_var.set("Export complete")
            elif kind == "start_server":
                outdir = payload
                assert isinstance(outdir, Path)
                self._start_server(outdir)
            elif kind == "progress":
                if not isinstance(payload, tuple) or len(payload) < 2:
                    continue
                percent, stage = payload[0], str(payload[1]).strip()
                detail = str(payload[2]).strip() if len(payload) >= 3 and payload[2] is not None else ""
                self._set_progress(percent, stage or None)
                if detail:
                    self._log(detail)
            elif kind == "log":
                self._log(str(payload))
            elif kind == "error":
                self._set_busy(False)
                details = str(payload)
                self._log("Export failed. See traceback in popup.")
                self.status_var.set("Export failed")
                messagebox.showerror("Export failed", details)

        self.after(120, self._poll_events)

    def _build_output_html_path(self, outdir: Path) -> Path:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = f"{self._OUTPUT_HTML_BASENAME}_{stamp}"
        candidate = outdir / f"{base}.html"
        counter = 1
        while candidate.exists():
            candidate = outdir / f"{base}_{counter:02d}.html"
            counter += 1
        return candidate

    def _resolve_viewer_html(self, outdir: Path) -> Path:
        pattern = f"{self._OUTPUT_HTML_BASENAME}_*.html"
        candidates = sorted(outdir.glob(pattern))
        if candidates:
            return candidates[-1]
        legacy_named = outdir / f"{self._OUTPUT_HTML_BASENAME}.html"
        if legacy_named.exists():
            return legacy_named
        legacy = outdir / "index.html"
        if legacy.exists():
            return legacy
        return outdir / pattern.replace("*", "YYYYMMDD_HHMMSS")

    def _resolve_viewer_from_state(self) -> Path | None:
        if self._last_output_html is not None and self._last_output_html.exists():
            return self._last_output_html
        outdir = self._last_outdir or (Path(self.outdir_var.get().strip()).expanduser() if self.outdir_var.get().strip() else None)
        if outdir is None:
            return None
        return self._resolve_viewer_html(outdir)

    def _start_server(self, target_path: Path | None = None) -> None:
        if target_path is None:
            viewer_file = self._resolve_viewer_from_state()
        else:
            viewer_file = target_path if target_path.suffix.lower() == ".html" else self._resolve_viewer_html(target_path)

        if viewer_file is None:
            messagebox.showinfo("No export", "Run an export first.")
            return
        target = viewer_file.parent

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
        self._log(f"Serving {target} at http://127.0.0.1:{port}/{viewer_file.name}")
        self.status_var.set(f"Serving on :{port}")
        webbrowser.open_new_tab(f"http://127.0.0.1:{port}/{viewer_file.name}")

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
        viewer = self._resolve_viewer_from_state()
        if viewer is None:
            messagebox.showinfo("No output", "Run an export first.")
            return

        if self._server is not None:
            try:
                port = int(self.port_var.get().strip())
            except ValueError:
                port = 8000
            webbrowser.open_new_tab(f"http://127.0.0.1:{port}/{viewer.name}")
            return

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
        self.unbind_all("<MouseWheel>")
        self.unbind_all("<Button-4>")
        self.unbind_all("<Button-5>")
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
