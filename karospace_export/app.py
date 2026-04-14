from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from functools import partial
from pathlib import Path
import http.server
import json
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
    from tkinter import filedialog, messagebox
    import customtkinter as ctk
except Exception as exc:  # pragma: no cover - platform/runtime guard
    tk = None
    filedialog = None
    messagebox = None
    ctk = None
    TK_IMPORT_ERROR = exc
else:
    TK_IMPORT_ERROR = None

_CTK_FRAME_BASE = ctk.CTkFrame if ctk is not None else object

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
    "background": "#f0f1f4",
    "text": "#1a1d23",
    "panel_bg": "#ffffff",
    "border": "#d8dbe1",
    "input_bg": "#f8f9fb",
    "muted": "#6b7280",
    "hover_bg": "#e6e8ed",
    "accent": _KI_COLORS["plum"],
    "accent_strong": _KI_COLORS["plum_dark"],
    "on_accent": "#ffffff",
    "hero_bg": "#faf7f9",
    "section_text": "#9ca3af",
}

_KAROSPACE_DARK_PALETTE = {
    "plum_dark": _KI_COLORS["plum_dark"],
    "orange": _KI_COLORS["orange"],
    "light_orange": _KI_COLORS["light_orange"],
    "light_blue": _KI_COLORS["light_blue"],
    "plum": _KI_COLORS["plum"],
    "background": "#101114",
    "text": "#d8dae0",
    "panel_bg": "#191b20",
    "border": "#2a2c35",
    "input_bg": "#222430",
    "muted": "#6e7283",
    "hover_bg": "#252730",
    "accent": _KI_COLORS["orange"],
    "accent_strong": _KI_COLORS["plum"],
    "on_accent": "#1a1a1a",
    "hero_bg": "#191b20",
    "section_text": "#575b6e",
}


def _palette_for_mode(mode: str) -> dict[str, str]:
    return (
        dict(_KAROSPACE_DARK_PALETTE)
        if str(mode).strip().lower() == "dark"
        else dict(_KAROSPACE_LIGHT_PALETTE)
    )


def _ui_font() -> str:
    """Return the preferred UI font family for the current platform."""
    if sys.platform == "darwin":
        return "Helvetica Neue"
    return "Segoe UI"


def _mono_font() -> str:
    """Return the preferred monospace font family for the current platform."""
    if sys.platform == "darwin":
        return "Menlo"
    return "Consolas"


def _ctk_theme_config(palette: dict[str, str]) -> dict[str, dict[str, object]]:
    ui = _ui_font()
    mono = _mono_font()
    section_color = palette.get("section_text", palette["muted"])
    return {
        "root": {"fg_color": palette["background"]},
        "root_frame": {"fg_color": palette["background"], "corner_radius": 0},
        "card_frame": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 12,
            "border_width": 1,
            "border_color": palette["border"],
        },
        "hero_card": {
            "fg_color": palette.get("hero_bg", palette["panel_bg"]),
            "corner_radius": 12,
            "border_width": 0,
            "border_color": palette["border"],
        },
        "highlight_card": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 10,
            "border_width": 1,
            "border_color": palette["border"],
        },
        "sub_frame": {"fg_color": "transparent", "corner_radius": 0},
        "section_label": {
            "font": (ui, 10),
            "text_color": section_color,
            "fg_color": "transparent",
            "anchor": "w",
        },
        "hero_label": {
            "font": (ui, 20, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "header_label": {
            "font": (ui, 14, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "subheader_label": {
            "font": (ui, 11),
            "text_color": palette["muted"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "field_label": {
            "font": (ui, 11),
            "text_color": palette["muted"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "body_label": {
            "font": (ui, 10),
            "text_color": palette["muted"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "primary_button": {
            "fg_color": palette["accent"],
            "hover_color": palette["accent_strong"],
            "text_color": palette.get("on_accent", "#ffffff"),
            "corner_radius": 6,
            "font": (ui, 11, "bold"),
            "height": 32,
            "border_width": 0,
        },
        "secondary_button": {
            "fg_color": palette["input_bg"],
            "hover_color": palette["hover_bg"],
            "text_color": palette["text"],
            "corner_radius": 6,
            "font": (ui, 10),
            "height": 30,
            "border_width": 1,
            "border_color": palette["border"],
        },
        "pill_label": {
            "font": (ui, 9, "bold"),
            "text_color": palette.get("on_accent", "#ffffff"),
            "fg_color": palette["accent"],
            "corner_radius": 999,
            "padx": 10,
            "pady": 3,
        },
        "muted_pill_label": {
            "font": (ui, 9),
            "text_color": palette["muted"],
            "fg_color": palette["hover_bg"],
            "corner_radius": 999,
            "padx": 10,
            "pady": 3,
        },
        "entry": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "placeholder_text_color": palette["muted"],
            "border_color": palette["border"],
            "border_width": 1,
            "corner_radius": 6,
            "height": 30,
            "font": (ui, 11),
        },
        "combo": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "button_color": palette["hover_bg"],
            "button_hover_color": palette["accent"],
            "dropdown_fg_color": palette["panel_bg"],
            "dropdown_text_color": palette["text"],
            "dropdown_hover_color": palette["hover_bg"],
            "corner_radius": 6,
            "font": (ui, 11),
        },
        "checkbox": {
            "fg_color": palette["accent"],
            "hover_color": palette["accent_strong"],
            "checkmark_color": palette.get("on_accent", "#ffffff"),
            "text_color": palette["text"],
            "border_color": palette["border"],
            "font": (ui, 11),
            "corner_radius": 4,
        },
        "tabview": {
            "fg_color": palette["panel_bg"],
            "segmented_button_fg_color": palette["background"],
            "segmented_button_selected_color": palette["hover_bg"],
            "segmented_button_selected_hover_color": palette["hover_bg"],
            "segmented_button_unselected_color": palette["background"],
            "segmented_button_unselected_hover_color": palette["hover_bg"],
            "text_color": palette["text"],
            "corner_radius": 10,
            "border_width": 0,
            "border_color": palette["border"],
        },
        "divider": {
            "fg_color": palette["border"],
            "corner_radius": 0,
        },
        "textbox": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "border_color": palette["border"],
            "corner_radius": 8,
            "border_width": 1,
            "font": (mono, 11),
            "scrollbar_button_color": palette["hover_bg"],
            "scrollbar_button_hover_color": palette["accent"],
        },
        "progress": {
            "fg_color": palette["background"],
            "progress_color": palette["accent"],
            "border_color": palette["background"],
            "corner_radius": 999,
            "height": 4,
        },
    }


def _get_anndata():
    import anndata as ad

    return ad


def _get_numpy():
    import numpy as np

    return np


class SearchableListEditor(_CTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        *,
        label: str,
        height: int = 8,
        help_text: str | None = None,
        palette: dict[str, str] | None = None,
    ) -> None:
        self._palette = dict(palette or _palette_for_mode("dark"))
        self._theme = _ctk_theme_config(self._palette)
        super().__init__(parent, **self._theme["sub_frame"])
        self._choices: list[str] = []
        self._input_var = tk.StringVar(value="")

        self.columnconfigure(0, weight=1)
        self.label_widget = ctk.CTkLabel(self, text=label, **self._theme["field_label"])
        self.label_widget.grid(row=0, column=0, sticky="w", pady=(0, 3))

        controls = ctk.CTkFrame(self, **self._theme["sub_frame"])
        controls.grid(row=1, column=0, sticky="ew")
        controls.columnconfigure(0, weight=1)

        self.entry = ctk.CTkComboBox(
            controls,
            variable=self._input_var,
            values=[],
            state="normal",
            **self._theme["combo"],
        )
        self.entry.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self.entry.bind("<KeyRelease>", self._on_search)
        self.entry.bind("<Return>", lambda _event: self.add_current())

        self.add_btn = ctk.CTkButton(controls, text="+ Add", command=self.add_current, width=88, **self._theme["secondary_button"])
        self.add_btn.grid(row=0, column=1, padx=(0, 6))
        self.remove_btn = ctk.CTkButton(
            controls,
            text="Remove",
            command=self.remove_selected,
            width=88,
            **self._theme["secondary_button"],
        )
        self.remove_btn.grid(row=0, column=2, padx=(0, 6))
        self.clear_btn = ctk.CTkButton(controls, text="Clear", command=self.clear, width=88, **self._theme["secondary_button"])
        self.clear_btn.grid(row=0, column=3)

        list_wrap = ctk.CTkFrame(self, **self._theme["sub_frame"])
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=height,
            selectmode="extended",
            activestyle="none",
            relief="flat",
            bd=0,
            highlightthickness=1,
            font=(_ui_font(), 10),
        )
        self.listbox.grid(row=0, column=0, sticky="ew")
        self.scroll = tk.Scrollbar(list_wrap, orient="vertical", command=self.listbox.yview)
        self.scroll.grid(row=0, column=1, sticky="ns")
        self.listbox.configure(yscrollcommand=self.scroll.set)
        self.listbox.bind("<Command-v>", self._on_paste)
        self.listbox.bind("<Control-v>", self._on_paste)

        self.help_label: ctk.CTkLabel | None = None
        if help_text:
            self.help_label = ctk.CTkLabel(self, text=help_text, **self._theme["subheader_label"])
            self.help_label.grid(row=3, column=0, sticky="w", pady=(6, 0))

        self.apply_palette(self._palette)

    def _on_paste(self, _event) -> None:
        try:
            text = self.clipboard_get()
        except Exception:
            return "break"
        existing = set(self.get_items())
        import re
        items = [v.strip() for v in re.split(r"[,\n\r\t]+", text) if v.strip()]
        added = 0
        for item in items:
            if item not in existing:
                self.listbox.insert("end", item)
                existing.add(item)
                added += 1
        return "break"

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
            self.entry.set("")
            return
        self.listbox.insert("end", value)
        self._input_var.set("")
        self.entry.set("")
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
        self._palette = dict(palette)
        self._theme = _ctk_theme_config(self._palette)
        self.label_widget.configure(**self._theme["field_label"])
        self.entry.configure(**self._theme["combo"])
        self.add_btn.configure(**self._theme["secondary_button"])
        self.remove_btn.configure(**self._theme["secondary_button"])
        self.clear_btn.configure(**self._theme["secondary_button"])
        if self.help_label is not None:
            self.help_label.configure(**self._theme["subheader_label"])
        self.listbox.configure(
            background=self._palette["input_bg"],
            foreground=self._palette["text"],
            selectbackground=self._palette["accent"],
            selectforeground=self._palette.get("on_accent", "#ffffff"),
            disabledforeground=self._palette["muted"],
            highlightbackground=self._palette["border"],
            highlightcolor=self._palette["accent"],
        )
        try:
            self.scroll.configure(
                background=self._palette["panel_bg"],
                troughcolor=self._palette["hover_bg"],
                activebackground=self._palette["accent"],
                highlightbackground=self._palette["border"],
            )
        except tk.TclError:
            try:
                self.scroll.configure(background=self._palette["panel_bg"], activebackground=self._palette["accent"])
            except tk.TclError:
                pass


class SearchableMultiSelectEditor(_CTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        *,
        label: str,
        height: int = 8,
        help_text: str | None = None,
        max_visible: int = 500,
        palette: dict[str, str] | None = None,
    ) -> None:
        self._palette = dict(palette or _palette_for_mode("dark"))
        self._theme = _ctk_theme_config(self._palette)
        super().__init__(parent, **self._theme["sub_frame"])
        self._choices: list[str] = []
        self._selected_values: set[str] = set()
        self._search_var = tk.StringVar(value="")
        self._info_var = tk.StringVar(value="")
        self._max_visible = max(50, int(max_visible))

        self.columnconfigure(0, weight=1)
        self.label_widget = ctk.CTkLabel(self, text=label, **self._theme["field_label"])
        self.label_widget.grid(row=0, column=0, sticky="w", pady=(0, 3))

        controls = ctk.CTkFrame(self, **self._theme["sub_frame"])
        controls.grid(row=1, column=0, sticky="ew")
        controls.columnconfigure(0, weight=1)

        self.search_entry = ctk.CTkEntry(controls, textvariable=self._search_var, **self._theme["entry"])
        self.search_entry.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self.search_entry.bind("<KeyRelease>", lambda _event: self._on_search_changed())

        self.select_all_btn = ctk.CTkButton(
            controls,
            text="Select all",
            command=self.select_all_matches,
            width=100,
            **self._theme["secondary_button"],
        )
        self.select_all_btn.grid(row=0, column=1, padx=(0, 6))
        self.clear_btn = ctk.CTkButton(
            controls,
            text="Clear",
            command=self.clear_selection,
            width=90,
            **self._theme["secondary_button"],
        )
        self.clear_btn.grid(row=0, column=2)

        list_wrap = ctk.CTkFrame(self, **self._theme["sub_frame"])
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(4, 0))
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=height,
            selectmode="extended",
            activestyle="none",
            relief="flat",
            bd=0,
            highlightthickness=1,
            font=(_ui_font(), 10),
        )
        self.listbox.grid(row=0, column=0, sticky="ew")
        self.listbox.bind("<<ListboxSelect>>", lambda _event: self._capture_visible_selection())
        self.scroll = tk.Scrollbar(list_wrap, orient="vertical", command=self.listbox.yview)
        self.scroll.grid(row=0, column=1, sticky="ns")
        self.listbox.configure(yscrollcommand=self.scroll.set)

        self.info_label = ctk.CTkLabel(self, textvariable=self._info_var, **self._theme["subheader_label"])
        self.info_label.grid(row=3, column=0, sticky="w", pady=(6, 0))

        self.help_label: ctk.CTkLabel | None = None
        if help_text:
            self.help_label = ctk.CTkLabel(self, text=help_text, **self._theme["subheader_label"])
            self.help_label.grid(row=4, column=0, sticky="w", pady=(4, 0))

        self.apply_palette(self._palette)

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
        self._palette = dict(palette)
        self._theme = _ctk_theme_config(self._palette)
        self.label_widget.configure(**self._theme["field_label"])
        self.search_entry.configure(**self._theme["entry"])
        self.select_all_btn.configure(**self._theme["secondary_button"])
        self.clear_btn.configure(**self._theme["secondary_button"])
        self.info_label.configure(**self._theme["subheader_label"])
        if self.help_label is not None:
            self.help_label.configure(**self._theme["subheader_label"])
        self.listbox.configure(
            background=self._palette["input_bg"],
            foreground=self._palette["text"],
            selectbackground=self._palette["accent"],
            selectforeground=self._palette.get("on_accent", "#ffffff"),
            disabledforeground=self._palette["muted"],
            highlightbackground=self._palette["border"],
            highlightcolor=self._palette["accent"],
        )
        try:
            self.scroll.configure(
                background=self._palette["panel_bg"],
                troughcolor=self._palette["hover_bg"],
                activebackground=self._palette["accent"],
                highlightbackground=self._palette["border"],
            )
        except tk.TclError:
            try:
                self.scroll.configure(background=self._palette["panel_bg"], activebackground=self._palette["accent"])
            except tk.TclError:
                pass


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
    enable_numba_jit: bool
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
    interaction_markers_method: str
    interaction_markers_layer: str | None
    # Output format
    output_format: str  # "html" | "karospace"
    # Gene storage / encoding
    gene_encoding: str
    gene_value_encoding: str
    gene_sidecar_format: str
    gene_storage: str
    gene_aux_path: str | None
    gene_sidecar_shard_size: int
    gene_sparse_zero_threshold: float
    pack_arrays: bool
    pack_arrays_min_len: int
    # Cluster DE
    cluster_de_enabled: bool
    cluster_de_groupby: list[str] | None
    cluster_de_top_n: int
    cluster_de_method: str
    cluster_de_layer: str | None
    cluster_de_min_cells: int
    # Gene analysis
    gene_correlation_top_n: int
    cluster_means_n_genes: int
    spatial_variable_genes_n: int
    # Viewer
    section_rotations: dict[str, float] | None
    scalebar_unit: str
    vmin: float | None
    vmax: float | None
    viewer_info_html: str | None
    # load_spatial_data metadata
    group_order: list[str] | None
    metadata_columns: list[str] | None
    metadata_value_order: dict[str, list[str]] | None
    metadata_max_columns: int | None


class _ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


class ExportApp(ctk.CTk if ctk is not None else object):
    _OUTPUT_HTML_BASENAME = "KaroSpace"
    _NEIGHBOR_MATRIX_BUDGET_BYTES = 512 * 1024 * 1024
    _MARKER_GROUPBY_MAX_UNIQUE = 2048
    _INTERACTION_GROUPBY_MAX_UNIQUE = 512

    def __init__(self) -> None:
        super().__init__()
        self.title("KaroSpaceBuilder")
        self.geometry("1280x820")
        self.minsize(1060, 720)

        self._queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self._export_thread: threading.Thread | None = None
        self._server: _ThreadingHTTPServer | None = None
        self._server_thread: threading.Thread | None = None
        self._last_outdir: Path | None = None
        self._last_output_html: Path | None = None
        self._inspected_h5ad_path: Path | None = None
        self._inspected_var_name_set: set[str] | None = None
        self._inspected_coords_mode: str | None = None
        self._themed_widgets: dict[str, list[tk.Widget]] = {}

        self._build_style()
        self._build_variables()
        self._build_layout()
        self.after(120, self._poll_events)
        self.protocol("WM_DELETE_WINDOW", self._on_close)

    def _build_style(self) -> None:
        self._themed_widgets = {}
        self._app_palette = _palette_for_mode("dark")
        self._theme = _ctk_theme_config(self._app_palette)
        ctk.set_appearance_mode("dark")
        self.configure(**self._theme["root"])

    def _apply_app_theme(self, mode: str) -> None:
        palette = _palette_for_mode(mode)
        self._app_palette = palette
        self._theme = _ctk_theme_config(palette)
        ctk.set_appearance_mode("dark" if mode.strip().lower() == "dark" else "light")
        self.configure(**self._theme["root"])

        for role, widgets in self._themed_widgets.items():
            role_style = self._theme.get(role)
            if not role_style:
                continue
            for widget in widgets:
                if widget is None:
                    continue
                try:
                    widget.configure(**role_style)
                except Exception:
                    continue

        if hasattr(self, "help_text"):
            self.help_text.configure(**self._theme["textbox"])
        if hasattr(self, "log_text"):
            self.log_text.configure(**self._theme["textbox"])
        if hasattr(self, "progress"):
            self.progress.configure(**self._theme["progress"])
        if hasattr(self, "main_scroll_frame"):
            self.main_scroll_frame.configure(
                scrollbar_button_color=palette["hover_bg"],
                scrollbar_button_hover_color=palette["accent"],
            )
        self._sync_runtime_chip()

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

    def _register_theme_widget(self, role: str, widget: tk.Widget) -> tk.Widget:
        self._themed_widgets.setdefault(role, []).append(widget)
        return widget

    def _header_label(self, parent: tk.Widget, text: str) -> ctk.CTkLabel:
        label = ctk.CTkLabel(parent, text=text, **self._theme["header_label"])
        self._register_theme_widget("header_label", label)
        return label

    def _hero_label(self, parent: tk.Widget, text: str) -> ctk.CTkLabel:
        label = ctk.CTkLabel(parent, text=text, **self._theme["hero_label"])
        self._register_theme_widget("hero_label", label)
        return label

    def _subheader_label(self, parent: tk.Widget, text: str | None = None, textvariable: tk.StringVar | None = None) -> ctk.CTkLabel:
        kwargs: dict[str, object] = dict(self._theme["subheader_label"])
        if text is not None:
            kwargs["text"] = text
        if textvariable is not None:
            kwargs["textvariable"] = textvariable
        label = ctk.CTkLabel(parent, **kwargs)
        self._register_theme_widget("subheader_label", label)
        return label

    def _section_label(self, parent: tk.Widget, text: str) -> ctk.CTkLabel:
        label = ctk.CTkLabel(parent, text=text, **self._theme["section_label"])
        self._register_theme_widget("section_label", label)
        return label

    def _field_label(self, parent: tk.Widget, text: str) -> ctk.CTkLabel:
        label = ctk.CTkLabel(parent, text=text, **self._theme["field_label"])
        self._register_theme_widget("field_label", label)
        return label

    def _body_label(self, parent: tk.Widget, text: str) -> ctk.CTkLabel:
        label = ctk.CTkLabel(parent, text=text, **self._theme["body_label"])
        self._register_theme_widget("body_label", label)
        return label

    def _secondary_button(self, parent: tk.Widget, text: str, command, width: int | None = None) -> ctk.CTkButton:
        kwargs: dict[str, object] = {"text": text, "command": command, **self._theme["secondary_button"]}
        if width is not None:
            kwargs["width"] = width
        button = ctk.CTkButton(parent, **kwargs)
        self._register_theme_widget("secondary_button", button)
        return button

    def _primary_button(self, parent: tk.Widget, text: str, command, width: int | None = None) -> ctk.CTkButton:
        kwargs: dict[str, object] = {"text": text, "command": command, **self._theme["primary_button"]}
        if width is not None:
            kwargs["width"] = width
        button = ctk.CTkButton(parent, **kwargs)
        self._register_theme_widget("primary_button", button)
        return button

    def _pill_label(
        self,
        parent: tk.Widget,
        *,
        text: str | None = None,
        textvariable: tk.StringVar | None = None,
        muted: bool = False,
    ) -> ctk.CTkLabel:
        role = "muted_pill_label" if muted else "pill_label"
        kwargs: dict[str, object] = dict(self._theme[role])
        if text is not None:
            kwargs["text"] = text
        if textvariable is not None:
            kwargs["textvariable"] = textvariable
        label = ctk.CTkLabel(parent, **kwargs)
        self._register_theme_widget(role, label)
        return label

    def _divider(self, parent: tk.Widget, *, height: int = 1) -> ctk.CTkFrame:
        frame = ctk.CTkFrame(parent, height=height, **self._theme["divider"])
        self._register_theme_widget("divider", frame)
        return frame

    def _make_card_frame(self, parent: tk.Widget, *, padding: int = 0) -> ctk.CTkFrame:
        frame = ctk.CTkFrame(parent, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", frame)
        if padding > 0:
            inner = ctk.CTkFrame(frame, **self._theme["sub_frame"])
            self._register_theme_widget("sub_frame", inner)
            inner.pack(fill="both", expand=True, padx=padding, pady=padding)
            return inner
        return frame

    def _make_sub_frame(self, parent: tk.Widget) -> ctk.CTkFrame:
        frame = ctk.CTkFrame(parent, **self._theme["sub_frame"])
        self._register_theme_widget("sub_frame", frame)
        return frame

    def _register_entry_widget(self, widget: tk.Widget) -> None:
        self._register_theme_widget("entry", widget)

    def _register_combo_widget(self, widget: tk.Widget) -> None:
        self._register_theme_widget("combo", widget)

    def _register_checkbox_widget(self, widget: tk.Widget) -> None:
        self._register_theme_widget("checkbox", widget)

    def _sync_app_theme_to_viewer_setting(self) -> None:
        self._apply_app_theme(self.theme_var.get())

    def _on_theme_selected(self, selected_mode: str) -> None:
        mode = str(selected_mode).strip().lower()
        if mode not in {"light", "dark"}:
            mode = "dark"
        if self.theme_var.get().strip().lower() != mode:
            self.theme_var.set(mode)
        self._apply_app_theme(mode)

    def _sync_runtime_chip(self) -> None:
        if not hasattr(self, "runtime_chip_label"):
            return
        status = self.status_var.get().strip().lower()
        if status.startswith("export running"):
            text = "RUNNING"
            fg = self._app_palette["accent"]
            tc = self._app_palette.get("on_accent", "#ffffff")
        elif status.startswith("export complete"):
            text = "COMPLETE"
            fg = "#1f8f5f"
            tc = "#ffffff"
        elif status.startswith("serving on"):
            text = "SERVING"
            fg = self._app_palette["accent_strong"]
            tc = self._app_palette.get("on_accent", "#ffffff")
        elif status.startswith("export failed"):
            text = "FAILED"
            fg = "#bf2f5e"
            tc = "#ffffff"
        else:
            text = "READY"
            fg = self._app_palette["hover_bg"]
            tc = self._app_palette["text"]
        self.runtime_chip_label.configure(text=text, fg_color=fg, text_color=tc)

    def _build_variables(self) -> None:
        self.h5ad_var = tk.StringVar()
        self.outdir_var = tk.StringVar()
        self.coords_var = tk.StringVar(value="auto")
        self.section_groupby_var = tk.StringVar(value="sample_id")
        self.initial_color_var = tk.StringVar(value="leiden")
        self.title_var = tk.StringVar(value="KaroSpace")
        self.theme_var = tk.StringVar(value="dark")
        self.numba_jit_var = tk.BooleanVar(value=False)
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
        self.interaction_markers_method_var = tk.StringVar(value="wilcoxon")
        self.interaction_markers_layer_var = tk.StringVar(value="normalized")
        self.selection_mode_var = tk.BooleanVar(value=False)

        self.downsample_var = tk.StringVar()

        # Output format
        self.output_format_var = tk.StringVar(value="embedded")

        # Gene storage / encoding
        self.gene_encoding_var = tk.StringVar(value="auto")
        self.gene_value_encoding_var = tk.StringVar(value="float32")
        self.gene_sidecar_format_var = tk.StringVar(value="json-v2")
        self.gene_storage_var = tk.StringVar(value="embedded")
        self.gene_aux_path_var = tk.StringVar()
        self.gene_sidecar_shard_size_var = tk.StringVar(value="256")
        self.gene_sparse_zero_threshold_var = tk.StringVar(value="0.8")
        self.pack_arrays_var = tk.BooleanVar(value=True)
        self.pack_arrays_min_len_var = tk.StringVar(value="1024")

        # Cluster DE
        self.cluster_de_enabled_var = tk.BooleanVar(value=False)
        self.cluster_de_top_n_var = tk.StringVar(value="20")
        self.cluster_de_method_var = tk.StringVar(value="wilcoxon")
        self.cluster_de_layer_var = tk.StringVar(value="normalized")
        self.cluster_de_min_cells_var = tk.StringVar(value="20")

        # Gene analysis
        self.gene_correlation_top_n_var = tk.StringVar(value="10")
        self.cluster_means_n_genes_var = tk.StringVar(value="500")
        self.spatial_variable_genes_n_var = tk.StringVar(value="200")

        # Viewer
        self.section_rotations_var = tk.StringVar()
        self.scalebar_unit_var = tk.StringVar(value="\u03bcm")
        self.vmin_var = tk.StringVar()
        self.vmax_var = tk.StringVar()
        self.viewer_info_html_var = tk.StringVar()

        # load_spatial_data metadata
        self.group_order_var = tk.StringVar()
        self.metadata_columns_var = tk.StringVar()
        self.metadata_value_order_var = tk.StringVar()
        self.metadata_max_columns_var = tk.StringVar()

        self.serve_var = tk.BooleanVar(value=False)
        self.port_var = tk.StringVar(value="8000")

        self.status_var = tk.StringVar(value="Ready")

    def _build_layout(self) -> None:
        shell = ctk.CTkFrame(self, **self._theme["root_frame"])
        self._register_theme_widget("root_frame", shell)
        shell.pack(fill="both", expand=True)
        shell.columnconfigure(0, weight=1)
        shell.rowconfigure(0, weight=1)

        root = ctk.CTkScrollableFrame(
            shell,
            **self._theme["root_frame"],
            scrollbar_button_color=self._app_palette["hover_bg"],
            scrollbar_button_hover_color=self._app_palette["accent"],
        )
        self.main_scroll_frame = root
        self._register_theme_widget("root_frame", root)
        root.grid(row=0, column=0, sticky="nsew", padx=10, pady=10)
        root.columnconfigure(0, weight=5)
        root.columnconfigure(1, weight=2)
        root.rowconfigure(0, weight=1)

        controls = ctk.CTkFrame(root, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", controls)
        controls.grid(row=0, column=0, sticky="nsew", padx=(0, 8))

        side = ctk.CTkFrame(root, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", side)
        side.grid(row=0, column=1, sticky="nsew")

        controls_inner = ctk.CTkFrame(controls, **self._theme["sub_frame"])
        self._register_theme_widget("sub_frame", controls_inner)
        controls_inner.pack(fill="both", expand=True, padx=14, pady=14)
        controls = controls_inner

        side_inner = ctk.CTkFrame(side, **self._theme["sub_frame"])
        self._register_theme_widget("sub_frame", side_inner)
        side_inner.pack(fill="both", expand=True, padx=14, pady=14)
        side = side_inner

        controls.columnconfigure(1, weight=1)
        side.columnconfigure(0, weight=1)
        side.rowconfigure(6, weight=1)

        hero = ctk.CTkFrame(controls, **self._theme["hero_card"])
        self._register_theme_widget("hero_card", hero)
        hero.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        hero.columnconfigure(0, weight=1)

        hero_inner = self._make_sub_frame(hero)
        hero_inner.grid(row=0, column=0, sticky="ew", padx=14, pady=10)
        hero_inner.columnconfigure(0, weight=1)

        self._hero_label(hero_inner, "KaroSpaceBuilder").pack(anchor="w")
        self._subheader_label(
            hero_inner,
            "Build interactive KaroSpace viewers from AnnData.",
        ).pack(anchor="w", pady=(1, 0))

        self._divider(controls, height=1).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 8))

        notebook = ctk.CTkTabview(controls, **self._theme["tabview"])
        self._register_theme_widget("tabview", notebook)
        notebook.grid(row=2, column=0, columnspan=3, sticky="nsew")
        controls.rowconfigure(2, weight=1)

        notebook.add("Basic")
        notebook.add("Colors & Genes")
        notebook.add("Export Format")
        notebook.add("Advanced")
        notebook.add("Help")
        basic_tab = notebook.tab("Basic")
        colors_tab = notebook.tab("Colors & Genes")
        advanced_tab = notebook.tab("Advanced")
        export_fmt_tab = notebook.tab("Export Format")
        help_tab = notebook.tab("Help")
        for tab in (basic_tab, colors_tab, advanced_tab, export_fmt_tab, help_tab):
            self._register_theme_widget("sub_frame", tab)

        basic_tab.columnconfigure(1, weight=1)
        self._section_label(basic_tab, "INPUT & VIEWER SETUP").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        self._divider(basic_tab, height=1).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        row = 2
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
        runtime_mode_row = self._make_sub_frame(basic_tab)
        self.numba_jit_check = ctk.CTkCheckBox(
            runtime_mode_row,
            text="Performance mode (enable numba JIT)",
            variable=self.numba_jit_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.numba_jit_check)
        self.numba_jit_check.pack(side="left")
        row = self._option_row(
            basic_tab,
            row,
            "Runtime mode",
            widget=runtime_mode_row,
            hint="Faster on some datasets. If unstable in the desktop app, turn this off.",
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
        downsample_container = self._make_sub_frame(basic_tab)
        downsample_entry = ctk.CTkEntry(downsample_container, textvariable=self.downsample_var, width=100, **self._theme["entry"])
        self._register_entry_widget(downsample_entry)
        downsample_entry.pack(side="left")
        self._body_label(downsample_container, "cells per section (blank = all)").pack(side="left", padx=(8, 0))
        row = self._option_row(
            basic_tab,
            row,
            "Downsample",
            widget=downsample_container,
            hint="Maps to export_to_html(downsample=...).",
        )

        self._divider(basic_tab, height=1).grid(row=row, column=0, columnspan=3, sticky="ew", pady=(6, 10))
        row += 1
        self._section_label(basic_tab, "METADATA LOADING").grid(row=row, column=0, columnspan=3, sticky="w", pady=(0, 6))
        row += 1
        self._subheader_label(
            basic_tab,
            "Optional parameters for load_spatial_data() — affects section grouping and metadata columns.",
        ).grid(row=row, column=1, columnspan=2, sticky="w", pady=(0, 8))
        row += 1
        row = self._option_row(
            basic_tab,
            row,
            "Group order",
            widget=self._entry(basic_tab, self.group_order_var),
            hint="Comma-separated section order (e.g. sample_A,sample_B). Blank = natural order.",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Metadata columns",
            widget=self._entry(basic_tab, self.metadata_columns_var),
            hint="Comma-separated extra obs columns to include. Blank = auto.",
        )
        row = self._option_row(
            basic_tab,
            row,
            "Metadata max cols",
            widget=self._entry(basic_tab, self.metadata_max_columns_var),
            hint="Max metadata columns to load. Blank = unlimited.",
        )

        colors_tab.columnconfigure(0, weight=1)
        self._section_label(colors_tab, "COLOR FIELDS & GENE SOURCES").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._subheader_label(
            colors_tab,
            "Build color menus and gene features from inspected obs/var fields.",
        ).grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._divider(colors_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        self.additional_colors_editor = SearchableListEditor(
            colors_tab,
            label="additional_colors (obs columns)",
            height=4,
            help_text="These become color options in KaroSpace. Use Inspect to load obs columns.",
            palette=self._app_palette,
        )
        self.additional_colors_editor.grid(row=3, column=0, sticky="ew", pady=(0, 10))

        self.groupby_editor = SearchableListEditor(
            colors_tab,
            label="groupby lists (marker/neighbor/interaction/cluster DE)",
            height=4,
            help_text="Used for marker_genes_groupby, cluster_de_groupby, and optionally neighbor/interaction groupby.",
            palette=self._app_palette,
        )
        self.groupby_editor.grid(row=4, column=0, sticky="ew", pady=(0, 4))

        autofill_row = self._make_sub_frame(colors_tab)
        autofill_row.grid(row=5, column=0, sticky="w", pady=(0, 12))
        self._secondary_button(
            autofill_row, "Auto-fill groupby from colors", self._autofill_groupby_from_colors, width=230,
        ).pack(side="left")

        genes_card = ctk.CTkFrame(colors_tab, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", genes_card)
        genes_card.grid(row=6, column=0, sticky="ew", pady=(0, 0))
        genes_card.columnconfigure(0, weight=1)
        genes_inner = self._make_sub_frame(genes_card)
        genes_inner.grid(row=0, column=0, sticky="ew", padx=12, pady=12)
        genes_inner.columnconfigure(0, weight=1)
        self._field_label(genes_inner, "Gene Selection").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._genes_mode_row(genes_inner).grid(row=1, column=0, sticky="ew")
        self._subheader_label(
            genes_inner,
            (
                "hvgs:N maps to use_hvgs=True/hvg_limit=N. top_mean/list_file/manual_list map to explicit genes list "
                "with use_hvgs=False."
            ),
        ).grid(row=2, column=0, sticky="w", pady=(4, 8))
        self.manual_genes_editor = SearchableListEditor(
            genes_inner,
            label="genes list (manual_list)",
            height=5,
            help_text="Search var_names and build the genes list with + Add / Remove.",
            palette=self._app_palette,
        )
        self.manual_genes_editor.grid(row=3, column=0, sticky="ew")

        selection_card = ctk.CTkFrame(colors_tab, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", selection_card)
        selection_card.grid(row=7, column=0, sticky="ew", pady=(12, 0))
        selection_card.columnconfigure(0, weight=1)

        selection_inner = self._make_sub_frame(selection_card)
        selection_inner.grid(row=0, column=0, sticky="ew", padx=12, pady=12)
        selection_inner.columnconfigure(0, weight=1)
        self.selection_mode_check = ctk.CTkCheckBox(
            selection_inner,
            text="Enable tick selection mode (from Inspect H5AD)",
            variable=self.selection_mode_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.selection_mode_check)
        self.selection_mode_check.grid(row=0, column=0, sticky="w")
        self._subheader_label(
            selection_inner,
            (
                "When enabled, selected items below are used during export for additional colors, "
                "groupby lists, and manual genes."
            ),
        ).grid(row=1, column=0, sticky="w", pady=(4, 8))

        self.selection_mode_content = self._make_sub_frame(selection_inner)
        self.selection_mode_content.grid(row=2, column=0, sticky="ew")
        self.selection_mode_content.columnconfigure(0, weight=1)
        self.selection_mode_content.columnconfigure(1, weight=1)

        self.selection_additional_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick additional_colors (obs)",
            height=5,
            help_text="Multi-select obs fields to include as additional viewer colors.",
            palette=self._app_palette,
        )
        self.selection_additional_picker.grid(row=0, column=0, sticky="ew", padx=(0, 8))

        self.selection_groupby_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick groupby lists (obs)",
            height=5,
            help_text="Multi-select obs columns for marker/neighbor/interaction groupby lists.",
            palette=self._app_palette,
        )
        self.selection_groupby_picker.grid(row=0, column=1, sticky="ew")

        self.selection_genes_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick genes (manual_list mode)",
            height=5,
            help_text="Multi-select genes from var_names. Used when genes mode is manual_list.",
            palette=self._app_palette,
        )
        self.selection_genes_picker.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(10, 0))

        selection_actions = self._make_sub_frame(self.selection_mode_content)
        selection_actions.grid(row=2, column=0, columnspan=2, sticky="w", pady=(10, 0))
        self.selection_apply_btn = self._secondary_button(
            selection_actions,
            "Apply selections to inputs",
            self._apply_tick_selections_to_inputs,
            width=190,
        )
        self.selection_apply_btn.pack(side="left")
        self.selection_sync_btn = self._secondary_button(
            selection_actions,
            "Sync from current inputs",
            self._load_inputs_into_tick_selection,
            width=185,
        )
        self.selection_sync_btn.pack(side="left", padx=(8, 0))

        advanced_tab.columnconfigure(0, weight=1)
        self._section_label(advanced_tab, "ADVANCED ANALYTICS").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._subheader_label(
            advanced_tab,
            "Fine tune marker, neighbor, interaction, and preview server behavior.",
        ).grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._divider(advanced_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        adv_header = self._make_sub_frame(advanced_tab)
        adv_header.grid(row=3, column=0, sticky="ew")
        self.advanced_toggle_btn = self._secondary_button(adv_header, "Show Advanced Options", self._toggle_advanced, width=190)
        self.advanced_toggle_btn.pack(anchor="w")
        self._subheader_label(
            adv_header,
            "Analytics and rendering parameters from KaroSpace export_to_html.",
        ).pack(anchor="w", pady=(6, 0))

        self.advanced_content = self._make_sub_frame(advanced_tab)
        self.advanced_content.grid(row=4, column=0, sticky="ew", pady=(10, 0))
        self.advanced_content.columnconfigure(1, weight=1)

        min_panel_row = self._make_sub_frame(self.advanced_content)
        min_panel_entry = ctk.CTkEntry(min_panel_row, textvariable=self.min_panel_size_var, width=90, **self._theme["entry"])
        self._register_entry_widget(min_panel_entry)
        min_panel_entry.pack(side="left")
        self._body_label(min_panel_row, "px").pack(side="left", padx=(8, 0))
        self._option_row(
            self.advanced_content,
            0,
            "Min panel size",
            widget=min_panel_row,
            hint="Minimum section panel width in exported KaroSpace HTML.",
        )

        spot_row = self._make_sub_frame(self.advanced_content)
        self.spot_size_combo = ctk.CTkComboBox(
            spot_row,
            variable=self.spot_size_var,
            values=["auto", "adaptive", "density"],
            width=120,
            state="normal",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.spot_size_combo)
        self.spot_size_combo.pack(side="left")
        self._body_label(spot_row, "or numeric value").pack(side="left", padx=(8, 0))
        self._option_row(
            self.advanced_content,
            2,
            "Spot size",
            widget=spot_row,
            hint="Use auto/adaptive/density or a positive number.",
        )

        marker_row = self._make_sub_frame(self.advanced_content)
        marker_top_n_entry = ctk.CTkEntry(marker_row, textvariable=self.marker_genes_top_n_var, width=90, **self._theme["entry"])
        self._register_entry_widget(marker_top_n_entry)
        marker_top_n_entry.pack(side="left")
        self._body_label(marker_row, "top genes per group").pack(side="left", padx=(8, 0))
        self._option_row(
            self.advanced_content,
            4,
            "Marker genes top N",
            widget=marker_row,
            hint="Maps to marker_genes_top_n.",
        )

        neighbor_row = self._make_sub_frame(self.advanced_content)
        self.neighbor_auto_check = ctk.CTkCheckBox(
            neighbor_row,
            text="Auto neighbor groupby = [initial color]",
            variable=self.neighbor_auto_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.neighbor_auto_check)
        self.neighbor_auto_check.pack(side="left")
        self._body_label(neighbor_row, "Permutations").pack(side="left", padx=(12, 4))
        neighbor_perm_entry = ctk.CTkEntry(neighbor_row, textvariable=self.neighbor_permutations_var, width=76, **self._theme["entry"])
        self._register_entry_widget(neighbor_perm_entry)
        neighbor_perm_entry.pack(side="left")
        self._body_label(neighbor_row, "Seed").pack(side="left", padx=(12, 4))
        neighbor_seed_entry = ctk.CTkEntry(neighbor_row, textvariable=self.neighbor_stats_seed_var, width=76, **self._theme["entry"])
        self._register_entry_widget(neighbor_seed_entry)
        neighbor_seed_entry.pack(side="left")
        self._option_row(
            self.advanced_content,
            6,
            "Neighbor stats",
            widget=neighbor_row,
            hint="Permutations accepts integer or 'auto'.",
        )

        interaction_row_1 = self._make_sub_frame(self.advanced_content)
        self.interaction_enabled_check = ctk.CTkCheckBox(
            interaction_row_1,
            text="Enable interaction markers (uses groupby list)",
            variable=self.interaction_markers_enabled_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.interaction_enabled_check)
        self.interaction_enabled_check.pack(side="left")
        self._option_row(self.advanced_content, 8, "Interaction markers", widget=interaction_row_1)

        interaction_row_2 = self._make_sub_frame(self.advanced_content)
        self._body_label(interaction_row_2, "Top targets").pack(side="left")
        interaction_top_targets_entry = ctk.CTkEntry(
            interaction_row_2,
            textvariable=self.interaction_markers_top_targets_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(interaction_top_targets_entry)
        interaction_top_targets_entry.pack(side="left", padx=(4, 10))
        self._body_label(interaction_row_2, "Top genes").pack(side="left")
        interaction_top_genes_entry = ctk.CTkEntry(
            interaction_row_2,
            textvariable=self.interaction_markers_top_genes_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(interaction_top_genes_entry)
        interaction_top_genes_entry.pack(side="left", padx=(4, 10))
        self._body_label(interaction_row_2, "Min cells").pack(side="left")
        interaction_min_cells_entry = ctk.CTkEntry(
            interaction_row_2,
            textvariable=self.interaction_markers_min_cells_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(interaction_min_cells_entry)
        interaction_min_cells_entry.pack(side="left", padx=(4, 10))
        self._body_label(interaction_row_2, "Min neighbors").pack(side="left")
        interaction_min_neighbors_entry = ctk.CTkEntry(
            interaction_row_2,
            textvariable=self.interaction_markers_min_neighbors_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(interaction_min_neighbors_entry)
        interaction_min_neighbors_entry.pack(side="left", padx=(4, 0))
        self._option_row(
            self.advanced_content,
            9,
            "Interaction limits",
            widget=interaction_row_2,
            hint="Maps to interaction_markers_top_targets/top_genes/min_cells/min_neighbors.",
        )

        interaction_method_row = self._make_sub_frame(self.advanced_content)
        self._body_label(interaction_method_row, "Method").pack(side="left")
        self.interaction_method_combo = ctk.CTkComboBox(
            interaction_method_row,
            variable=self.interaction_markers_method_var,
            values=["wilcoxon", "t-test", "logreg"],
            width=110,
            state="normal",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.interaction_method_combo)
        self.interaction_method_combo.pack(side="left", padx=(4, 12))
        self._body_label(interaction_method_row, "Layer").pack(side="left")
        interaction_layer_entry = ctk.CTkEntry(
            interaction_method_row, textvariable=self.interaction_markers_layer_var, width=100, **self._theme["entry"],
        )
        self._register_entry_widget(interaction_layer_entry)
        interaction_layer_entry.pack(side="left", padx=(4, 0))
        self._option_row(
            self.advanced_content,
            10,
            "Interaction method",
            widget=interaction_method_row,
            hint="DE method and data layer for interaction markers.",
        )

        serve_row = self._make_sub_frame(self.advanced_content)
        self.serve_check = ctk.CTkCheckBox(serve_row, text="Serve after export", variable=self.serve_var, **self._theme["checkbox"])
        self._register_checkbox_widget(self.serve_check)
        self.serve_check.pack(side="left")
        self._body_label(serve_row, "Port").pack(side="left", padx=(12, 4))
        serve_port_entry = ctk.CTkEntry(serve_row, textvariable=self.port_var, width=90, **self._theme["entry"])
        self._register_entry_widget(serve_port_entry)
        serve_port_entry.pack(side="left")
        self._option_row(
            self.advanced_content,
            12,
            "Preview server",
            widget=serve_row,
            hint="Optional local server to open the latest generated KaroSpace_*.html file.",
        )

        # --- Cluster DE section ---
        self._divider(self.advanced_content, height=1).grid(row=14, column=0, columnspan=3, sticky="ew", pady=(10, 6))
        self._section_label(self.advanced_content, "CLUSTER DE").grid(row=15, column=0, columnspan=3, sticky="w", pady=(0, 6))

        cluster_de_row_1 = self._make_sub_frame(self.advanced_content)
        self.cluster_de_enabled_check = ctk.CTkCheckBox(
            cluster_de_row_1,
            text="Enable cluster DE (uses groupby list)",
            variable=self.cluster_de_enabled_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.cluster_de_enabled_check)
        self.cluster_de_enabled_check.pack(side="left")
        self._option_row(self.advanced_content, 16, "Cluster DE", widget=cluster_de_row_1)

        cluster_de_row_2 = self._make_sub_frame(self.advanced_content)
        self._body_label(cluster_de_row_2, "Top N").pack(side="left")
        cluster_de_top_n_entry = ctk.CTkEntry(
            cluster_de_row_2, textvariable=self.cluster_de_top_n_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(cluster_de_top_n_entry)
        cluster_de_top_n_entry.pack(side="left", padx=(4, 12))
        self._body_label(cluster_de_row_2, "Method").pack(side="left")
        self.cluster_de_method_combo = ctk.CTkComboBox(
            cluster_de_row_2,
            variable=self.cluster_de_method_var,
            values=["wilcoxon", "t-test", "logreg"],
            width=110,
            state="normal",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.cluster_de_method_combo)
        self.cluster_de_method_combo.pack(side="left", padx=(4, 12))
        self._body_label(cluster_de_row_2, "Layer").pack(side="left")
        cluster_de_layer_entry = ctk.CTkEntry(
            cluster_de_row_2, textvariable=self.cluster_de_layer_var, width=100, **self._theme["entry"],
        )
        self._register_entry_widget(cluster_de_layer_entry)
        cluster_de_layer_entry.pack(side="left", padx=(4, 12))
        self._body_label(cluster_de_row_2, "Min cells").pack(side="left")
        cluster_de_min_cells_entry = ctk.CTkEntry(
            cluster_de_row_2, textvariable=self.cluster_de_min_cells_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(cluster_de_min_cells_entry)
        cluster_de_min_cells_entry.pack(side="left", padx=(4, 0))
        self._option_row(
            self.advanced_content,
            17,
            "Cluster DE params",
            widget=cluster_de_row_2,
            hint="Pairwise cluster differential expression. Produces volcano plots in viewer.",
        )

        # --- Gene analysis section ---
        self._divider(self.advanced_content, height=1).grid(row=19, column=0, columnspan=3, sticky="ew", pady=(10, 6))
        self._section_label(self.advanced_content, "GENE ANALYSIS").grid(row=20, column=0, columnspan=3, sticky="w", pady=(0, 6))

        gene_analysis_row = self._make_sub_frame(self.advanced_content)
        self._body_label(gene_analysis_row, "Correlation top N").pack(side="left")
        gene_corr_entry = ctk.CTkEntry(
            gene_analysis_row, textvariable=self.gene_correlation_top_n_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(gene_corr_entry)
        gene_corr_entry.pack(side="left", padx=(4, 12))
        self._body_label(gene_analysis_row, "Cluster means genes").pack(side="left")
        cluster_means_entry = ctk.CTkEntry(
            gene_analysis_row, textvariable=self.cluster_means_n_genes_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(cluster_means_entry)
        cluster_means_entry.pack(side="left", padx=(4, 12))
        self._body_label(gene_analysis_row, "Spatial var genes").pack(side="left")
        spatial_var_entry = ctk.CTkEntry(
            gene_analysis_row, textvariable=self.spatial_variable_genes_n_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(spatial_var_entry)
        spatial_var_entry.pack(side="left", padx=(4, 0))
        self._option_row(
            self.advanced_content,
            21,
            "Gene analysis",
            widget=gene_analysis_row,
            hint="Gene correlation, cluster mean expression, and Moran's I spatial variable genes.",
        )

        # --- Viewer rendering section ---
        self._divider(self.advanced_content, height=1).grid(row=23, column=0, columnspan=3, sticky="ew", pady=(10, 6))
        self._section_label(self.advanced_content, "VIEWER RENDERING").grid(row=24, column=0, columnspan=3, sticky="w", pady=(0, 6))

        scalebar_row = self._make_sub_frame(self.advanced_content)
        self._body_label(scalebar_row, "Scalebar unit").pack(side="left")
        scalebar_entry = ctk.CTkEntry(
            scalebar_row, textvariable=self.scalebar_unit_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(scalebar_entry)
        scalebar_entry.pack(side="left", padx=(4, 12))
        self._body_label(scalebar_row, "vmin").pack(side="left")
        vmin_entry = ctk.CTkEntry(
            scalebar_row, textvariable=self.vmin_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(vmin_entry)
        vmin_entry.pack(side="left", padx=(4, 12))
        self._body_label(scalebar_row, "vmax").pack(side="left")
        vmax_entry = ctk.CTkEntry(
            scalebar_row, textvariable=self.vmax_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(vmax_entry)
        vmax_entry.pack(side="left", padx=(4, 0))
        self._option_row(
            self.advanced_content,
            25,
            "Scale & color range",
            widget=scalebar_row,
            hint="Scalebar unit label. vmin/vmax set continuous color range (blank = auto).",
        )

        rotations_row = self._make_sub_frame(self.advanced_content)
        rotations_entry = ctk.CTkEntry(
            rotations_row, textvariable=self.section_rotations_var, width=350, **self._theme["entry"],
        )
        self._register_entry_widget(rotations_entry)
        rotations_entry.pack(side="left")
        self._option_row(
            self.advanced_content,
            27,
            "Section rotations",
            widget=rotations_row,
            hint='JSON mapping section_id to degrees, e.g. {"sample_1": 90, "sample_2": 180}.',
        )

        self._set_advanced_visible(False)

        # ===================== Export Format tab =====================
        export_fmt_tab.columnconfigure(0, weight=1)
        export_fmt_tab.columnconfigure(1, weight=1)
        efmt_row = 0
        self._section_label(export_fmt_tab, "OUTPUT FORMAT").grid(row=efmt_row, column=0, columnspan=3, sticky="w", pady=(0, 6))
        efmt_row += 1
        self._subheader_label(
            export_fmt_tab,
            "Choose HTML (single file) or .karospace (sidecar package). Karospace requires sidecar gene storage.",
        ).grid(row=efmt_row, column=0, columnspan=3, sticky="w", pady=(0, 8))
        efmt_row += 1
        self._divider(export_fmt_tab, height=1).grid(row=efmt_row, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        efmt_row += 1

        output_fmt_row = self._make_sub_frame(export_fmt_tab)
        self.output_format_combo = ctk.CTkOptionMenu(
            output_fmt_row,
            variable=self.output_format_var,
            values=["embedded", "sidecar", "karospace", "sidecar + karospace"],
            **self._theme["combo"],
        )
        self._register_combo_widget(self.output_format_combo)
        self.output_format_combo.pack(side="left")
        efmt_row = self._option_row(
            export_fmt_tab,
            efmt_row,
            "Output format",
            widget=output_fmt_row,
            hint=(
                "embedded = single HTML. sidecar = HTML + gene shards (server hosting). "
                "karospace = .karospace + .loader.html (offline). "
                "sidecar + karospace = both outputs."
            ),
        )

        # --- Gene Storage & Encoding ---
        self._divider(export_fmt_tab, height=1).grid(row=efmt_row, column=0, columnspan=3, sticky="ew", pady=(6, 10))
        efmt_row += 1
        self._section_label(export_fmt_tab, "GENE STORAGE & ENCODING").grid(
            row=efmt_row, column=0, columnspan=3, sticky="w", pady=(0, 6),
        )
        efmt_row += 1

        storage_row = self._make_sub_frame(export_fmt_tab)
        self._body_label(storage_row, "Storage").pack(side="left")
        self.gene_storage_combo = ctk.CTkOptionMenu(
            storage_row,
            variable=self.gene_storage_var,
            values=["embedded", "sidecar"],
            **self._theme["combo"],
        )
        self._register_combo_widget(self.gene_storage_combo)
        self.gene_storage_combo.pack(side="left", padx=(4, 12))
        self._body_label(storage_row, "Encoding").pack(side="left")
        self.gene_encoding_combo = ctk.CTkOptionMenu(
            storage_row,
            variable=self.gene_encoding_var,
            values=["auto", "dense", "sparse"],
            **self._theme["combo"],
        )
        self._register_combo_widget(self.gene_encoding_combo)
        self.gene_encoding_combo.pack(side="left", padx=(4, 0))
        efmt_row = self._option_row(
            export_fmt_tab,
            efmt_row,
            "Gene storage",
            widget=storage_row,
            hint="embedded = genes in HTML. sidecar = external shards (required for .karospace).",
        )

        sidecar_row = self._make_sub_frame(export_fmt_tab)
        self._body_label(sidecar_row, "Value encoding").pack(side="left")
        self.gene_value_encoding_combo = ctk.CTkOptionMenu(
            sidecar_row,
            variable=self.gene_value_encoding_var,
            values=["float32", "uint16", "uint8"],
            **self._theme["combo"],
        )
        self._register_combo_widget(self.gene_value_encoding_combo)
        self.gene_value_encoding_combo.pack(side="left", padx=(4, 12))
        self._body_label(sidecar_row, "Sidecar format").pack(side="left")
        self.gene_sidecar_format_combo = ctk.CTkOptionMenu(
            sidecar_row,
            variable=self.gene_sidecar_format_var,
            values=["json-v2", "binary-v1"],
            **self._theme["combo"],
        )
        self._register_combo_widget(self.gene_sidecar_format_combo)
        self.gene_sidecar_format_combo.pack(side="left", padx=(4, 0))
        efmt_row = self._option_row(
            export_fmt_tab,
            efmt_row,
            "Value / sidecar fmt",
            widget=sidecar_row,
            hint="float32 = full precision. uint16/uint8 = quantized (requires binary-v1 format).",
        )

        shard_row = self._make_sub_frame(export_fmt_tab)
        self._body_label(shard_row, "Shard size").pack(side="left")
        shard_size_entry = ctk.CTkEntry(
            shard_row, textvariable=self.gene_sidecar_shard_size_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(shard_size_entry)
        shard_size_entry.pack(side="left", padx=(4, 12))
        self._body_label(shard_row, "Sparse threshold").pack(side="left")
        sparse_thresh_entry = ctk.CTkEntry(
            shard_row, textvariable=self.gene_sparse_zero_threshold_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(sparse_thresh_entry)
        sparse_thresh_entry.pack(side="left", padx=(4, 0))
        efmt_row = self._option_row(
            export_fmt_tab,
            efmt_row,
            "Shard / sparse",
            widget=shard_row,
            hint="Genes per shard file (sidecar). Sparse threshold = zero fraction for auto encoding.",
        )

        self.sidecar_aux_frame = self._make_sub_frame(export_fmt_tab)
        aux_entry = ctk.CTkEntry(
            self.sidecar_aux_frame, textvariable=self.gene_aux_path_var, width=300, **self._theme["entry"],
        )
        self._register_entry_widget(aux_entry)
        aux_entry.pack(side="left")
        self._body_label(self.sidecar_aux_frame, "(filename for .karospace, path otherwise)").pack(
            side="left", padx=(8, 0),
        )
        efmt_row = self._option_row(
            export_fmt_tab,
            efmt_row,
            "Gene aux path",
            widget=self.sidecar_aux_frame,
            hint="Sidecar gene manifest path. Leave blank for auto-generated name.",
        )

        pack_row = self._make_sub_frame(export_fmt_tab)
        self.pack_arrays_check = ctk.CTkCheckBox(
            pack_row, text="Pack arrays (base64)", variable=self.pack_arrays_var, **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.pack_arrays_check)
        self.pack_arrays_check.pack(side="left")
        self._body_label(pack_row, "Min length").pack(side="left", padx=(12, 4))
        pack_min_entry = ctk.CTkEntry(
            pack_row, textvariable=self.pack_arrays_min_len_var, width=70, **self._theme["entry"],
        )
        self._register_entry_widget(pack_min_entry)
        pack_min_entry.pack(side="left")
        efmt_row = self._option_row(
            export_fmt_tab,
            efmt_row,
            "Array packing",
            widget=pack_row,
            hint="Base64 packing for coordinates/colors. Min length = cells threshold for packing.",
        )

        help_tab.columnconfigure(0, weight=1)
        help_tab.rowconfigure(2, weight=1)
        self._section_label(help_tab, "GUIDE & WORKFLOW").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._divider(help_tab, height=1).grid(row=1, column=0, sticky="ew", pady=(0, 10))
        self.help_text = ctk.CTkTextbox(help_tab, wrap="word", **self._theme["textbox"])
        self.help_text.grid(row=2, column=0, sticky="nsew")
        self.help_text.insert(
            "1.0",
            "Basic tab\n"
            "- Input .h5ad: absolute path to your AnnData file.\n"
            "- Output directory: folder where KaroSpace_YYYYMMDD_HHMMSS.html is written.\n"
            "- Coordinates: auto/obsm/obs-centroid modes are converted to KaroSpace spatial input.\n"
            "- Section groupby: section split column used by load_spatial_data.\n"
            "- Initial color/theme/outline/title map directly to export_to_html.\n"
            "- Runtime mode: optional numba JIT performance mode (can be less stable in frozen app).\n"
            "- Downsample: integer cells per section (blank keeps all).\n"
            "- Metadata loading: group_order, metadata_columns, metadata_max_columns for load_spatial_data.\n\n"
            "Colors & Genes tab\n"
            "- additional_colors: obs columns offered as categorical coloring fields.\n"
            "- groupby lists: columns used for marker/neighbor/interaction/cluster DE analytics.\n"
            "- Tick selection mode: optional inspected multi-select mode.\n"
            "- genes mode:\n"
            "  hvgs -> use_hvgs=True with hvg_limit from count\n"
            "  top_mean -> compute top_mean genes list from adata\n"
            "  list_file -> choose a text file with one gene name per line\n"
            "  manual_list -> build list from var_names picker\n"
            "  (non-hvgs modes set use_hvgs=False)\n\n"
            "Advanced tab\n"
            "- Min panel size, spot size, marker top N.\n"
            "- Neighbor stats permutations/seed and auto groupby mode.\n"
            "- Interaction markers: limits, method (wilcoxon/t-test/logreg), layer.\n"
            "- Cluster DE: pairwise differential expression with volcano plots.\n"
            "- Gene analysis: correlation top N, cluster means N genes, spatial variable genes N.\n"
            "- Viewer rendering: scalebar unit, vmin/vmax color range, section rotations (JSON).\n"
            "- Serve after export starts a local preview server.\n\n"
            "Export Format tab\n"
            "- Output format:\n"
            "  embedded = single self-contained HTML file\n"
            "  sidecar = HTML + .genes.json + shard dir (for server hosting)\n"
            "  karospace = .karospace ZIP + .loader.html (for offline sharing)\n"
            "  sidecar + karospace = both sidecar and karospace outputs\n"
            "- Gene storage: embedded (in HTML) or sidecar (external shards).\n"
            "- Gene encoding: auto/dense/sparse. Value encoding: float32/uint16/uint8.\n"
            "- Sidecar format: json-v2 or binary-v1 (KSB1 compressed shards).\n"
            "- Shard size, sparse threshold, gene aux path.\n"
            "- Array packing: base64 optimization for coordinates/colors.\n"
            "- Note: .karospace output requires sidecar gene storage.\n\n"
            "Tip: click Inspect H5AD to load searchable dropdown choices from adata.obs and adata.var_names."
        )
        self.help_text.configure(state="disabled")

        self._divider(controls, height=1).grid(row=3, column=0, columnspan=3, sticky="ew", pady=(14, 0))

        button_row = self._make_sub_frame(controls)
        button_row.grid(row=4, column=0, columnspan=3, sticky="ew", pady=(12, 0))

        self.export_btn = self._primary_button(button_row, "Build Viewer", self._on_export, width=150)
        self.export_btn.pack(side="left")

        self.inspect_btn = self._secondary_button(button_row, "Inspect Dataset", self._inspect_h5ad, width=145)
        self.inspect_btn.pack(side="left", padx=(10, 0))

        self.stop_server_btn = self._secondary_button(button_row, "Stop Preview", self._stop_server, width=125)
        self.stop_server_btn.pack(side="left", padx=(10, 0))

        runtime_top = self._make_sub_frame(side)
        runtime_top.grid(row=0, column=0, sticky="ew")
        self._section_label(runtime_top, "STATUS").pack(side="left")
        self.runtime_chip_label = self._pill_label(runtime_top, text="READY", muted=True)
        self.runtime_chip_label.pack(side="right")

        self.progress = ctk.CTkProgressBar(side, mode="determinate", **self._theme["progress"])
        self.progress.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        self.progress.set(0.0)
        self._subheader_label(side, textvariable=self.status_var).grid(row=2, column=0, sticky="w", pady=(6, 14))

        launch_row = self._make_sub_frame(side)
        launch_row.grid(row=3, column=0, sticky="ew", pady=(0, 14))
        self._secondary_button(launch_row, "Open Folder", self._open_output_folder, width=130).pack(side="left")
        self._secondary_button(launch_row, "Open Viewer", self._open_viewer, width=120).pack(side="left", padx=(8, 0))

        self._divider(side, height=1).grid(row=4, column=0, sticky="ew", pady=(0, 14))

        self._section_label(side, "EVENT LOG").grid(row=5, column=0, sticky="w", pady=(0, 6))
        self.log_text = ctk.CTkTextbox(side, wrap="word", **self._theme["textbox"])
        self.log_text.grid(row=6, column=0, sticky="nsew")
        side.rowconfigure(6, weight=1)
        self.log_text.configure(state="disabled")

        self.genes_mode_var.trace_add("write", lambda *_: self._update_genes_mode_visibility())
        self.neighbor_auto_var.trace_add("write", lambda *_: self._update_neighbor_groupby_state())
        self.selection_mode_var.trace_add("write", lambda *_: self._update_selection_mode_visibility())
        self.status_var.trace_add("write", lambda *_: self._sync_runtime_chip())
        self.output_format_var.trace_add("write", lambda *_: self._update_output_format())
        self.gene_storage_var.trace_add("write", lambda *_: self._update_gene_storage_visibility())
        self._apply_preset("default", log=False)
        self._sync_app_theme_to_viewer_setting()
        self._sync_runtime_chip()
        self._update_genes_mode_visibility()
        self._update_neighbor_groupby_state()
        self._set_selection_mode_visible(False)

    def _path_field(
        self,
        parent: tk.Widget,
        row: int,
        label: str,
        variable: tk.StringVar,
        choose_file: bool,
        optional: bool = False,
    ) -> int:
        self._field_label(parent, label).grid(row=row, column=0, sticky="nw", pady=(0, 8), padx=(0, 14))
        entry = self._entry(parent, variable)
        entry.grid(row=row, column=1, sticky="ew", pady=(0, 8))

        if choose_file:
            button = self._secondary_button(parent, "Browse", lambda: self._choose_file(variable, optional=optional), width=96)
        else:
            button = self._secondary_button(parent, "Browse", lambda: self._choose_dir(variable), width=96)
        button.grid(row=row, column=2, sticky="e", pady=(0, 8), padx=(10, 0))
        return row + 1

    def _option_row(
        self,
        parent: tk.Widget,
        row: int,
        label: str,
        widget: tk.Widget,
        hint: str | None = None,
    ) -> int:
        self._field_label(parent, label).grid(row=row, column=0, sticky="nw", pady=(0, 2), padx=(0, 10))
        widget.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(0, 2))
        row += 1
        if hint:
            self._subheader_label(parent, hint).grid(row=row, column=1, columnspan=2, sticky="w", pady=(0, 6))
            row += 1
        return row

    def _entry(self, parent: tk.Widget, variable: tk.StringVar) -> ctk.CTkEntry:
        entry = ctk.CTkEntry(parent, textvariable=variable, **self._theme["entry"])
        self._register_entry_widget(entry)
        return entry

    def _coords_dropdown(self, parent: tk.Widget) -> ctk.CTkOptionMenu:
        combo = ctk.CTkOptionMenu(
            parent,
            variable=self.coords_var,
            values=["auto", "obsm:spatial", "obs:centroid_x_y"],
            **self._theme["combo"],
        )
        self._register_combo_widget(combo)
        return combo

    def _groupby_dropdown(self, parent: tk.Widget) -> ctk.CTkComboBox:
        self.groupby_combo = ctk.CTkComboBox(
            parent,
            variable=self.section_groupby_var,
            values=[],
            state="normal",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.groupby_combo)
        return self.groupby_combo

    def _color_dropdown(self, parent: tk.Widget) -> ctk.CTkComboBox:
        self.color_combo = ctk.CTkComboBox(
            parent,
            variable=self.initial_color_var,
            values=[],
            state="normal",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.color_combo)
        return self.color_combo

    def _outline_dropdown(self, parent: tk.Widget) -> ctk.CTkComboBox:
        self.outline_combo = ctk.CTkComboBox(
            parent,
            variable=self.outline_by_var,
            values=[],
            state="normal",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.outline_combo)
        return self.outline_combo

    def _theme_dropdown(self, parent: tk.Widget) -> ctk.CTkOptionMenu:
        combo = ctk.CTkOptionMenu(
            parent,
            variable=self.theme_var,
            values=["light", "dark"],
            command=self._on_theme_selected,
            **self._theme["combo"],
        )
        self._register_combo_widget(combo)
        return combo

    def _genes_mode_row(self, parent: tk.Widget) -> ctk.CTkFrame:
        wrap = self._make_sub_frame(parent)

        self.genes_mode_combo = ctk.CTkOptionMenu(
            wrap,
            variable=self.genes_mode_var,
            values=["hvgs", "top_mean", "list_file", "manual_list"],
            width=130,
            **self._theme["combo"],
        )
        self._register_combo_widget(self.genes_mode_combo)
        self.genes_mode_combo.pack(side="left")

        self.genes_count_entry = ctk.CTkEntry(wrap, textvariable=self.genes_count_var, width=92, **self._theme["entry"])
        self._register_entry_widget(self.genes_count_entry)
        self.genes_count_entry.pack(side="left", padx=(10, 0))

        self.genes_count_label = self._body_label(wrap, "N genes")
        self.genes_count_label.pack(side="left", padx=(8, 0))

        self.gene_list_entry = ctk.CTkEntry(wrap, textvariable=self.gene_list_path_var, width=260, **self._theme["entry"])
        self._register_entry_widget(self.gene_list_entry)
        self.gene_list_button = self._secondary_button(
            wrap,
            "Gene List",
            lambda: self._choose_file(self.gene_list_path_var, optional=False),
            width=96,
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

    def _update_output_format(self) -> None:
        fmt = self.output_format_var.get().strip().lower()
        needs_sidecar = fmt in {"sidecar", "karospace", "sidecar + karospace"}
        if needs_sidecar:
            self.gene_storage_var.set("sidecar")
            self.gene_storage_combo.configure(state="disabled")
        else:
            self.gene_storage_combo.configure(state="normal")

    def _update_gene_storage_visibility(self) -> None:
        storage = self.gene_storage_var.get().strip().lower()
        sidecar = storage == "sidecar"
        state = "normal" if sidecar else "disabled"
        self.gene_value_encoding_combo.configure(state=state)
        self.gene_sidecar_format_combo.configure(state=state)

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

    def _autofill_groupby_from_colors(self) -> None:
        primary = self.initial_color_var.get().strip()
        additional = self.additional_colors_editor.get_items()
        sources = ([primary] if primary else []) + additional
        if not sources:
            self._log("No colors to auto-fill from. Set initial color or add additional colors first.")
            return
        merged = self._merge_unique(sources)
        self.groupby_editor.set_items(merged)
        self._log(f"Auto-filled groupby with {len(merged)} columns: {', '.join(merged)}")

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
        self.numba_jit_var.set(False)
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
        self.interaction_markers_method_var.set("wilcoxon")
        self.interaction_markers_layer_var.set("normalized")

        # Output format defaults
        self.output_format_var.set("embedded")
        self.gene_encoding_var.set("auto")
        self.gene_value_encoding_var.set("float32")
        self.gene_sidecar_format_var.set("json-v2")
        self.gene_storage_var.set("embedded")
        self.gene_aux_path_var.set("")
        self.gene_sidecar_shard_size_var.set("256")
        self.gene_sparse_zero_threshold_var.set("0.8")
        self.pack_arrays_var.set(True)
        self.pack_arrays_min_len_var.set("1024")

        # Cluster DE defaults
        self.cluster_de_enabled_var.set(False)
        self.cluster_de_top_n_var.set("20")
        self.cluster_de_method_var.set("wilcoxon")
        self.cluster_de_layer_var.set("normalized")
        self.cluster_de_min_cells_var.set("20")

        # Gene analysis defaults
        self.gene_correlation_top_n_var.set("10")
        self.cluster_means_n_genes_var.set("500")
        self.spatial_variable_genes_n_var.set("200")

        # Viewer defaults
        self.section_rotations_var.set("")
        self.scalebar_unit_var.set("\u03bcm")
        self.vmin_var.set("")
        self.vmax_var.set("")
        self.viewer_info_html_var.set("")

        # load_spatial_data metadata defaults
        self.group_order_var.set("")
        self.metadata_columns_var.set("")
        self.metadata_value_order_var.set("")
        self.metadata_max_columns_var.set("")

        self.status_var.set("Ready")

        self._update_genes_mode_visibility()
        self._update_neighbor_groupby_state()
        self._update_gene_storage_visibility()
        self._load_inputs_into_tick_selection(log=False)
        if log:
            self._log("Applied default input values.")

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
        interaction_method = self.interaction_markers_method_var.get().strip() or "wilcoxon"
        interaction_layer_text = self.interaction_markers_layer_var.get().strip()
        interaction_layer = interaction_layer_text if interaction_layer_text else "normalized"

        # Output format
        output_format = self.output_format_var.get().strip().lower() or "embedded"
        if output_format not in {"embedded", "sidecar", "karospace", "sidecar + karospace"}:
            raise ValueError("Output format must be embedded, sidecar, karospace, or sidecar + karospace.")

        # Gene storage / encoding
        gene_storage = self.gene_storage_var.get().strip().lower() or "embedded"
        if output_format in {"sidecar", "karospace", "sidecar + karospace"} and gene_storage != "sidecar":
            gene_storage = "sidecar"
        gene_encoding = self.gene_encoding_var.get().strip().lower() or "auto"
        gene_value_encoding = self.gene_value_encoding_var.get().strip().lower() or "float32"
        gene_sidecar_format = self.gene_sidecar_format_var.get().strip() or "json-v2"
        gene_aux_path_text = self.gene_aux_path_var.get().strip() or None
        gene_sidecar_shard_size = self._parse_positive_int(
            "Sidecar shard size", self.gene_sidecar_shard_size_var.get() or "256",
        )
        try:
            gene_sparse_zero_threshold = float(self.gene_sparse_zero_threshold_var.get() or "0.8")
        except ValueError:
            raise ValueError("Gene sparse zero threshold must be a number between 0 and 1.")
        pack_arrays = bool(self.pack_arrays_var.get())
        pack_arrays_min_len = self._parse_positive_int(
            "Pack arrays min len", self.pack_arrays_min_len_var.get() or "1024",
        )

        # Cluster DE
        cluster_de_enabled = bool(self.cluster_de_enabled_var.get())
        cluster_de_groupby = (groupby_lists or None) if cluster_de_enabled else None
        cluster_de_top_n = self._parse_positive_int(
            "Cluster DE top N", self.cluster_de_top_n_var.get() or "20",
        )
        cluster_de_method = self.cluster_de_method_var.get().strip() or "wilcoxon"
        cluster_de_layer_text = self.cluster_de_layer_var.get().strip()
        cluster_de_layer = cluster_de_layer_text if cluster_de_layer_text else "normalized"
        cluster_de_min_cells = self._parse_positive_int(
            "Cluster DE min cells", self.cluster_de_min_cells_var.get() or "20",
        )

        # Gene analysis
        gene_correlation_top_n = self._parse_positive_int(
            "Gene correlation top N", self.gene_correlation_top_n_var.get() or "10",
        )
        cluster_means_n_genes = self._parse_positive_int(
            "Cluster means N genes", self.cluster_means_n_genes_var.get() or "500",
        )
        spatial_variable_genes_n = self._parse_positive_int(
            "Spatial variable genes N", self.spatial_variable_genes_n_var.get() or "200",
        )

        # Section rotations
        rotations_text = self.section_rotations_var.get().strip()
        section_rotations = None
        if rotations_text:
            try:
                section_rotations = json.loads(rotations_text)
                if not isinstance(section_rotations, dict):
                    raise ValueError("Section rotations must be a JSON object.")
            except json.JSONDecodeError as exc:
                raise ValueError(f"Section rotations must be valid JSON: {exc}") from exc

        # Scalebar, vmin/vmax
        scalebar_unit = self.scalebar_unit_var.get().strip() or "\u03bcm"
        vmin_text = self.vmin_var.get().strip()
        vmin = float(vmin_text) if vmin_text else None
        vmax_text = self.vmax_var.get().strip()
        vmax = float(vmax_text) if vmax_text else None
        viewer_info_html = self.viewer_info_html_var.get().strip() or None

        # load_spatial_data metadata params
        group_order_text = self.group_order_var.get().strip()
        group_order = [s.strip() for s in group_order_text.split(",") if s.strip()] if group_order_text else None
        metadata_columns_text = self.metadata_columns_var.get().strip()
        metadata_columns = [s.strip() for s in metadata_columns_text.split(",") if s.strip()] if metadata_columns_text else None
        metadata_max_columns_text = self.metadata_max_columns_var.get().strip()
        metadata_max_columns = (
            self._parse_positive_int("Metadata max columns", metadata_max_columns_text)
            if metadata_max_columns_text
            else None
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
            enable_numba_jit=bool(self.numba_jit_var.get()),
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
            interaction_markers_method=interaction_method,
            interaction_markers_layer=interaction_layer,
            output_format=output_format,
            gene_encoding=gene_encoding,
            gene_value_encoding=gene_value_encoding,
            gene_sidecar_format=gene_sidecar_format,
            gene_storage=gene_storage,
            gene_aux_path=gene_aux_path_text,
            gene_sidecar_shard_size=gene_sidecar_shard_size,
            gene_sparse_zero_threshold=gene_sparse_zero_threshold,
            pack_arrays=pack_arrays,
            pack_arrays_min_len=pack_arrays_min_len,
            cluster_de_enabled=cluster_de_enabled,
            cluster_de_groupby=cluster_de_groupby,
            cluster_de_top_n=cluster_de_top_n,
            cluster_de_method=cluster_de_method,
            cluster_de_layer=cluster_de_layer,
            cluster_de_min_cells=cluster_de_min_cells,
            gene_correlation_top_n=gene_correlation_top_n,
            cluster_means_n_genes=cluster_means_n_genes,
            spatial_variable_genes_n=spatial_variable_genes_n,
            section_rotations=section_rotations,
            scalebar_unit=scalebar_unit,
            vmin=vmin,
            vmax=vmax,
            viewer_info_html=viewer_info_html,
            group_order=group_order,
            metadata_columns=metadata_columns,
            metadata_value_order=None,
            metadata_max_columns=metadata_max_columns,
        )

    @staticmethod
    def _import_karospace_api(*, enable_numba_jit: bool = False):
        import importlib
        import importlib.metadata as importlib_metadata
        import inspect

        # scanpy imports numba with cache=True in some code paths; in certain
        # environments this crashes during import. Disable JIT by default for
        # robust GUI startup/export unless the user explicitly enables performance mode.
        os.environ["NUMBA_DISABLE_JIT"] = "0" if bool(enable_numba_jit) else "1"

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
            pkg_fn = getattr(module, "package_sidecar_viewer", None)
            return module.load_spatial_data, module.export_to_html, pkg_fn
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
                    pkg_fn = getattr(module, "package_sidecar_viewer", None)
                    return module.load_spatial_data, module.export_to_html, pkg_fn
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
            self.numba_jit_check,
            self.selection_mode_check,
            self.selection_apply_btn,
            self.selection_sync_btn,
            self.output_format_combo,
            self.gene_storage_combo,
            self.gene_encoding_combo,
            self.gene_value_encoding_combo,
            self.gene_sidecar_format_combo,
            self.pack_arrays_check,
            self.cluster_de_enabled_check,
            self.cluster_de_method_combo,
            self.interaction_method_combo,
            self.interaction_enabled_check,
        ]
        for widget in widgets:
            widget.configure(state="disabled" if busy else "normal")
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
        self.progress.set(percent / 100.0)
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
            f"interaction_enabled={config.interaction_markers_enabled}, "
            f"cluster_de_enabled={config.cluster_de_enabled}."
        )
        self._log(
            "Export format: "
            f"mode={config.output_format}, "
            f"gene_storage={config.gene_storage}, "
            f"gene_encoding={config.gene_encoding}, "
            f"value_encoding={config.gene_value_encoding}, "
            f"sidecar_format={config.gene_sidecar_format}."
        )
        self._log(
            "Gene analysis: "
            f"correlation_top_n={config.gene_correlation_top_n}, "
            f"spatial_var_genes={config.spatial_variable_genes_n}, "
            f"cluster_means={config.cluster_means_n_genes}."
        )
        self._log(
            "Runtime mode: "
            f"numba_jit={'on' if config.enable_numba_jit else 'off (safe mode)'}."
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
            load_spatial_data, export_to_html, package_sidecar_viewer = self._import_karospace_api(enable_numba_jit=config.enable_numba_jit)

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
            import inspect as _inspect_load
            _load_accepted = set(_inspect_load.signature(load_spatial_data).parameters)
            load_kwargs = {
                "groupby": config.section_groupby,
                "spatial_key": spatial_key,
                "group_order": config.group_order,
                "metadata_columns": config.metadata_columns,
                "metadata_value_order": config.metadata_value_order,
                "metadata_max_columns": config.metadata_max_columns,
            }
            load_kwargs = {k: v for k, v in load_kwargs.items() if k in _load_accepted}
            dataset = load_spatial_data(str(input_path), **load_kwargs)
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

            # Determine export mode.
            # "karospace" mode: export_to_html handles .karospace directly.
            # "sidecar + karospace": export as sidecar first, then package.
            # "sidecar" / "embedded": export as .html with gene_storage accordingly.
            fmt = config.output_format
            export_as_karospace_directly = fmt == "karospace"
            package_after_sidecar = fmt == "sidecar + karospace"

            if export_as_karospace_directly:
                output_path = self._build_output_html_path(config.outdir, output_format="karospace")
            else:
                output_path = self._build_output_html_path(config.outdir, output_format="html")

            fmt_label = {
                "embedded": "embedded HTML viewer",
                "sidecar": "sidecar viewer",
                "karospace": "karospace package",
                "sidecar + karospace": "sidecar viewer (+ karospace packaging)",
            }.get(fmt, "viewer")
            emit_progress(76, "Writing viewer", f"Exporting {fmt_label} to {output_path}.")
            export_started = time.perf_counter()

            # Build kwargs dict, then filter to only params the installed
            # karospace version accepts. This keeps Builder compatible with
            # older karospace releases that lack newer parameters.
            import inspect as _inspect
            _accepted = set(_inspect.signature(export_to_html).parameters)

            export_kwargs = {
                "output_path": str(output_path),
                "color": config.initial_color,
                "title": config.title,
                "min_panel_size": config.min_panel_size,
                "spot_size": config.spot_size,
                "downsample": config.downsample,
                "theme": config.theme,
                "outline_by": config.outline_by,
                "additional_colors": config.additional_colors,
                "genes": config.genes,
                "use_hvgs": config.use_hvgs,
                "hvg_limit": config.hvg_limit,
                "marker_genes_groupby": marker_groupby,
                "marker_genes_top_n": config.marker_genes_top_n,
                "neighbor_stats_groupby": neighbor_groupby,
                "neighbor_stats_permutations": neighbor_permutations,
                "neighbor_stats_seed": config.neighbor_stats_seed,
                "interaction_markers_groupby": interaction_groupby if config.interaction_markers_enabled else None,
                "interaction_markers_top_targets": config.interaction_markers_top_targets,
                "interaction_markers_top_genes": config.interaction_markers_top_genes,
                "interaction_markers_min_cells": config.interaction_markers_min_cells,
                "interaction_markers_min_neighbors": config.interaction_markers_min_neighbors,
                "interaction_markers_method": config.interaction_markers_method,
                "interaction_markers_layer": config.interaction_markers_layer,
                "gene_encoding": config.gene_encoding,
                "gene_value_encoding": config.gene_value_encoding,
                "gene_sidecar_format": config.gene_sidecar_format,
                "gene_storage": config.gene_storage,
                "gene_aux_path": config.gene_aux_path,
                "gene_sidecar_shard_size": config.gene_sidecar_shard_size,
                "gene_sparse_zero_threshold": config.gene_sparse_zero_threshold,
                "pack_arrays": config.pack_arrays,
                "pack_arrays_min_len": config.pack_arrays_min_len,
                "cluster_de_groupby": config.cluster_de_groupby,
                "cluster_de_top_n": config.cluster_de_top_n,
                "cluster_de_method": config.cluster_de_method,
                "cluster_de_layer": config.cluster_de_layer,
                "cluster_de_min_cells": config.cluster_de_min_cells,
                "gene_correlation_top_n": config.gene_correlation_top_n,
                "cluster_means_n_genes": config.cluster_means_n_genes,
                "spatial_variable_genes_n": config.spatial_variable_genes_n,
                "section_rotations": config.section_rotations,
                "scalebar_unit": config.scalebar_unit,
                "vmin": config.vmin,
                "vmax": config.vmax,
                "viewer_info_html": config.viewer_info_html,
            }
            skipped = sorted(set(export_kwargs) - _accepted)
            if skipped:
                self._queue.put((
                    "log",
                    f"Note: installed karospace does not support: {', '.join(skipped)}. "
                    "Update karospace to use these features.",
                ))
            filtered_kwargs = {k: v for k, v in export_kwargs.items() if k in _accepted}

            output_html_path = Path(
                export_to_html(dataset, **filtered_kwargs)
            ).expanduser()
            export_elapsed = time.perf_counter() - export_started
            emit_progress(90, "Finalizing", f"Viewer exported in {export_elapsed:.1f}s: {output_html_path}")

            # If "sidecar + karospace", package the sidecar output into .karospace
            if package_after_sidecar:
                if package_sidecar_viewer is None:
                    self._queue.put((
                        "log",
                        "Warning: package_sidecar_viewer not available in this karospace version. "
                        "Sidecar files exported, but .karospace packaging skipped. "
                        "Update karospace to enable packaging.",
                    ))
                else:
                    emit_progress(92, "Packaging", "Creating .karospace archive from sidecar output.")
                    package_started = time.perf_counter()
                    karospace_path = package_sidecar_viewer(str(output_html_path))
                    package_elapsed = time.perf_counter() - package_started
                    self._queue.put((
                        "log",
                        f"Packaged .karospace + .loader.html in {package_elapsed:.1f}s: {karospace_path}",
                    ))

            result = AppResult(
                outdir=output_html_path.parent,
                n_cells=int(dataset.n_cells),
                n_sections=int(dataset.n_sections),
                output_html=output_html_path,
            )
            total_elapsed = time.perf_counter() - total_started
            emit_progress(100, "Complete", f"Total export time: {total_elapsed:.1f}s.")
            self._queue.put(("done", result))
            if serve_after_export and output_html_path.suffix.lower() == ".html":
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

    def _build_output_html_path(self, outdir: Path, output_format: str = "html") -> Path:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base = f"{self._OUTPUT_HTML_BASENAME}_{stamp}"
        suffix = ".karospace" if output_format == "karospace" else ".html"
        candidate = outdir / f"{base}{suffix}"
        counter = 1
        while candidate.exists():
            candidate = outdir / f"{base}_{counter:02d}{suffix}"
            counter += 1
        return candidate

    def _resolve_viewer_html(self, outdir: Path) -> Path:
        for ext in ("*.html", "*.karospace"):
            pattern = f"{self._OUTPUT_HTML_BASENAME}_{ext}"
            candidates = sorted(outdir.glob(pattern))
            if candidates:
                return candidates[-1]
        legacy_named = outdir / f"{self._OUTPUT_HTML_BASENAME}.html"
        if legacy_named.exists():
            return legacy_named
        legacy = outdir / "index.html"
        if legacy.exists():
            return legacy
        return outdir / f"{self._OUTPUT_HTML_BASENAME}_YYYYMMDD_HHMMSS.html"

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
