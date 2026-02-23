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


def _ctk_theme_config(palette: dict[str, str]) -> dict[str, dict[str, object]]:
    return {
        "root": {"fg_color": palette["background"]},
        "root_frame": {"fg_color": palette["background"], "corner_radius": 0},
        "card_frame": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 14,
            "border_width": 1,
            "border_color": palette["border"],
        },
        "hero_card": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 18,
            "border_width": 1,
            "border_color": palette["accent"],
        },
        "highlight_card": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 12,
            "border_width": 1,
            "border_color": palette["accent_strong"],
        },
        "sub_frame": {"fg_color": "transparent", "corner_radius": 0},
        "section_label": {
            "font": ("Segoe UI", 10, "bold"),
            "text_color": palette["accent"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "hero_label": {
            "font": ("Segoe UI", 32, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "header_label": {
            "font": ("Segoe UI", 20, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "subheader_label": {
            "font": ("Segoe UI", 11),
            "text_color": palette["muted"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "field_label": {
            "font": ("Segoe UI", 11, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "body_label": {
            "font": ("Segoe UI", 10),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "primary_button": {
            "fg_color": palette["accent"],
            "hover_color": palette["accent_strong"],
            "text_color": palette.get("on_accent", "#ffffff"),
            "corner_radius": 12,
            "font": ("Segoe UI", 11, "bold"),
            "height": 38,
            "border_width": 0,
        },
        "secondary_button": {
            "fg_color": palette["input_bg"],
            "hover_color": palette["hover_bg"],
            "text_color": palette["text"],
            "corner_radius": 12,
            "font": ("Segoe UI", 10),
            "height": 36,
            "border_width": 1,
            "border_color": palette["border"],
        },
        "pill_label": {
            "font": ("Segoe UI", 9, "bold"),
            "text_color": palette.get("on_accent", "#ffffff"),
            "fg_color": palette["accent"],
            "corner_radius": 999,
            "padx": 10,
            "pady": 4,
        },
        "muted_pill_label": {
            "font": ("Segoe UI", 9, "bold"),
            "text_color": palette["text"],
            "fg_color": palette["hover_bg"],
            "corner_radius": 999,
            "padx": 10,
            "pady": 4,
        },
        "entry": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "placeholder_text_color": palette["muted"],
            "border_color": palette["border"],
            "corner_radius": 10,
        },
        "combo": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "button_color": palette["hover_bg"],
            "button_hover_color": palette["accent"],
            "dropdown_fg_color": palette["panel_bg"],
            "dropdown_text_color": palette["text"],
            "dropdown_hover_color": palette["hover_bg"],
            "corner_radius": 10,
        },
        "checkbox": {
            "fg_color": palette["accent"],
            "hover_color": palette["accent_strong"],
            "checkmark_color": palette.get("on_accent", "#ffffff"),
            "text_color": palette["text"],
            "border_color": palette["border"],
            "font": ("Segoe UI", 10),
        },
        "tabview": {
            "fg_color": palette["panel_bg"],
            "segmented_button_fg_color": palette["hover_bg"],
            "segmented_button_selected_color": palette["accent"],
            "segmented_button_selected_hover_color": palette["accent_strong"],
            "segmented_button_unselected_color": palette["hover_bg"],
            "segmented_button_unselected_hover_color": palette["border"],
            "text_color": palette["text"],
            "corner_radius": 12,
            "border_width": 0,
            "border_color": palette["border"],
        },
        "divider": {
            "fg_color": palette["border"],
            "corner_radius": 999,
        },
        "textbox": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "border_color": palette["border"],
            "corner_radius": 10,
            "border_width": 1,
            "scrollbar_button_color": palette["hover_bg"],
            "scrollbar_button_hover_color": palette["accent"],
        },
        "progress": {
            "fg_color": palette["hover_bg"],
            "progress_color": palette["accent"],
            "border_color": palette["border"],
            "corner_radius": 10,
            "height": 14,
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
        self.label_widget.grid(row=0, column=0, sticky="w", pady=(0, 6))

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
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=height,
            selectmode="extended",
            activestyle="none",
            relief="solid",
            bd=1,
            highlightthickness=1,
            font=("Segoe UI", 10),
        )
        self.listbox.grid(row=0, column=0, sticky="ew")
        self.scroll = tk.Scrollbar(list_wrap, orient="vertical", command=self.listbox.yview)
        self.scroll.grid(row=0, column=1, sticky="ns")
        self.listbox.configure(yscrollcommand=self.scroll.set)

        self.help_label: ctk.CTkLabel | None = None
        if help_text:
            self.help_label = ctk.CTkLabel(self, text=help_text, **self._theme["subheader_label"])
            self.help_label.grid(row=3, column=0, sticky="w", pady=(6, 0))

        self.apply_palette(self._palette)

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
        self.label_widget.grid(row=0, column=0, sticky="w", pady=(0, 6))

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
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(8, 0))
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=height,
            selectmode="extended",
            activestyle="none",
            relief="solid",
            bd=1,
            highlightthickness=1,
            font=("Segoe UI", 10),
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


class _ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


class ExportApp(ctk.CTk if ctk is not None else object):
    _OUTPUT_HTML_BASENAME = "KaroSpace"
    _NEIGHBOR_MATRIX_BUDGET_BYTES = 512 * 1024 * 1024
    _MARKER_GROUPBY_MAX_UNIQUE = 96
    _INTERACTION_GROUPBY_MAX_UNIQUE = 96
    _QUALITATIVE_MAX_UNIQUE = 96
    _QUALITATIVE_MAX_FRACTION = 0.25

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
            try:
                progress_value = float(self.progress.get())
            except Exception:
                progress_value = 0.0
            if progress_value <= 0.0:
                self.progress.configure(progress_color=palette["hover_bg"])
        if hasattr(self, "main_scroll_frame"):
            self.main_scroll_frame.configure(
                scrollbar_button_color=palette["hover_bg"],
                scrollbar_button_hover_color=palette["accent"],
            )
        self._sync_theme_toggle_icon()
        self._update_action_states()
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

    def _toggle_theme(self) -> None:
        current = self.theme_var.get().strip().lower() or "dark"
        next_mode = "light" if current == "dark" else "dark"
        self._on_theme_selected(next_mode)

    def _sync_theme_toggle_icon(self) -> None:
        if not hasattr(self, "theme_toggle_btn"):
            return
        mode = self.theme_var.get().strip().lower() or "dark"
        if mode == "dark":
            icon, label = "☀", "Light mode"
        else:
            icon, label = "☾", "Dark mode"
        self.theme_toggle_btn.configure(text=icon)
        if hasattr(self, "theme_toggle_label"):
            self.theme_toggle_label.configure(text=label)

    def _set_button_variant(self, button: ctk.CTkButton, *, primary: bool, enabled: bool) -> None:
        role = "primary_button" if primary else "secondary_button"
        button.configure(**self._theme[role], state="normal" if enabled else "disabled")

    def _resolve_h5ad_input_path(self) -> Path | None:
        raw = self.h5ad_var.get().strip()
        if not raw:
            return None
        try:
            return Path(raw).expanduser().resolve()
        except Exception:
            return None

    def _is_current_input_inspected(self) -> bool:
        current = self._resolve_h5ad_input_path()
        if current is None or self._inspected_h5ad_path is None:
            return False
        return current == self._inspected_h5ad_path

    def _set_dataset_controls_visible(self, inspected_ready: bool) -> None:
        if hasattr(self, "colors_locked_hint"):
            if inspected_ready:
                self.colors_locked_hint.grid_remove()
            else:
                self.colors_locked_hint.grid()
        if hasattr(self, "colors_dataset_frame"):
            if inspected_ready:
                self.colors_dataset_frame.grid()
            else:
                self.colors_dataset_frame.grid_remove()

        for attr in ("groupby_combo", "color_combo", "outline_combo"):
            widget = getattr(self, attr, None)
            if widget is not None:
                widget.configure(state="normal" if inspected_ready else "disabled")

    def _update_action_states(self) -> None:
        if not hasattr(self, "inspect_btn") or not hasattr(self, "export_btn"):
            return
        busy = bool(self._export_thread and self._export_thread.is_alive())
        has_input = bool(self.h5ad_var.get().strip())
        inspected_ready = self._is_current_input_inspected()

        inspect_enabled = has_input and not busy
        export_enabled = inspected_ready and not busy

        self._set_button_variant(self.inspect_btn, primary=has_input, enabled=inspect_enabled)
        self._set_button_variant(self.export_btn, primary=inspected_ready, enabled=export_enabled)
        self._set_dataset_controls_visible(inspected_ready)

    def _on_h5ad_input_changed(self, *_args) -> None:
        if not self._is_current_input_inspected() and self.status_var.get().strip().lower() == "inspection complete":
            self.status_var.set("Ready")
        self._update_action_states()

    @staticmethod
    def _widget_descends_from(widget: tk.Widget | None, ancestor: tk.Widget | None) -> bool:
        if widget is None or ancestor is None:
            return False
        target = str(ancestor)
        current = widget
        while current is not None:
            if str(current) == target:
                return True
            parent_name = str(current.winfo_parent())
            if not parent_name:
                break
            try:
                current = current.nametowidget(parent_name)
            except Exception:
                break
        return False

    @staticmethod
    def _widget_handles_own_wheel(widget: tk.Widget | None) -> bool:
        if widget is None:
            return False
        class_name = str(widget.winfo_class()).lower()
        return any(token in class_name for token in ("listbox", "text", "entry", "scrollbar"))

    def _install_main_scroll_wheel_bindings(self) -> None:
        if not hasattr(self, "main_scroll_frame"):
            return
        canvas = getattr(self.main_scroll_frame, "_parent_canvas", None)
        if canvas is None:
            return
        canvas.configure(yscrollincrement=18)

        def _on_wheel(event) -> str | None:
            widget = getattr(event, "widget", None)
            if not self._widget_descends_from(widget, self.main_scroll_frame):
                return None
            if self._widget_handles_own_wheel(widget):
                return None
            delta = getattr(event, "delta", 0)
            if delta:
                step = -1 if int(delta) > 0 else 1
                canvas.yview_scroll(step, "units")
                return "break"
            num = getattr(event, "num", 0)
            if int(num) == 4:
                canvas.yview_scroll(-1, "units")
                return "break"
            if int(num) == 5:
                canvas.yview_scroll(1, "units")
                return "break"
            return None

        self.bind_all("<MouseWheel>", _on_wheel, add="+")
        self.bind_all("<Button-4>", _on_wheel, add="+")
        self.bind_all("<Button-5>", _on_wheel, add="+")

    def _qualitative_obs_columns(self, obs) -> tuple[list[str], dict[str, int]]:
        try:
            from pandas.api.types import (
                is_bool_dtype,
                is_categorical_dtype,
                is_integer_dtype,
                is_object_dtype,
                is_string_dtype,
            )
        except Exception:
            return [str(c) for c in getattr(obs, "columns", [])], {}

        columns = [str(c) for c in getattr(obs, "columns", [])]
        n_obs = int(getattr(obs, "shape", [0])[0]) if getattr(obs, "shape", None) else 0
        qualitative: list[str] = []
        counts: dict[str, int] = {}
        max_ratio_count = max(2, int(round(n_obs * self._QUALITATIVE_MAX_FRACTION))) if n_obs else 2

        for column in columns:
            try:
                series = obs[column]
            except Exception:
                continue

            try:
                unique_count = int(series.nunique(dropna=True))
            except Exception:
                unique_count = -1
            if unique_count >= 0:
                counts[column] = unique_count
                if unique_count < 2:
                    continue
                if unique_count > self._QUALITATIVE_MAX_UNIQUE:
                    continue
                if unique_count > max_ratio_count and unique_count > 24:
                    continue

            dtype = getattr(series, "dtype", None)
            is_qualitative = (
                is_categorical_dtype(dtype)
                or is_string_dtype(dtype)
                or is_object_dtype(dtype)
                or is_bool_dtype(dtype)
                or (is_integer_dtype(dtype) and (unique_count < 0 or unique_count <= 24))
            )
            if is_qualitative:
                qualitative.append(column)

        return qualitative, counts

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
        self.selection_mode_var = tk.BooleanVar(value=False)

        self.downsample_var = tk.StringVar()

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
        root.grid(row=0, column=0, sticky="nsew", padx=12, pady=12)
        root.columnconfigure(0, weight=3)
        root.columnconfigure(1, weight=2)
        root.rowconfigure(0, weight=1)

        controls = ctk.CTkFrame(root, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", controls)
        controls.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        side = ctk.CTkFrame(root, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", side)
        side.grid(row=0, column=1, sticky="nsew")

        controls_inner = ctk.CTkFrame(controls, **self._theme["sub_frame"])
        self._register_theme_widget("sub_frame", controls_inner)
        controls_inner.pack(fill="both", expand=True, padx=18, pady=18)
        controls = controls_inner

        side_inner = ctk.CTkFrame(side, **self._theme["sub_frame"])
        self._register_theme_widget("sub_frame", side_inner)
        side_inner.pack(fill="both", expand=True, padx=18, pady=18)
        side = side_inner

        controls.columnconfigure(1, weight=1)
        side.columnconfigure(0, weight=1)
        side.rowconfigure(6, weight=1)

        hero = ctk.CTkFrame(controls, **self._theme["hero_card"])
        self._register_theme_widget("hero_card", hero)
        hero.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        hero.columnconfigure(0, weight=1)

        hero_inner = self._make_sub_frame(hero)
        hero_inner.grid(row=0, column=0, sticky="ew", padx=14, pady=14)
        hero_inner.columnconfigure(0, weight=1)

        hero_left = self._make_sub_frame(hero_inner)
        hero_left.grid(row=0, column=0, sticky="w")
        self._section_label(hero_left, "DESKTOP BUILDER").pack(anchor="w")
        self._hero_label(hero_left, "KaroSpaceBuilder").pack(anchor="w", pady=(2, 0))
        self._subheader_label(
            hero_left,
            "Export AnnData into a static KaroSpace viewer bundle with guided inputs and inspected field pickers.",
        ).pack(anchor="w", pady=(2, 0))

        self._divider(controls, height=1).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 12))

        notebook_wrap = self._make_sub_frame(controls)
        notebook_wrap.grid(row=2, column=0, columnspan=3, sticky="nsew")
        notebook_wrap.columnconfigure(0, weight=1)
        notebook_wrap.rowconfigure(0, weight=1)
        controls.rowconfigure(2, weight=1)

        notebook = ctk.CTkTabview(notebook_wrap, **self._theme["tabview"])
        self._register_theme_widget("tabview", notebook)
        notebook.grid(row=0, column=0, sticky="nsew")

        theme_toggle_wrap = self._make_sub_frame(notebook_wrap)
        theme_toggle_wrap.grid(row=0, column=1, sticky="ne", padx=(8, 0))
        self.theme_toggle_btn = self._secondary_button(theme_toggle_wrap, "☀", self._toggle_theme, width=34)
        self.theme_toggle_btn.pack(side="top", anchor="e")
        self.theme_toggle_label = self._subheader_label(theme_toggle_wrap, "Light mode")
        self.theme_toggle_label.pack(side="top", anchor="e", pady=(4, 0))

        notebook.add("Basic")
        notebook.add("Colors & Genes")
        notebook.add("Advanced")
        notebook.add("Help")
        basic_tab = notebook.tab("Basic")
        colors_tab = notebook.tab("Colors & Genes")
        advanced_tab = notebook.tab("Advanced")
        help_tab = notebook.tab("Help")
        for tab in (basic_tab, colors_tab, advanced_tab, help_tab):
            self._register_theme_widget("sub_frame", tab)

        basic_tab.columnconfigure(1, weight=1)
        self._section_label(basic_tab, "Input & Viewer Setup").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        self._divider(basic_tab, height=1).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 10))
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
            "Performance mode",
            widget=runtime_mode_row,
            hint="Optional speed-up for heavy datasets. Turn off if exports become unstable.",
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
        self._option_row(
            basic_tab,
            row,
            "Downsample",
            widget=downsample_container,
            hint="Maps to export_to_html(downsample=...).",
        )

        colors_tab.columnconfigure(0, weight=1)
        self._section_label(colors_tab, "Color Fields & Gene Sources").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._subheader_label(
            colors_tab,
            "Build color menus and gene features from inspected obs/var fields.",
        ).grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._divider(colors_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        self.colors_locked_hint = self._subheader_label(
            colors_tab,
            "Inspect Dataset first to unlock color, groupby, and gene parameters from this .h5ad.",
        )
        self.colors_locked_hint.grid(row=3, column=0, sticky="w", pady=(0, 8))

        self.colors_dataset_frame = self._make_sub_frame(colors_tab)
        self.colors_dataset_frame.grid(row=4, column=0, sticky="ew")
        self.colors_dataset_frame.columnconfigure(0, weight=1)
        self.additional_colors_editor = SearchableListEditor(
            self.colors_dataset_frame,
            label="additional_colors (obs columns)",
            height=6,
            help_text="These become color options in KaroSpace. Use Inspect to load obs columns.",
            palette=self._app_palette,
        )
        self.additional_colors_editor.grid(row=0, column=0, sticky="ew", pady=(0, 10))

        self.groupby_editor = SearchableListEditor(
            self.colors_dataset_frame,
            label="groupby lists (marker/neighbor/interaction)",
            height=6,
            help_text="Used for marker_genes_groupby and optionally neighbor/interaction groupby.",
            palette=self._app_palette,
        )
        self.groupby_editor.grid(row=1, column=0, sticky="ew", pady=(0, 12))

        genes_card = ctk.CTkFrame(self.colors_dataset_frame, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", genes_card)
        genes_card.grid(row=2, column=0, sticky="ew")
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
            height=8,
            help_text="Search var_names and build the genes list with + Add / Remove.",
            palette=self._app_palette,
        )
        self.manual_genes_editor.grid(row=3, column=0, sticky="ew")

        selection_card = ctk.CTkFrame(self.colors_dataset_frame, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", selection_card)
        selection_card.grid(row=3, column=0, sticky="ew", pady=(12, 0))
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
            height=7,
            help_text="Multi-select obs fields to include as additional viewer colors.",
            palette=self._app_palette,
        )
        self.selection_additional_picker.grid(row=0, column=0, sticky="ew", padx=(0, 8))

        self.selection_groupby_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick groupby lists (obs)",
            height=7,
            help_text="Multi-select obs columns for marker/neighbor/interaction groupby lists.",
            palette=self._app_palette,
        )
        self.selection_groupby_picker.grid(row=0, column=1, sticky="ew")

        self.selection_genes_picker = SearchableMultiSelectEditor(
            self.selection_mode_content,
            label="Tick genes (manual_list mode)",
            height=8,
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
        self._section_label(advanced_tab, "Advanced Analytics").grid(row=0, column=0, sticky="w", pady=(0, 6))
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
            11,
            "Preview server",
            widget=serve_row,
            hint="Optional local server to open the latest generated KaroSpace_*.html file.",
        )
        self._set_advanced_visible(False)

        help_tab.columnconfigure(0, weight=1)
        help_tab.rowconfigure(2, weight=1)
        self._section_label(help_tab, "Guide & Workflow").grid(row=0, column=0, sticky="w", pady=(0, 6))
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
            "- Initial color/outline/title map directly to export_to_html.\n"
            "- Theme toggle: sun/moon button beside tabs switches exported viewer and app theme.\n"
            "- Performance mode: optional numba JIT speed-up (can be less stable in frozen app).\n"
            "- Downsample: integer cells per section (blank keeps all).\n\n"
            "Colors & Genes tab\n"
            "- Locked until Inspect H5AD is complete for the selected input file.\n"
            "- additional_colors: qualitative obs columns offered as categorical coloring fields.\n"
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
            "Tip: click Inspect H5AD to load searchable dropdown choices from adata.obs and adata.var_names."
        )
        self.help_text.configure(state="disabled")

        action_card = ctk.CTkFrame(controls, **self._theme["highlight_card"])
        self._register_theme_widget("highlight_card", action_card)
        action_card.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(12, 0))
        button_row = self._make_sub_frame(action_card)
        button_row.pack(fill="x", padx=12, pady=10)
        self._section_label(button_row, "Actions").pack(side="left", padx=(0, 10))

        self.inspect_btn = self._secondary_button(button_row, "Inspect Dataset", self._inspect_h5ad, width=150)
        self.inspect_btn.pack(side="left")

        self.export_btn = self._secondary_button(button_row, "Build Viewer", self._on_export, width=140)
        self.export_btn.pack(side="left", padx=(10, 0))

        self.stop_server_btn = self._secondary_button(button_row, "Stop Preview", self._stop_server, width=130)
        self.stop_server_btn.pack(side="left", padx=(10, 0))

        runtime_top = self._make_sub_frame(side)
        runtime_top.grid(row=0, column=0, sticky="ew")
        self._section_label(runtime_top, "Runtime").pack(side="left")
        self.runtime_chip_label = self._pill_label(runtime_top, text="READY", muted=True)
        self.runtime_chip_label.pack(side="right")
        self._header_label(side, "Activity").grid(row=1, column=0, sticky="w", pady=(4, 0))
        self._subheader_label(side, textvariable=self.status_var).grid(row=2, column=0, sticky="w", pady=(2, 12))

        self.progress = ctk.CTkProgressBar(side, mode="determinate", **self._theme["progress"])
        self.progress.grid(row=3, column=0, sticky="ew")
        self.progress.set(0.0)
        self._subheader_label(side, "Live export progress").grid(row=4, column=0, sticky="w", pady=(6, 10))

        launch_card = ctk.CTkFrame(side, **self._theme["highlight_card"])
        self._register_theme_widget("highlight_card", launch_card)
        launch_card.grid(row=5, column=0, sticky="ew", pady=(0, 10))
        launch_row = self._make_sub_frame(launch_card)
        launch_row.pack(fill="x", padx=10, pady=10)
        self._secondary_button(launch_row, "Open Output Folder", self._open_output_folder, width=160).pack(side="left")
        self._secondary_button(launch_row, "Open Viewer", self._open_viewer, width=120).pack(side="left", padx=(10, 0))

        log_card = ctk.CTkFrame(side, **self._theme["highlight_card"])
        self._register_theme_widget("highlight_card", log_card)
        log_card.grid(row=6, column=0, sticky="nsew")
        log_wrap = self._make_sub_frame(log_card)
        log_wrap.pack(fill="both", expand=True, padx=10, pady=10)
        log_wrap.columnconfigure(0, weight=1)
        log_wrap.rowconfigure(0, weight=1)

        self.log_text = ctk.CTkTextbox(log_wrap, wrap="word", font=("Consolas", 10), **self._theme["textbox"])
        self.log_text.grid(row=0, column=0, sticky="nsew")
        self.log_text.configure(state="disabled")

        self.genes_mode_var.trace_add("write", lambda *_: self._update_genes_mode_visibility())
        self.neighbor_auto_var.trace_add("write", lambda *_: self._update_neighbor_groupby_state())
        self.selection_mode_var.trace_add("write", lambda *_: self._update_selection_mode_visibility())
        self.status_var.trace_add("write", lambda *_: self._sync_runtime_chip())
        self.h5ad_var.trace_add("write", self._on_h5ad_input_changed)
        self._apply_preset("default", log=False)
        self._install_main_scroll_wheel_bindings()
        self._sync_app_theme_to_viewer_setting()
        self._sync_theme_toggle_icon()
        self._sync_runtime_chip()
        self._update_genes_mode_visibility()
        self._update_neighbor_groupby_state()
        self._set_selection_mode_visible(False)
        self._update_action_states()

    def _path_field(
        self,
        parent: tk.Widget,
        row: int,
        label: str,
        variable: tk.StringVar,
        choose_file: bool,
        optional: bool = False,
    ) -> int:
        self._field_label(parent, label).grid(row=row, column=0, sticky="nw", pady=(0, 6), padx=(0, 12))
        entry = self._entry(parent, variable)
        entry.grid(row=row, column=1, sticky="ew", pady=(0, 6))

        if choose_file:
            button = self._secondary_button(parent, "Browse", lambda: self._choose_file(variable, optional=optional), width=96)
        else:
            button = self._secondary_button(parent, "Browse", lambda: self._choose_dir(variable), width=96)
        button.grid(row=row, column=2, sticky="e", pady=(0, 6), padx=(10, 0))
        return row + 1

    def _option_row(
        self,
        parent: tk.Widget,
        row: int,
        label: str,
        widget: tk.Widget,
        hint: str | None = None,
    ) -> int:
        self._field_label(parent, label).grid(row=row, column=0, sticky="nw", pady=(0, 2), padx=(0, 12))
        widget.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(0, 2))
        row += 1
        if hint:
            self._subheader_label(parent, hint).grid(row=row, column=1, columnspan=2, sticky="w", pady=(0, 8))
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
            self.h5ad_var.set("")
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
        self.status_var.set("Ready")

        self._update_genes_mode_visibility()
        self._update_neighbor_groupby_state()
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
        path = self._resolve_h5ad_input_path()
        if path is None:
            messagebox.showerror("Missing input", "Pick an input .h5ad first.")
            return

        if not path.exists():
            messagebox.showerror("Missing file", f"Input file not found:\n{path}")
            return

        self._inspected_h5ad_path = None
        self._inspected_var_name_set = None
        self._inspected_coords_mode = None
        self._update_action_states()
        self._log(f"Inspecting {path}")
        adata = None
        try:
            ad_mod = _get_anndata()
            try:
                adata = ad_mod.read_h5ad(path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(path)

            all_obs_cols = [str(c) for c in adata.obs.columns]
            obs_cols, unique_counts = self._qualitative_obs_columns(adata.obs)
            if not obs_cols:
                raise ValueError(
                    "No qualitative obs columns were detected after filtering. "
                    "Check adata.obs and convert cluster/group labels to categorical columns."
                )
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

            existing_groupby = [name for name in self.groupby_editor.get_items() if name in obs_col_set]
            if not existing_groupby:
                existing_groupby = [c for c in obs_cols if c in {"sample_id", "sample", "condition", "batch", "donor"}]
                if not existing_groupby and obs_cols:
                    existing_groupby = [obs_cols[0]]
            self.groupby_editor.set_items(existing_groupby)
            self.additional_colors_editor.set_items(self._merge_unique(existing_additional, existing_groupby))

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

            if all_obs_cols:
                filtered_out = len(all_obs_cols) - len(obs_cols)
                if filtered_out > 0:
                    self._log(
                        f"Loaded {len(obs_cols)} qualitative obs columns (filtered out {filtered_out} high-cardinality/non-categorical fields)."
                    )
                else:
                    self._log(f"Loaded {len(obs_cols)} qualitative obs columns into additional_colors/groupby pickers.")
            if var_names:
                self._log(f"Loaded {len(var_names)} genes into manual genes picker.")
            largest_groups = sorted(unique_counts.items(), key=lambda item: item[1], reverse=True)[:3]
            if largest_groups:
                preview = ", ".join(f"{name}={count}" for name, count in largest_groups)
                self._log(f"Largest accepted group counts: {preview}")

            has_spatial = "spatial" in adata.obsm
            has_centroid = {"centroid_x", "centroid_y"}.issubset(set(all_obs_cols))
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
                f"obs columns: {len(obs_cols)} qualitative / {len(all_obs_cols)} total | cells: {adata.n_obs} | genes: {adata.n_vars} | "
                f"coords: {'obsm:spatial' if has_spatial else 'obs centroids' if has_centroid else 'not detected'}"
            )
            self.status_var.set("Inspection complete")
            self._update_action_states()
        except Exception as exc:
            friendly = self._format_anndata_read_error(exc)
            messagebox.showerror("Inspect failed", friendly)
            self._log(f"Inspect failed: {friendly}")
            self.status_var.set("Inspection failed")
            self._update_action_states()
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

    @staticmethod
    def _format_anndata_read_error(exc: Exception) -> str:
        message = str(exc)
        if "encoding_type='null'" in message and "No read method registered for IOSpec" in message:
            return (
                "This .h5ad was written with a newer AnnData encoding than this environment supports.\n\n"
                "If you are on Python 3.10, create a Python 3.11+ environment for this dataset and install:\n"
                "python -m pip install --upgrade \"anndata>=0.12\" \"h5py>=3.10\"\n\n"
                "Then restart KaroSpaceBuilder from that environment and inspect again."
            )
        return message

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
        if not self._matches_inspected_h5ad(h5ad_path):
            raise ValueError("Inspect Dataset first for the selected input .h5ad before building the viewer.")

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

        additional_colors = self._merge_unique(additional_colors, groupby_lists)

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
        try:
            from pandas.api.types import (
                is_bool_dtype,
                is_categorical_dtype,
                is_integer_dtype,
                is_object_dtype,
                is_string_dtype,
            )
        except Exception:
            is_bool_dtype = None
            is_categorical_dtype = None
            is_integer_dtype = None
            is_object_dtype = None
            is_string_dtype = None

        def unique_count(column: str) -> int | None:
            if column in unique_cache:
                return unique_cache[column]
            try:
                count = int(obs[column].nunique(dropna=True))
            except Exception:
                count = None
            unique_cache[column] = count
            return count

        def is_qualitative(column: str, count: int | None) -> bool:
            if is_categorical_dtype is None:
                return True
            try:
                series = obs[column]
                dtype = getattr(series, "dtype", None)
            except Exception:
                return True
            if is_categorical_dtype(dtype) or is_string_dtype(dtype) or is_object_dtype(dtype) or is_bool_dtype(dtype):
                return True
            if is_integer_dtype(dtype) and (count is None or count <= 24):
                return True
            return False

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
            if not is_qualitative(column, count):
                warnings.append(f"Skipping marker groupby '{column}': non-qualitative dtype.")
                continue
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
            if not is_qualitative(column, count):
                warnings.append(f"Skipping interaction groupby '{column}': non-qualitative dtype.")
                continue
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
            if not is_qualitative(column, count):
                warnings.append(f"Skipping neighbor groupby '{column}': non-qualitative dtype.")
                continue
            if count is not None and count < 2:
                warnings.append(f"Skipping neighbor groupby '{column}': only {count} unique value.")
                continue
            if count is not None and count > self._QUALITATIVE_MAX_UNIQUE:
                warnings.append(
                    f"Skipping neighbor groupby '{column}': {count:,} unique values exceeds "
                    f"{self._QUALITATIVE_MAX_UNIQUE:,}."
                )
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
        self._update_action_states()

    @staticmethod
    def _coerce_progress_value(value: object) -> int:
        try:
            percent = int(round(float(value)))
        except (TypeError, ValueError):
            percent = 0
        return max(0, min(100, percent))

    def _set_progress(self, value: object, stage: str | None = None) -> None:
        percent = self._coerce_progress_value(value)
        progress_color = self._app_palette["hover_bg"] if percent <= 0 else self._app_palette["accent"]
        self.progress.configure(progress_color=progress_color)
        self.progress.set(percent / 100.0)
        if stage:
            self.status_var.set(f"Export running... {percent}% | {stage}")

    def _on_export(self) -> None:
        if self._export_thread and self._export_thread.is_alive():
            messagebox.showinfo("Export running", "An export is already running.")
            return
        if not self._is_current_input_inspected():
            messagebox.showinfo("Inspection required", "Inspect Dataset first for the selected input .h5ad.")
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
            load_spatial_data, export_to_html = self._import_karospace_api(enable_numba_jit=config.enable_numba_jit)

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
