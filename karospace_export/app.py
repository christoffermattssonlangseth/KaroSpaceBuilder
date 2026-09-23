from __future__ import annotations

import contextlib
from dataclasses import dataclass
from datetime import datetime
from functools import partial
from pathlib import Path
import colorsys
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
    "background": "#ffffff",
    "text": "#1a1d23",
    "panel_bg": "#ffffff",
    "border": "#d8dbe1",
    "input_bg": "#f8f9fb",
    "muted": "#6b7280",
    "hover_bg": "#e6e8ed",
    "accent": _KI_COLORS["plum"],
    "accent_strong": _KI_COLORS["plum_dark"],
    "secondary": _KI_COLORS["plum"],
    "on_secondary": "#ffffff",
    "on_secondary_idle": _KI_COLORS["plum"],
    "danger": "#c9252d",
    "danger_hover": "#a51f25",
    "on_danger": "#ffffff",
    "on_accent": "#ffffff",
    "hero_bg": "#ffffff",
}

_KAROSPACE_DARK_PALETTE = {
    "plum_dark": _KI_COLORS["plum_dark"],
    "orange": _KI_COLORS["orange"],
    "light_orange": _KI_COLORS["light_orange"],
    "light_blue": _KI_COLORS["light_blue"],
    "plum": _KI_COLORS["plum"],
    "background": "#000000",
    "text": "#e8eaef",
    "panel_bg": "#1c1d22",
    "border": "#2e3038",
    "input_bg": "#24262c",
    "muted": "#9ca3af",
    "hover_bg": "#282a31",
    "accent": _KI_COLORS["orange"],
    "accent_strong": _KI_COLORS["plum"],
    "secondary": _KI_COLORS["orange"],
    "on_secondary": "#1a1d23",
    "on_secondary_idle": "#ffffff",
    "danger": "#e04b53",
    "danger_hover": "#b9333a",
    "on_danger": "#ffffff",
    "on_accent": "#1a1a1a",
    "hero_bg": "#000000",
}


def _palette_for_mode(mode: str) -> dict[str, str]:
    return (
        dict(_KAROSPACE_DARK_PALETTE)
        if str(mode).strip().lower() == "dark"
        else dict(_KAROSPACE_LIGHT_PALETTE)
    )


def _shift_hex_luminance(hex_color: str, points: float) -> str:
    return _transform_hex_luminance(hex_color, points=points)


def _set_hex_luminance(hex_color: str, luminance_percent: float) -> str:
    return _transform_hex_luminance(hex_color, absolute=luminance_percent)


def _is_dark_palette(palette: dict[str, str]) -> bool:
    return str(palette.get("background", "")).strip().lower() == "#000000"


def _secondary_idle_color(palette: dict[str, str]) -> str:
    return _set_hex_luminance(palette["secondary"], 5 if _is_dark_palette(palette) else 99)


def _transform_hex_luminance(
    hex_color: str,
    *,
    points: float | None = None,
    absolute: float | None = None,
) -> str:
    cleaned = hex_color.strip().lstrip("#")
    if len(cleaned) != 6:
        return hex_color
    try:
        red = int(cleaned[0:2], 16) / 255
        green = int(cleaned[2:4], 16) / 255
        blue = int(cleaned[4:6], 16) / 255
    except ValueError:
        return hex_color

    hue, luminance, saturation = colorsys.rgb_to_hls(red, green, blue)
    if absolute is not None:
        luminance = absolute / 100.0
    elif points is not None:
        luminance = luminance + (points / 100.0)
    luminance = max(0.0, min(1.0, luminance))
    red, green, blue = colorsys.hls_to_rgb(hue, luminance, saturation)
    return f"#{round(red * 255):02x}{round(green * 255):02x}{round(blue * 255):02x}"


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
    return {
        "root": {"fg_color": palette["background"]},
        "root_frame": {"fg_color": palette["background"], "corner_radius": 0},
        "card_frame": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 16,
            "border_width": 1,
            "border_color": palette["border"],
        },
        "hero_card": {
            "fg_color": palette.get("hero_bg", palette["panel_bg"]),
            "corner_radius": 20,
            "border_width": 2,
            "border_color": palette["secondary"],
        },
        "highlight_card": {
            "fg_color": palette["panel_bg"],
            "corner_radius": 14,
            "border_width": 1,
            "border_color": palette["secondary"],
        },
        "sub_frame": {"fg_color": "transparent", "corner_radius": 0},
        "section_label": {
            "font": (ui, 11, "bold"),
            "text_color": palette["secondary"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "hero_label": {
            "font": (ui, 30, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "header_label": {
            "font": (ui, 18, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "subheader_label": {
            "font": (ui, 12),
            "text_color": palette["muted"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "field_label": {
            "font": (ui, 12, "bold"),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "body_label": {
            "font": (ui, 11),
            "text_color": palette["text"],
            "fg_color": "transparent",
            "anchor": "w",
        },
        "primary_button": {
            "fg_color": palette["accent"],
            "hover_color": palette["accent_strong"],
            "text_color": palette.get("on_accent", "#ffffff"),
            "corner_radius": 12,
            "font": (ui, 12, "bold"),
            "height": 42,
            "border_width": 0,
        },
        "secondary_button": {
            "fg_color": _secondary_idle_color(palette),
            "hover_color": palette["secondary"],
            "text_color": palette.get("on_secondary_idle", palette["on_secondary"]),
            "corner_radius": 10,
            "font": (ui, 11),
            "height": 40,
            "border_width": 1,
            "border_color": palette["secondary"],
        },
        "pill_label": {
            "font": (ui, 10, "bold"),
            "text_color": palette.get("on_accent", "#ffffff"),
            "fg_color": palette["accent"],
            "corner_radius": 999,
            "padx": 12,
            "pady": 6,
        },
        "muted_pill_label": {
            "font": (ui, 10, "bold"),
            "text_color": palette["text"],
            "fg_color": palette["hover_bg"],
            "corner_radius": 999,
            "padx": 12,
            "pady": 6,
        },
        "entry": {
            "fg_color": palette["input_bg"],
            "text_color": palette["text"],
            "placeholder_text_color": palette["muted"],
            "border_color": palette["border"],
            "border_width": 1,
            "corner_radius": 10,
            "height": 40,
            "font": (ui, 12),
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
            "font": (ui, 12),
        },
        "checkbox": {
            "fg_color": palette["accent"],
            "hover_color": palette["accent_strong"],
            "checkmark_color": palette.get("on_accent", "#ffffff"),
            "text_color": palette["text"],
            "border_color": palette["border"],
            "font": (ui, 11),
            "corner_radius": 6,
        },
        "tabview": {
            "fg_color": palette["panel_bg"],
            "segmented_button_fg_color": palette["hover_bg"],
            "segmented_button_selected_color": palette["accent"],
            "segmented_button_selected_hover_color": palette["accent"],
            "segmented_button_unselected_color": palette["hover_bg"],
            "segmented_button_unselected_hover_color": _set_hex_luminance(palette["secondary"], 98),
            "text_color": palette["text"],
            "corner_radius": 0,
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
            "border_color": palette["secondary"],
            "corner_radius": 10,
            "border_width": 1,
            "font": (mono, 10),
            "scrollbar_button_color": palette["hover_bg"],
            "scrollbar_button_hover_color": palette["accent"],
        },
        "progress": {
            "fg_color": palette["hover_bg"],
            "progress_color": palette["secondary"],
            "border_color": palette["border"],
            "corner_radius": 10,
            "height": 18,
        },
    }


def _get_anndata():
    import anndata as ad

    return ad


def _get_numpy():
    import numpy as np

    return np


def _style_secondary_button_for_palette(button: tk.Widget, palette: dict[str, str], *, borderless: bool = False) -> None:
    try:
        button.configure(
            fg_color=_secondary_idle_color(palette),
            hover_color=palette["secondary"],
            text_color=palette.get("on_secondary_idle", palette["on_secondary"]),
            border_width=0 if borderless else 1,
            border_color=palette["secondary"],
        )
    except Exception:
        pass


def _bind_secondary_button_feedback(button: tk.Widget, palette_getter, *, reset_callback=None) -> None:
    def pointer_inside() -> bool:
        try:
            pointer_x, pointer_y = button.winfo_pointerxy()
            root_x = button.winfo_rootx()
            root_y = button.winfo_rooty()
            return root_x <= pointer_x <= root_x + button.winfo_width() and root_y <= pointer_y <= root_y + button.winfo_height()
        except Exception:
            return False

    def reset() -> None:
        if reset_callback is not None:
            reset_callback()
        else:
            _style_secondary_button_for_palette(button, palette_getter())

    def on_enter(_event: object) -> None:
        try:
            if button.cget("state") == "disabled":
                return
            palette = palette_getter()
            button.configure(fg_color=palette["secondary"], text_color=palette["on_secondary"])
        except Exception:
            pass

    def on_press(_event: object) -> None:
        try:
            if button.cget("state") == "disabled":
                return
            palette = palette_getter()
            button.configure(
                fg_color=_shift_hex_luminance(palette["secondary"], -10),
                text_color=palette["on_secondary"],
            )
        except Exception:
            pass

    def on_release(_event: object) -> None:
        try:
            if button.cget("state") == "disabled":
                return
            if pointer_inside():
                palette = palette_getter()
                button.configure(fg_color=palette["secondary"], text_color=palette["on_secondary"])
            else:
                reset()
        except Exception:
            pass

    def on_leave(_event: object) -> None:
        try:
            if button.cget("state") == "disabled":
                return
            reset()
        except Exception:
            pass

    button.bind("<Enter>", on_enter, add="+")
    button.bind("<ButtonPress-1>", on_press, add="+")
    button.bind("<ButtonRelease-1>", on_release, add="+")
    button.bind("<Leave>", on_leave, add="+")


class SearchableListEditor(_CTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        *,
        label: str,
        height: int = 8,
        help_text: str | None = None,
        palette: dict[str, str] | None = None,
        on_change=None,
        stack_controls: bool = False,
        allow_select_all: bool = False,
    ) -> None:
        self._palette = dict(palette or _palette_for_mode("dark"))
        self._theme = _ctk_theme_config(self._palette)
        self._on_change = on_change
        self._suspend_change_notification = False
        self._allow_select_all = allow_select_all
        super().__init__(parent, **self._theme["sub_frame"])
        self._choices: list[str] = []
        self._required_items: list[str] = []
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
        if stack_controls:
            controls.columnconfigure(1, weight=1)
            controls.columnconfigure(2, weight=1)
            self.entry.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        else:
            self.entry.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        self.entry.bind("<KeyRelease>", self._on_search)
        self.entry.bind("<Return>", lambda _event: self.add_current())

        self.add_btn = ctk.CTkButton(controls, text="+ Add", command=self.add_current, width=88, **self._theme["secondary_button"])
        _bind_secondary_button_feedback(self.add_btn, lambda: self._palette)
        if stack_controls:
            self.add_btn.grid(row=1, column=0, sticky="ew", padx=(0, 6))
        else:
            self.add_btn.grid(row=0, column=1, padx=(0, 6))
        self.remove_btn = ctk.CTkButton(
            controls,
            text="Remove",
            command=self.remove_selected,
            width=88,
            **self._theme["secondary_button"],
        )
        _bind_secondary_button_feedback(self.remove_btn, lambda: self._palette)
        if stack_controls:
            self.remove_btn.grid(row=1, column=1, sticky="ew", padx=(0, 6))
        else:
            self.remove_btn.grid(row=0, column=2, padx=(0, 6))
        self.clear_btn = ctk.CTkButton(controls, text="Clear", command=self.clear, width=88, **self._theme["secondary_button"])
        _bind_secondary_button_feedback(self.clear_btn, lambda: self._palette)
        if stack_controls:
            self.clear_btn.grid(row=1, column=2, sticky="ew")
        else:
            self.clear_btn.grid(row=0, column=3)

        self.select_all_btn: ctk.CTkButton | None = None
        if self._allow_select_all:
            self.select_all_btn = ctk.CTkButton(
                controls,
                text="Select all",
                command=self.add_all,
                width=96,
                **self._theme["secondary_button"],
            )
            _bind_secondary_button_feedback(self.select_all_btn, lambda: self._palette)
            if stack_controls:
                controls.columnconfigure(3, weight=1)
                self.select_all_btn.grid(row=1, column=3, sticky="ew", padx=(6, 0))
            else:
                self.select_all_btn.grid(row=0, column=4, padx=(6, 0))

        list_wrap = ctk.CTkFrame(self, **self._theme["sub_frame"])
        list_wrap.grid(row=2, column=0, sticky="ew", pady=(8, 0))
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

        self.help_label: ctk.CTkLabel | None = None
        if help_text:
            self.help_label = ctk.CTkLabel(self, text=help_text, **self._theme["subheader_label"])
            self.help_label.configure(wraplength=520, justify="left")
            self.help_label.grid(row=3, column=0, sticky="ew", pady=(6, 0))

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
        self._choices = self._dedupe_items(values)
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
        self._notify_change()

    def remove_selected(self) -> None:
        required = set(self._required_items)
        removed = False
        for idx in reversed(self.listbox.curselection()):
            if str(self.listbox.get(idx)) in required:
                continue
            self.listbox.delete(idx)
            removed = True
        if removed:
            self._notify_change()

    def clear(self) -> None:
        required = set(self._required_items)
        removed = False
        for idx in range(self.listbox.size() - 1, -1, -1):
            if str(self.listbox.get(idx)) in required:
                continue
            self.listbox.delete(idx)
            removed = True
        if removed:
            self._notify_change()

    def add_all(self) -> None:
        if not self._choices:
            return
        self.set_items([*self.get_items(), *self._choices])

    def set_items(self, values: list[str] | tuple[str, ...]) -> None:
        self._suspend_change_notification = True
        self.listbox.delete(0, "end")
        for value in self._with_required_items(values):
            self.listbox.insert("end", value)
        self._suspend_change_notification = False
        self._notify_change()

    def get_items(self) -> list[str]:
        return [str(v) for v in self.listbox.get(0, "end")]

    def set_required_items(self, values: list[str] | tuple[str, ...]) -> None:
        previous_required = set(self._required_items)
        required = self._dedupe_items(values)
        required_set = set(required)
        existing = [item for item in self.get_items() if item not in previous_required or item in required_set]
        self._required_items = required
        self.set_items(existing)

    @staticmethod
    def _dedupe_items(values: list[str] | tuple[str, ...]) -> list[str]:
        seen: set[str] = set()
        ordered: list[str] = []
        for raw in values:
            value = str(raw).strip()
            if not value or value in seen:
                continue
            seen.add(value)
            ordered.append(value)
        return ordered

    def _with_required_items(self, values: list[str] | tuple[str, ...]) -> list[str]:
        return self._dedupe_items([*self._required_items, *values])

    def _notify_change(self) -> None:
        if self._suspend_change_notification:
            return
        if self._on_change is None:
            return
        try:
            self._on_change()
        except Exception:
            pass

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


class SectionOrderEditor(_CTK_FRAME_BASE):
    def __init__(
        self,
        parent,
        *,
        variable: tk.StringVar,
        palette: dict[str, str] | None = None,
        on_change=None,
    ) -> None:
        self._palette = dict(palette or _palette_for_mode("dark"))
        self._theme = _ctk_theme_config(self._palette)
        self._variable = variable
        self._on_change = on_change
        self._drag_index: int | None = None
        self._syncing = False
        self._enabled = True
        super().__init__(parent, **self._theme["sub_frame"])

        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        list_wrap = ctk.CTkFrame(self, **self._theme["sub_frame"])
        list_wrap.grid(row=0, column=0, sticky="ew")
        list_wrap.columnconfigure(0, weight=1)

        self.listbox = tk.Listbox(
            list_wrap,
            height=6,
            selectmode="browse",
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
        self.listbox.bind("<Button-1>", self._start_drag)
        self.listbox.bind("<B1-Motion>", self._drag_to)
        self.listbox.bind("<ButtonRelease-1>", self._end_drag)

        controls = ctk.CTkFrame(self, **self._theme["sub_frame"])
        controls.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        self.up_btn = ctk.CTkButton(controls, text="Move up", command=lambda: self._move_selected(-1), width=96, **self._theme["secondary_button"])
        self.down_btn = ctk.CTkButton(controls, text="Move down", command=lambda: self._move_selected(1), width=112, **self._theme["secondary_button"])
        self.reset_btn = ctk.CTkButton(controls, text="Reset order", command=self.reset_to_detected_order, width=112, **self._theme["secondary_button"])
        for button in (self.up_btn, self.down_btn, self.reset_btn):
            _bind_secondary_button_feedback(button, lambda: self._palette)
        self.up_btn.pack(side="left")
        self.down_btn.pack(side="left", padx=(8, 0))
        self.reset_btn.pack(side="left", padx=(8, 0))

        self.info_label = ctk.CTkLabel(
            self,
            text="Drag section values to reorder. Export uses this order.",
            **self._theme["subheader_label"],
        )
        self.info_label.configure(wraplength=520, justify="left")
        self.info_label.grid(row=2, column=0, sticky="ew", pady=(6, 0))

        self._detected_values: list[str] = []
        self._variable.trace_add("write", lambda *_: self._sync_from_variable())
        self.apply_palette(self._palette)

    def set_values(self, values: list[str] | tuple[str, ...], *, preserve_current: bool = True) -> None:
        detected = self._dedupe(values)
        current = self.get_items() if preserve_current else []
        current = [value for value in current if value in set(detected)]
        merged = [*current, *[value for value in detected if value not in set(current)]]
        self._detected_values = detected
        self._set_items(merged)

    def clear(self) -> None:
        self._detected_values = []
        self._set_items([])

    def reset_to_detected_order(self) -> None:
        self._set_items(self._detected_values)

    def get_items(self) -> list[str]:
        return [str(value) for value in self.listbox.get(0, "end")]

    def set_enabled(self, enabled: bool) -> None:
        self._enabled = bool(enabled)
        state = "normal" if enabled else "disabled"
        self.listbox.configure(state=state)
        self.up_btn.configure(state=state)
        self.down_btn.configure(state=state)
        self.reset_btn.configure(state=state)

    def apply_palette(self, palette: dict[str, str]) -> None:
        self._palette = dict(palette)
        self._theme = _ctk_theme_config(self._palette)
        self.up_btn.configure(**self._theme["secondary_button"])
        self.down_btn.configure(**self._theme["secondary_button"])
        self.reset_btn.configure(**self._theme["secondary_button"])
        self.info_label.configure(**self._theme["subheader_label"])
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

    @staticmethod
    def _dedupe(values: list[str] | tuple[str, ...]) -> list[str]:
        seen: set[str] = set()
        out: list[str] = []
        for raw in values:
            value = str(raw).strip()
            if not value or value in seen:
                continue
            seen.add(value)
            out.append(value)
        return out

    def _parse_variable(self) -> list[str]:
        return self._dedupe(tuple(item.strip() for item in self._variable.get().split(",")))

    def _set_items(self, values: list[str]) -> None:
        self._syncing = True
        restore_state = str(self.listbox.cget("state"))
        if restore_state == "disabled":
            self.listbox.configure(state="normal")
        self.listbox.delete(0, "end")
        for value in values:
            self.listbox.insert("end", value)
        self._write_variable()
        self.listbox.configure(state=restore_state)
        self._syncing = False
        self._notify_change()

    def _write_variable(self) -> None:
        self._variable.set(",".join(self.get_items()))

    def _sync_from_variable(self) -> None:
        if self._syncing:
            return
        values = self._parse_variable()
        if values == self.get_items():
            return
        self._set_items(values)

    def _notify_change(self) -> None:
        if self._on_change is None:
            return
        try:
            self._on_change()
        except Exception:
            pass

    def _start_drag(self, event) -> None:
        index = self.listbox.nearest(event.y)
        self._drag_index = index if 0 <= index < self.listbox.size() else None

    def _drag_to(self, event) -> None:
        if self._drag_index is None:
            return
        target = self.listbox.nearest(event.y)
        if target == self._drag_index or target < 0 or target >= self.listbox.size():
            return
        value = self.listbox.get(self._drag_index)
        self.listbox.delete(self._drag_index)
        self.listbox.insert(target, value)
        self.listbox.selection_clear(0, "end")
        self.listbox.selection_set(target)
        self._drag_index = target
        self._write_variable()

    def _end_drag(self, _event) -> None:
        if self._drag_index is None:
            return
        self._drag_index = None
        self._write_variable()
        self._notify_change()

    def _move_selected(self, delta: int) -> None:
        selection = self.listbox.curselection()
        if not selection:
            return
        index = int(selection[0])
        target = index + int(delta)
        if target < 0 or target >= self.listbox.size():
            return
        value = self.listbox.get(index)
        self.listbox.delete(index)
        self.listbox.insert(target, value)
        self.listbox.selection_set(target)
        self.listbox.see(target)
        self._write_variable()
        self._notify_change()


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
    output_html_path: Path
    coords_mode: str | None
    spatial_key: str
    spatial_columns: tuple[str, str] | None
    spatialdata_table: str | None
    section_groupby: str
    section_order: list[str] | None
    section_metadata: list[str] | None
    section_metadata_extra: list[str] | None
    metadata_value_order: dict[str, list[str]] | None
    metadata_max_columns: int | None
    initial_color: str
    title: str
    enable_numba_jit: bool
    outline_by: str | None
    metadata_labels: dict[str, str] | None
    viewer_info_html: str | None
    tutorial: bool
    embed_reproducibility_info: bool
    source_input_path: str | None
    min_panel_size: int
    spot_size: float | str | None
    downsample: int | None
    additional_colors: list[str] | None
    genes: list[str] | None
    features_list: str | None
    feature_encoding: str
    feature_value_encoding: str
    feature_storage: str
    also_export_karospace: bool
    feature_manifest_path: str | None
    feature_sidecar_shard_size: int
    feature_sparse_zero_threshold: float
    modalities: list[str] | None
    neighbor_stats_groupby: list[str] | None
    neighbor_stats_permutations: int | None
    neighbor_stats_seed: int
    statistics_additional_annotations: list[str] | None
    statistics_modalities: list[str] | str | None
    statistics_contrast_categories: dict[str, list[str] | None] | list[str] | None
    statistics_counts_layer: str | None
    statistics_normalization: str
    statistics_scale_factor: float
    statistics_normalized_layer: str | None
    statistics_min_cell_counts: int
    statistics_min_feature_counts: int
    statistics_n_cpus: int
    wilcoxon: str
    wilcoxon_runtime_limit: str
    wilcoxon_min_cells_per_group: int
    wilcoxon_min_pct_expressed: float
    wilcoxon_p_adjust_method: str
    wilcoxon_padj_cutoff: float
    wilcoxon_log2fc_cutoff: float
    wilcoxon_embed_top_n_per_comparison: int
    wilcoxon_top_n_per_category: int
    pseudobulk: str | None
    pseudobulk_replicate_annotation: str | None
    pseudobulk_min_cells_per_pseudobulk: int
    pseudobulk_min_replicates: int
    pseudobulk_min_pct_expressed: float
    pseudobulk_p_adjust_method: str
    pseudobulk_padj_cutoff: float
    pseudobulk_log2fc_cutoff: float
    pseudobulk_deseq2_fit_type: str
    pseudobulk_embed_top_n_per_comparison: int
    pathway_gmt: list[str] | None
    pathway: str | None
    pathway_organism: str
    pathway_top_n: int
    pathway_min_overlap: int
    pathway_gsea_permutations: int
    interaction_markers: str | None
    interaction_markers_top_targets: int
    interaction_markers_top_features: int
    interaction_markers_min_cells: int
    interaction_markers_min_neighbors: int
    section_rotations: dict[str, float] | None
    deconvolutions: dict[str, str] | None
    feature_correlation_top_n: int
    spatial_variable_features_n: int
    scalebar_unit: str
    section_images: dict[str, object] | None
    section_images_max_px: int


class _ThreadingHTTPServer(socketserver.ThreadingMixIn, http.server.HTTPServer):
    daemon_threads = True


class _ExportCancelled(Exception):
    pass


class _EventLogTee:
    def __init__(self, stream, event_queue: queue.Queue[tuple[str, object]], *, prefix: str = "") -> None:
        self._stream = stream
        self._queue = event_queue
        self._prefix = prefix
        self._buffer = ""

    @property
    def encoding(self):
        return getattr(self._stream, "encoding", None)

    def isatty(self) -> bool:
        try:
            return bool(self._stream.isatty())
        except Exception:
            return False

    def write(self, text) -> int:
        if text is None:
            return 0
        value = str(text)
        try:
            self._stream.write(value)
        except Exception:
            pass

        self._buffer += value.replace("\r", "\n")
        while "\n" in self._buffer:
            line, self._buffer = self._buffer.split("\n", 1)
            self._emit(line)
        return len(value)

    def flush(self) -> None:
        try:
            self._stream.flush()
        except Exception:
            pass
        if self._buffer.strip():
            self._emit(self._buffer)
        self._buffer = ""

    def _emit(self, line: str) -> None:
        message = line.strip()
        if not message:
            return
        self._queue.put(("log", f"{self._prefix}{message}"))


class ExportApp(ctk.CTk if ctk is not None else object):
    _OUTPUT_HTML_BASENAME = "KaroSpace"
    _HTML_SIZE_MB_PER_ELEMENT = 2.313e-6
    _HTML_SIZE_WARNING_ELEMENTS = 150_000_000
    _PARAMETER_STATE_SCHEMA = "karospace_builder_parameters"
    _PARAMETER_LIST_EDITORS = (
        "modalities_editor",
        "additional_colors_editor",
        "statistics_annotations_editor",
        "statistics_modalities_editor",
        "section_metadata_editor",
        "section_metadata_extra_editor",
        "manual_genes_editor",
    )

    def __init__(self) -> None:
        super().__init__()
        self.title("KaroSpaceBuilder")
        self.geometry("1280x820")
        self.minsize(1060, 720)

        self._queue: queue.Queue[tuple[str, object]] = queue.Queue()
        self._export_thread: threading.Thread | None = None
        self._cancel_requested = threading.Event()
        self._server: _ThreadingHTTPServer | None = None
        self._server_thread: threading.Thread | None = None
        self._last_outdir: Path | None = None
        self._last_output_html: Path | None = None
        self._has_valid_input_file = False
        self._has_inspected_input_file = False
        self._inspect_locked_after_press = False
        self._inspected_h5ad_path: Path | None = None
        self._inspected_coords_mode: str | None = None
        self._inspected_n_cells: int | None = None
        self._inspected_n_genes: int | None = None
        self._inspected_obs_cols: set[str] = set()
        self._inspected_section_key_cols: list[str] = []
        self._inspected_feature_names_by_modality: dict[str, list[str]] = {}
        self._inspected_obsm_keys: list[str] = []
        self._inspected_spatial_obs_cols: list[str] = []
        self._inspected_layer_keys: list[str] = []
        self._inspected_section_counts_by_column: dict[str, list[int]] = {}
        self._inspected_section_values_by_column: dict[str, list[str]] = {}
        self._feature_items_by_modality: dict[str, list[str]] = {}
        self._active_feature_modality = "rna"
        self._spatialdata_table_choices: list[str] = []
        self._spatialdata_table_widgets: list[tk.Widget] = []
        self._spatialdata_table_source_path: Path | None = None
        self._syncing_spatialdata_table = False
        self._loading_parameter_state = False
        self._runtime_chip_animation_after_id: str | None = None
        self._runtime_chip_animation_step = 0
        self._themed_widgets: dict[str, list[tk.Widget]] = {}
        self._inspection_gated_widgets: list[tk.Widget] = []
        self._right_panel_widgets: list[tk.Widget] = []
        self._secondary_action_buttons: list[tk.Widget] = []
        self._tab_hover_bound: set[str] = set()
        self._placeholder_refreshers: list[object] = []
        self._hidden_widget_layouts: dict[tk.Widget, tuple[str, dict[str, object]]] = {}

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
        if hasattr(self, "theme_toggle_btn"):
            self.theme_toggle_btn.configure(**self._theme["secondary_button"])
            self._sync_theme_toggle()
        if hasattr(self, "main_scroll_frame"):
            self.main_scroll_frame.configure(
                scrollbar_button_color=palette["hover_bg"],
                scrollbar_button_hover_color=palette["accent"],
            )
        self._sync_runtime_chip()
        for attr in (
            "additional_colors_editor",
            "statistics_annotations_editor",
            "statistics_modalities_editor",
            "section_metadata_editor",
            "section_metadata_extra_editor",
            "manual_genes_editor",
            "modalities_editor",
            "section_order_editor",
        ):
            widget = getattr(self, attr, None)
            if widget is not None and hasattr(widget, "apply_palette"):
                widget.apply_palette(palette)
        for button in self._secondary_action_buttons:
            self._style_secondary_action_button(button)
        for refresh_placeholder in self._placeholder_refreshers:
            try:
                refresh_placeholder()
            except Exception:
                pass
        self._sync_tab_button_styles()
        if hasattr(self, "inspect_btn"):
            self._refresh_input_gate()
        if hasattr(self, "export_estimate_summary_label"):
            self._update_export_estimate()

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
        self._bind_secondary_press_feedback(button)
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

    def _toggle_theme(self) -> None:
        current = self.theme_var.get().strip().lower()
        next_mode = "light" if current == "dark" else "dark"
        self.theme_var.set(next_mode)
        self._apply_app_theme(next_mode)

    def _sync_theme_toggle(self) -> None:
        if not hasattr(self, "theme_toggle_btn"):
            return
        mode = self.theme_var.get().strip().lower()
        self.theme_toggle_btn.configure(text="☀" if mode == "dark" else "☾")

    @staticmethod
    def _guide_text() -> str:
        return (
            "General tab\n"
            "- Input file: AnnData .h5ad file or SpatialData .zarr directory.\n"
            "- Parameter JSON saves the current form state or imports a saved form state.\n"
            "- SpatialData table: shown only for Zarr inputs; table choices are discovered automatically and the first table is selected by default.\n"
            "- Output path: maps to CLI --output. Relative paths are resolved inside Output directory; .karospace output requires sidecar feature storage.\n"
            "- Spatial coordinates: maps to --spatial-key/--spatial-x/--spatial-y.\n"
            "- Modality selects which detected modalities are embedded/exported and limits the Features tab modality dropdown.\n"
            "- Section key, section order, and SpatialData table map to load_spatial_data(...).\n"
            "- Runtime mode: optional numba JIT performance mode (can be less stable in frozen app).\n"
            "- Title maps to export_to_html.\n\n"
            "Metadata tab\n"
            "- cell_annotations maps to additional cell annotation columns in the new KaroSpace API.\n"
            "- section_metadata and section_metadata_extra control section filters and stored section metadata.\n"
            "- Metadata value order, labels, and max columns map to the matching API/CLI metadata arguments.\n\n"
            "Features tab\n"
            "- features: build a var_names list manually with the searchable picker.\n"
            "- features_list points to a one-feature-per-line file and can be used alone or combined with features.\n"
            "- feature_storage, feature_encoding, feature value encoding, sidecar shard size, sparse threshold, and modalities map directly to KaroSpace feature options.\n\n"
            "Statistics tab\n"
            "- Statistics annotations, modalities, contrast categories, count layers, normalization, filters, and CPU settings map to the current statistics_* API.\n"
            "- Enable Wilcoxon, pseudobulk, or pathway analysis with checkboxes to show each group's arguments.\n\n"
            "Neighborhoods tab\n"
            "- Neighbor statistics and interaction marker controls mirror the API/CLI neighborhood arguments.\n\n"
            "Connections tab\n"
            "- Discovery panel limits control optional feature correlation and spatial variable feature summaries.\n\n"
            "Viewer tab\n"
            "- Downsample: integer cells per section (blank keeps all).\n"
            "- Tutorial, reproducibility info, viewer info HTML, and viewer info HTML file map to viewer rendering arguments.\n\n"
            "Sections tab\n"
            "- Min panel size, spot size, scalebar unit, section rotations, section images, image size, and preview server settings live here.\n"
            "- Deconvolutions JSON maps to viewer deconvolution rendering arguments.\n"
            "- Serve after export starts a local preview server for the configured output.\n\n"
            "Tip: click Inspect Dataset to load searchable obs dropdown choices, SpatialData table columns, and feature pickers."
        )

    def _show_guide_popup(self) -> None:
        existing = getattr(self, "_guide_popup", None)
        if existing is not None and existing.winfo_exists():
            existing.destroy()
            return

        popup = ctk.CTkToplevel(self)
        self._guide_popup = popup
        popup.overrideredirect(True)
        popup.geometry("720x620")
        popup.minsize(560, 440)
        popup.configure(**self._theme["root"])
        popup.columnconfigure(0, weight=1)
        popup.rowconfigure(1, weight=1)
        popup.bind("<Escape>", lambda _event: popup.destroy())

        self.update_idletasks()
        x = self.winfo_rootx() + max(24, int((self.winfo_width() - 720) / 2))
        y = self.winfo_rooty() + 96
        popup.geometry(f"720x620+{x}+{y}")

        header = ctk.CTkFrame(popup, **self._theme["sub_frame"])
        header.grid(row=0, column=0, sticky="ew", padx=18, pady=(18, 10))
        header.columnconfigure(0, weight=1)
        ctk.CTkLabel(header, text="Guide & Workflow", **self._theme["header_label"]).grid(row=0, column=0, sticky="w")
        close_btn = ctk.CTkButton(header, text="×", command=popup.destroy, width=42, **self._theme["secondary_button"])
        close_btn.grid(row=0, column=1, sticky="e")

        text = ctk.CTkTextbox(popup, wrap="word", **self._theme["textbox"])
        text.grid(row=1, column=0, sticky="nsew", padx=18, pady=(0, 18))
        text.insert("1.0", self._guide_text())
        text.configure(state="disabled")
        popup.focus()

    def _on_tab_selected(self) -> None:
        if not hasattr(self, "notebook"):
            return
        selected = self.notebook.get()
        gated_tabs = {"Metadata", "Features", "Statistics", "Neighborhoods", "Connections", "Viewer", "Sections"}
        if selected in gated_tabs and not self._has_inspected_input_file:
            self.notebook.set(self._last_enabled_tab)
            self._sync_tab_button_styles()
            return
        self._last_enabled_tab = selected
        self._sync_tab_button_styles()

    def _valid_h5ad_path(self) -> Path | None:
        raw = self.h5ad_var.get().strip()
        if not raw:
            return None
        path = Path(raw).expanduser()
        try:
            resolved = path.resolve()
        except Exception:
            return None
        if resolved.exists() and (
            (resolved.is_file() and resolved.suffix.lower() == ".h5ad")
            or self._is_zarr_path(resolved)
        ):
            return resolved
        return None

    @staticmethod
    def _is_zarr_path(path: Path | None) -> bool:
        if path is None or not path.is_dir():
            return False
        return (
            path.suffix.lower() == ".zarr"
            or (path / "zarr.json").exists()
            or (path / ".zgroup").exists()
            or (path / "tables").is_dir()
        )

    def _sync_spatialdata_table_selector(self, valid_path: Path | None, *, busy: bool = False) -> None:
        is_zarr = self._is_zarr_path(valid_path)
        self._set_widgets_visible(self._spatialdata_table_widgets, is_zarr)
        if not hasattr(self, "spatialdata_table_combo"):
            return
        if not is_zarr:
            self._spatialdata_table_choices = []
            self._syncing_spatialdata_table = True
            try:
                self.spatialdata_table_var.set("")
                self.spatialdata_table_combo.configure(values=[], state="disabled")
            finally:
                self._syncing_spatialdata_table = False
            return

        values = list(self._spatialdata_table_choices)
        current = self.spatialdata_table_var.get().strip()
        if current and current not in values:
            values = [current, *values]
        state = "readonly" if values and not busy else "disabled"
        self.spatialdata_table_combo.configure(values=values, state=state)

    def _autoselect_spatialdata_table(self, valid_path: Path | None) -> None:
        if not self._is_zarr_path(valid_path):
            return
        assert valid_path is not None

        if self._spatialdata_table_source_path != valid_path:
            table_choices = self._discover_spatialdata_tables(valid_path)
            self._spatialdata_table_choices = table_choices
            self._spatialdata_table_source_path = valid_path
            current = self.spatialdata_table_var.get().strip()
            selected = current if current in table_choices else (table_choices[0] if table_choices else "")
            self._syncing_spatialdata_table = True
            try:
                self.spatialdata_table_var.set(selected)
            finally:
                self._syncing_spatialdata_table = False
            if table_choices:
                self._log(
                    f"Found {len(table_choices)} SpatialData table(s); selected '{selected}'."
                )
            else:
                self._log("No SpatialData tables were found in the selected Zarr directory.")
            return

        if self._spatialdata_table_choices and not self.spatialdata_table_var.get().strip():
            self._syncing_spatialdata_table = True
            try:
                self.spatialdata_table_var.set(self._spatialdata_table_choices[0])
            finally:
                self._syncing_spatialdata_table = False

    def _input_ready_for_inspect(self, valid_path: Path | None) -> bool:
        if valid_path is None:
            return False
        if self._is_zarr_path(valid_path):
            return bool(self.spatialdata_table_var.get().strip())
        return True

    def _on_spatialdata_table_selected(self, _choice: str | None = None) -> None:
        if self._syncing_spatialdata_table:
            return
        path = self._valid_h5ad_path()
        if not self._is_zarr_path(path):
            return
        if self._has_inspected_input_file:
            self._has_inspected_input_file = False
            self._clear_inspection_metadata()
            self._inspect_locked_after_press = False
            self._refresh_input_gate()
            self._log("SpatialData table changed. Inspect Dataset again to load fields from the selected table.")

    def _discover_spatialdata_tables(self, path: Path) -> list[str]:
        tables_dir = path / "tables"
        if tables_dir.is_dir():
            return sorted(child.name for child in tables_dir.iterdir() if child.is_dir())

        try:
            import importlib

            data_loader = importlib.import_module("karospace.data_loader")
            read_spatialdata_zarr = getattr(data_loader, "_read_spatialdata_zarr")
            sdata = read_spatialdata_zarr(str(path))
            tables = getattr(sdata, "tables", None)
            if tables is not None and hasattr(tables, "keys"):
                return [str(key) for key in tables.keys()]
        except Exception as exc:
            self._log(f"Could not pre-read SpatialData table choices: {exc}")
        return []

    def _configure_widget_state(self, widget: object | None, enabled: bool) -> None:
        if widget is None:
            return
        if hasattr(widget, "set_enabled"):
            try:
                widget.set_enabled(enabled)
            except Exception:
                pass
        try:
            widget.configure(state="normal" if enabled else "disabled")
        except Exception:
            pass
        if enabled:
            for role, widgets in self._themed_widgets.items():
                if widget in widgets:
                    style = self._theme.get(role)
                    if style:
                        try:
                            widget.configure(**style)
                        except Exception:
                            pass
                    break
            return

        disabled_styles = {
            "fg_color": self._app_palette["hover_bg"],
            "text_color": self._app_palette["muted"],
            "button_color": self._app_palette["hover_bg"],
            "button_hover_color": self._app_palette["hover_bg"],
            "hover_color": self._app_palette["hover_bg"],
            "border_color": self._app_palette["border"],
        }
        for key, value in disabled_styles.items():
            try:
                widget.configure(**{key: value})
            except Exception:
                pass

    def _set_widgets_visible(self, widgets: list[tk.Widget], visible: bool) -> None:
        for widget in widgets:
            try:
                if visible:
                    saved = self._hidden_widget_layouts.pop(widget, None)
                    if saved is None:
                        continue
                    manager, info = saved
                    if manager == "grid":
                        widget.grid(**info)
                    elif manager == "pack":
                        widget.pack(**info)
                    elif manager == "place":
                        widget.place(**info)
                else:
                    if widget in self._hidden_widget_layouts:
                        continue
                    manager = widget.winfo_manager()
                    if manager == "grid":
                        info = dict(widget.grid_info())
                        self._hidden_widget_layouts[widget] = (manager, info)
                        widget.grid_remove()
                    elif manager == "pack":
                        info = dict(widget.pack_info())
                        self._hidden_widget_layouts[widget] = (manager, info)
                        widget.pack_forget()
                    elif manager == "place":
                        info = dict(widget.place_info())
                        self._hidden_widget_layouts[widget] = (manager, info)
                        widget.place_forget()
            except Exception:
                pass

    def _refresh_tab_gate(self) -> None:
        if not hasattr(self, "notebook"):
            return
        segmented = getattr(self.notebook, "_segmented_button", None)
        buttons = getattr(segmented, "_buttons_dict", {}) if segmented is not None else {}
        for tab_name in ("Metadata", "Features", "Statistics", "Neighborhoods", "Connections", "Viewer", "Sections"):
            button = buttons.get(tab_name)
            self._configure_widget_state(button, self._has_inspected_input_file)
        if not self._has_inspected_input_file and self.notebook.get() in {"Metadata", "Features", "Statistics", "Neighborhoods", "Connections", "Viewer", "Sections"}:
            self.notebook.set("General")
            self._last_enabled_tab = "General"
        self._sync_tab_button_styles()

    def _sync_tab_button_styles(self) -> None:
        if not hasattr(self, "notebook"):
            return
        segmented = getattr(self.notebook, "_segmented_button", None)
        buttons = getattr(segmented, "_buttons_dict", {}) if segmented is not None else {}
        selected = self.notebook.get()
        if segmented is not None:
            try:
                segmented.grid_configure(sticky="ew", padx=0)
                for column_index in range(len(buttons)):
                    segmented.columnconfigure(column_index, weight=1)
                segmented.configure(
                    fg_color=self._app_palette["hover_bg"],
                    unselected_color=self._app_palette["hover_bg"],
                    unselected_hover_color=_set_hex_luminance(self._app_palette["secondary"], 98),
                    corner_radius=0,
                )
            except Exception:
                pass

        for tab_name, button in buttons.items():
            try:
                disabled = button.cget("state") == "disabled"
            except Exception:
                disabled = False

            try:
                button.configure(
                    border_width=0,
                    border_color=self._app_palette["hover_bg"],
                    corner_radius=10,
                    text_color=(
                        self._app_palette["muted"]
                        if disabled
                        else self._app_palette.get("on_accent", "#ffffff")
                        if tab_name == selected
                        else self._app_palette["text"]
                    ),
                )
            except Exception:
                pass

            if tab_name in self._tab_hover_bound:
                continue

            def on_enter(_event: object, name: str = tab_name, tab_button: tk.Widget = button) -> None:
                try:
                    if tab_button.cget("state") == "disabled":
                        return
                    if self.notebook.get() == name:
                        tab_button.configure(text_color=self._app_palette.get("on_accent", "#ffffff"))
                    else:
                        tab_button.configure(text_color=self._app_palette["secondary"])
                except Exception:
                    pass

            def on_leave(_event: object, name: str = tab_name, tab_button: tk.Widget = button) -> None:
                try:
                    if tab_button.cget("state") == "disabled":
                        tab_button.configure(text_color=self._app_palette["muted"])
                    elif self.notebook.get() == name:
                        tab_button.configure(text_color=self._app_palette.get("on_accent", "#ffffff"))
                    else:
                        tab_button.configure(text_color=self._app_palette["text"])
                except Exception:
                    pass

            button.bind("<Enter>", on_enter, add="+")
            button.bind("<Leave>", on_leave, add="+")
            self._tab_hover_bound.add(tab_name)

    def _refresh_inspect_button_state(self, *, busy: bool | None = None) -> None:
        if not hasattr(self, "inspect_btn"):
            return
        if busy is None:
            busy = bool(self._export_thread and self._export_thread.is_alive())

        valid_path = self._valid_h5ad_path()
        armed = (
            self._input_ready_for_inspect(valid_path)
            and not self._inspect_locked_after_press
            and not busy
        )
        self._configure_widget_state(self.inspect_btn, armed)
        if armed:
            self._style_secondary_action_button(self.inspect_btn)

    def _style_secondary_action_button(self, button: tk.Widget, *, borderless: bool = False) -> None:
        _style_secondary_button_for_palette(button, self._app_palette, borderless=borderless)

    def _export_is_running(self) -> bool:
        return bool(self._export_thread and self._export_thread.is_alive())

    def _is_export_stop_button(self, button: tk.Widget) -> bool:
        return button is getattr(self, "export_btn", None) and self._export_is_running()

    def _style_export_stop_button(self) -> None:
        button = getattr(self, "export_btn", None)
        if button is None:
            return
        try:
            button.configure(
                text="Stop Preview",
                command=self._cancel_operations,
                fg_color=self._app_palette["danger"],
                hover_color=self._app_palette["danger_hover"],
                text_color=self._app_palette["on_danger"],
                border_width=1,
                border_color=self._app_palette["danger"],
            )
        except Exception:
            pass

    def _bind_secondary_press_feedback(self, button: tk.Widget) -> None:
        def reset() -> None:
            if self._is_export_stop_button(button):
                self._style_export_stop_button()
            else:
                self._style_secondary_action_button(button)

        _bind_secondary_button_feedback(button, lambda: self._app_palette, reset_callback=reset)

    def _sync_right_panel_title_colors(self) -> None:
        for attr in ("runtime_title_label", "event_log_label"):
            widget = getattr(self, attr, None)
            if widget is None:
                continue
            try:
                widget.configure(state="normal", text_color=self._app_palette["secondary"])
            except Exception:
                pass

    def _set_inspect_loading(self, loading: bool) -> None:
        indicator = getattr(self, "inspect_loading_label", None)
        if indicator is None:
            return

        if loading:
            try:
                indicator.configure(text="⟳ Inspecting...", text_color=self._app_palette["secondary"])
                if not indicator.winfo_manager():
                    indicator.pack(side="left", padx=(12, 0))
                indicator.lift()
                self.status_var.set("Inspecting dataset...")
                self.update()
            except Exception:
                pass
            return

        try:
            indicator.pack_forget()
        except Exception:
            pass

    def _on_export_button(self) -> None:
        if self._export_thread and self._export_thread.is_alive():
            self._cancel_operations()
            return
        self._on_export()

    def _sync_export_button_state(self, *, busy: bool, enabled: bool) -> None:
        button = getattr(self, "export_btn", None)
        if button is None:
            return
        try:
            if busy:
                self._configure_widget_state(button, True)
                self._style_export_stop_button()
            else:
                button.configure(text="Build Viewer", command=self._on_export_button)
                self._configure_widget_state(button, enabled)
                if enabled:
                    self._style_secondary_action_button(button)
        except Exception:
            pass

    def _cancel_operations(self) -> None:
        self._cancel_requested.set()
        if self._server is not None:
            self._stop_server()
        if self._export_thread and self._export_thread.is_alive():
            self.status_var.set("Cancel requested")
            self._log("Cancel requested. Current export step will stop at the next safe checkpoint.")
        else:
            self.status_var.set("Ready")

    def _refresh_right_panel_gate(self) -> None:
        inspected = bool(self._has_inspected_input_file)
        side_outer = getattr(self, "side_outer", None)
        if side_outer is not None:
            try:
                if inspected:
                    side_outer.configure(**self._theme["card_frame"])
                else:
                    side_outer.configure(
                        fg_color=self._app_palette["hover_bg"],
                        border_color=self._app_palette["border"],
                    )
            except Exception:
                pass

        for widget in self._right_panel_widgets:
            self._configure_widget_state(widget, inspected)
        self._sync_right_panel_title_colors()

        if hasattr(self, "progress"):
            try:
                if inspected:
                    self.progress.configure(**self._theme["progress"])
                else:
                    self.progress.configure(
                        fg_color=self._app_palette["hover_bg"],
                        progress_color=self._app_palette["muted"],
                    )
            except Exception:
                pass
        if hasattr(self, "log_text"):
            try:
                self.log_text.configure(state="disabled")
                if inspected:
                    self.log_text.configure(**self._theme["textbox"])
                else:
                    self.log_text.configure(
                        fg_color=self._app_palette["hover_bg"],
                        text_color=self._app_palette["muted"],
                        border_color=self._app_palette["secondary"],
                    )
            except Exception:
                pass
        self._sync_runtime_chip()

    def _clear_inspection_metadata(self) -> None:
        self._inspected_h5ad_path = None
        self._inspected_coords_mode = None
        self._inspected_n_cells = None
        self._inspected_n_genes = None
        self._inspected_obs_cols = set()
        self._inspected_section_key_cols = []
        self._inspected_feature_names_by_modality = {}
        self._inspected_obsm_keys = []
        self._inspected_spatial_obs_cols = []
        self._inspected_layer_keys = []
        self._inspected_section_counts_by_column = {}
        self._inspected_section_values_by_column = {}
        self._feature_items_by_modality = {}
        self._active_feature_modality = "rna"
        if hasattr(self, "section_order_editor"):
            self.section_order_editor.clear()
        if hasattr(self, "feature_modality_combo"):
            self.feature_modality_var.set("rna")
            self.feature_modality_combo.configure(values=["rna"])
        if hasattr(self, "modalities_editor"):
            self.modalities_editor.set_choices([])
            self.modalities_editor.set_items([])
        if hasattr(self, "manual_genes_editor"):
            self.manual_genes_editor.set_choices([])
            self.manual_genes_editor.set_items([])
        if hasattr(self, "additional_colors_editor"):
            self.additional_colors_editor.set_choices([])
            self.additional_colors_editor.set_required_items([])
            self.additional_colors_editor.set_items([])
        if hasattr(self, "statistics_annotations_editor"):
            self.statistics_annotations_editor.set_choices([])
            self.statistics_annotations_editor.set_items([])
        if hasattr(self, "statistics_modalities_editor"):
            self.statistics_modalities_editor.set_choices([])
            self.statistics_modalities_editor.set_items([])
        if hasattr(self, "pseudobulk_replicate_combo"):
            self.pseudobulk_replicate_combo.configure(values=[""])
            self.pseudobulk_replicate_annotation_var.set("")
        self._sync_statistics_layer_choices()

    @staticmethod
    def _format_count(value: int) -> str:
        return f"{int(value):,}"

    @staticmethod
    def _limited_feature_names(values: object, max_features: int | None) -> list[str]:
        if values is None:
            return []
        try:
            iterable = values[:max_features] if max_features is not None else values
        except Exception:
            iterable = values
        names: list[str] = []
        seen: set[str] = set()
        for raw in iterable:
            name = str(raw).strip()
            if not name or name in seen:
                continue
            seen.add(name)
            names.append(name)
        return names

    @classmethod
    def _feature_names_from_modality_var(cls, var_obj: object, expected_rows: int, max_features: int | None) -> list[str]:
        if var_obj is None:
            return []
        try:
            import pandas as pd

            df = var_obj.copy() if hasattr(var_obj, "copy") else pd.DataFrame(var_obj)
            if int(getattr(df, "shape", [0])[0]) != int(expected_rows):
                return []
            if getattr(df.index, "name", None) is None:
                for candidate in ("protein", "gene", "feature", "name"):
                    if candidate in getattr(df, "columns", []):
                        df = df.set_index(candidate)
                        break
            return cls._limited_feature_names(df.index, max_features)
        except Exception:
            return []

    @staticmethod
    def _layer_keys_from_adata(adata) -> list[str]:
        layers = getattr(adata, "layers", None)
        if layers is None:
            return []
        try:
            keys = layers.keys()
        except Exception:
            return []
        layer_names: list[str] = []
        for key in keys:
            if key is None:
                continue
            value = str(key).strip()
            if value and value.lower() not in {"none", "off"}:
                layer_names.append(value)
        return layer_names

    @classmethod
    def _feature_names_by_modality_from_adata(
        cls,
        adata: object,
        *,
        max_features: int | None = 120000,
    ) -> dict[str, list[str]]:
        names_by_modality: dict[str, list[str]] = {
            "rna": cls._limited_feature_names(getattr(adata, "var_names", []), max_features)
        }
        obsm = getattr(adata, "obsm", None)
        uns = getattr(adata, "uns", {}) or {}
        if obsm is None or not hasattr(obsm, "keys"):
            return names_by_modality
        for raw_key in list(obsm.keys()):
            key = str(raw_key)
            if key in {"spatial", "deconvolutions"} or key.startswith("X_"):
                continue
            matrix = obsm[raw_key]
            shape = getattr(matrix, "shape", None)
            if shape is None or len(shape) != 2 or int(shape[1]) <= 0:
                continue
            modality_names = cls._feature_names_from_modality_var(
                uns.get(f"{key}_var"),
                int(shape[1]),
                max_features,
            )
            if modality_names:
                names_by_modality[key] = modality_names
        return names_by_modality

    def _save_active_feature_selection(self) -> None:
        if not hasattr(self, "manual_genes_editor"):
            return
        modality = (self._active_feature_modality or self.feature_modality_var.get().strip() or "rna").strip()
        self._feature_items_by_modality[modality] = self.manual_genes_editor.get_items()

    def _load_feature_selection_for_modality(self, modality: str) -> None:
        if not hasattr(self, "manual_genes_editor"):
            return
        choices = self._inspected_feature_names_by_modality.get(modality, [])
        self.manual_genes_editor.set_choices(choices)
        allowed = set(choices)
        saved = self._feature_items_by_modality.get(modality, [])
        if allowed:
            saved = [name for name in saved if name in allowed]
        self.manual_genes_editor.set_items(saved)

    def _selected_embed_modalities(self) -> list[str]:
        editor = getattr(self, "modalities_editor", None)
        if editor is not None and hasattr(editor, "get_items"):
            return self._merge_unique(editor.get_items())
        return []

    def _feature_modality_choices(self) -> list[str]:
        inspected = list(self._inspected_feature_names_by_modality)
        selected = self._selected_embed_modalities()
        if selected:
            if inspected:
                allowed = set(inspected)
                choices = [modality for modality in selected if modality in allowed]
            else:
                choices = selected
        else:
            choices = inspected
        return choices or (["rna"] if not inspected else inspected[:1])

    def _sync_feature_modality_choices(self) -> None:
        if not hasattr(self, "feature_modality_combo"):
            return
        self._save_active_feature_selection()
        choices = self._feature_modality_choices()
        current = self.feature_modality_var.get().strip()
        if current not in choices:
            current = "rna" if "rna" in choices else choices[0]
        self._active_feature_modality = current
        self.feature_modality_var.set(current)
        self.feature_modality_combo.configure(values=choices)
        self._load_feature_selection_for_modality(current)
        self._sync_statistics_modality_choices()
        self._update_export_estimate()

    def _sync_statistics_modality_choices(self) -> None:
        editor = getattr(self, "statistics_modalities_editor", None)
        if editor is None:
            return
        choices = self._feature_modality_choices()
        allowed = set(choices)
        editor.set_choices(choices)
        editor.set_items([modality for modality in editor.get_items() if modality in allowed])

    def _on_embed_modalities_changed(self) -> None:
        self._sync_feature_modality_choices()

    def _set_feature_modalities(self, names_by_modality: dict[str, list[str]]) -> None:
        self._save_active_feature_selection()
        cleaned = {
            str(name).strip(): self._limited_feature_names(features, None)
            for name, features in names_by_modality.items()
            if str(name).strip()
        }
        if not cleaned:
            cleaned = {"rna": []}
        choices = list(cleaned.keys())
        self._inspected_feature_names_by_modality = cleaned
        for modality, selected in list(self._feature_items_by_modality.items()):
            allowed = set(cleaned.get(modality, []))
            if not allowed:
                self._feature_items_by_modality.pop(modality, None)
            else:
                self._feature_items_by_modality[modality] = [name for name in selected if name in allowed]
        if hasattr(self, "modalities_editor"):
            self.modalities_editor.set_choices(choices)
            existing = [modality for modality in self.modalities_editor.get_items() if modality in set(choices)]
            if not existing:
                existing = ["rna"] if "rna" in cleaned else [choices[0]]
            self.modalities_editor.set_items(existing)
        self._sync_feature_modality_choices()

    def _on_feature_modality_changed(self) -> None:
        next_modality = self.feature_modality_var.get().strip() or "rna"
        if next_modality == self._active_feature_modality:
            return
        self._save_active_feature_selection()
        self._active_feature_modality = next_modality
        self._load_feature_selection_for_modality(next_modality)
        self._update_export_estimate()

    def _on_also_karospace_toggled(self) -> None:
        # A .karospace can only be built from a sidecar bundle, so enabling the
        # package implies sidecar storage.
        if not hasattr(self, "also_karospace_var"):
            return
        if self.also_karospace_var.get() and self.feature_storage_var.get().strip().lower() != "sidecar":
            self.feature_storage_var.set("sidecar")

    def _on_feature_storage_changed(self) -> None:
        # Embedded storage produces no sidecar, so the .karospace package is not
        # available; keep the checkbox consistent with the chosen storage.
        if not hasattr(self, "also_karospace_var"):
            return
        if self.feature_storage_var.get().strip().lower() != "sidecar":
            self.also_karospace_var.set(False)

    def _selected_features_by_modality(self) -> dict[str, list[str]]:
        self._save_active_feature_selection()
        selected: dict[str, list[str]] = {}
        for modality in self._feature_modality_choices():
            values = self._merge_unique(self._feature_items_by_modality.get(modality, []))
            if values:
                selected[modality] = values
        return selected

    @staticmethod
    def _section_values_from_adata(adata: object, column: str) -> list[str]:
        obs = getattr(adata, "obs", None)
        if obs is None or column not in getattr(obs, "columns", []):
            return []
        try:
            raw_values = obs[column].dropna().drop_duplicates().tolist()
        except Exception:
            raw_values = obs[column].tolist()
        seen: set[str] = set()
        values: list[str] = []
        for raw in raw_values:
            value = str(raw).strip()
            if not value or value.lower() == "nan" or value in seen:
                continue
            seen.add(value)
            values.append(value)
        return values

    @classmethod
    def _eligible_section_key_columns(
        cls,
        adata: object,
        obs_cols: list[str] | tuple[str, ...],
        *,
        max_unique: int = 500,
    ) -> list[str]:
        obs = getattr(adata, "obs", None)
        if obs is None:
            return []
        eligible: list[str] = []
        for column in obs_cols:
            name = str(column)
            if name not in getattr(obs, "columns", []):
                continue
            try:
                n_unique = int(obs[name].dropna().nunique())
            except Exception:
                values = cls._section_values_from_adata(adata, name)
                n_unique = len(values)
            if 0 < n_unique < max_unique:
                eligible.append(name)
        return eligible

    @staticmethod
    def _obsm_keys_from_adata(adata: object) -> list[str]:
        obsm = getattr(adata, "obsm", None)
        keys = []
        if obsm is not None and hasattr(obsm, "keys"):
            keys = [str(key) for key in obsm.keys()]
        seen: set[str] = set()
        ordered: list[str] = []
        for key in keys:
            if not key or key in seen:
                continue
            seen.add(key)
            ordered.append(key)
        if "spatial" in ordered:
            ordered = ["spatial", *[key for key in ordered if key != "spatial"]]
        return ordered or ["spatial"]

    @staticmethod
    def _spatial_obs_columns_from_adata(adata: object, obs_cols: list[str] | tuple[str, ...]) -> list[str]:
        obs = getattr(adata, "obs", None)
        if obs is None:
            return []
        columns: list[str] = []
        for column in obs_cols:
            name = str(column)
            if name not in getattr(obs, "columns", []):
                continue
            dtype = getattr(obs[name], "dtype", None)
            if getattr(dtype, "kind", "") in {"i", "u", "f"}:
                columns.append(name)
        return columns

    def _get_section_values_for_order(self, column: str) -> list[str]:
        if not column or column not in self._inspected_obs_cols:
            return []
        if self._inspected_section_key_cols and column not in self._inspected_section_key_cols:
            return []
        cached = self._inspected_section_values_by_column.get(column)
        if cached is not None:
            return cached
        path = self._inspected_h5ad_path
        if path is None or not path.exists():
            return []

        adata = None
        try:
            adata = self._read_adata_for_feature_ops(path)
            values = self._section_values_from_adata(adata, column)
            self._inspected_section_values_by_column[column] = values
            return values
        except Exception:
            return []
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    def _refresh_section_order_values(self, *, preserve_current: bool = True) -> None:
        editor = getattr(self, "section_order_editor", None)
        if editor is None:
            return
        if not self._has_inspected_input_file:
            editor.clear()
            return
        values = self._get_section_values_for_order(self.section_groupby_var.get().strip())
        editor.set_values(values, preserve_current=preserve_current)

    def _on_section_groupby_changed(self) -> None:
        self._sync_required_section_metadata()
        self._refresh_section_order_values(preserve_current=False)
        self._update_export_estimate()

    def _required_cell_annotation_items(self) -> list[str]:
        return self._merge_unique([self.initial_color_var.get().strip()])

    def _sync_required_cell_annotations(self) -> None:
        editor = getattr(self, "additional_colors_editor", None)
        if editor is None:
            return
        editor.set_required_items(self._required_cell_annotation_items())

    def _required_section_metadata_items(self) -> list[str]:
        required: list[str] = []
        outline_by = self.outline_by_var.get().strip()
        if outline_by:
            required.append(outline_by)
        return self._merge_unique(required)

    def _sync_required_section_metadata(self) -> None:
        editor = getattr(self, "section_metadata_editor", None)
        if editor is None:
            return
        editor.set_required_items(self._required_section_metadata_items())

    def _sync_statistics_layer_choices(self, *, preserve_current: bool = False) -> None:
        layer_keys = self._inspected_layer_keys
        counts_values = self._merge_unique(["none"], layer_keys)
        normalized_values = self._merge_unique(["off"], layer_keys)

        counts_current = self.statistics_counts_layer_var.get().strip()
        if preserve_current and counts_current and counts_current not in counts_values:
            counts_values = self._merge_unique(counts_values, [counts_current])
        if hasattr(self, "statistics_counts_layer_combo"):
            self.statistics_counts_layer_combo.configure(values=counts_values)
        if not preserve_current and counts_current not in counts_values:
            self.statistics_counts_layer_var.set("counts" if "counts" in layer_keys else "none")

        normalized_current = self.statistics_normalized_layer_var.get().strip()
        if preserve_current and normalized_current and normalized_current not in normalized_values:
            normalized_values = self._merge_unique(normalized_values, [normalized_current])
        if hasattr(self, "statistics_normalized_layer_combo"):
            self.statistics_normalized_layer_combo.configure(values=normalized_values)
        if not preserve_current and normalized_current not in normalized_values:
            self.statistics_normalized_layer_var.set("off")

    def _estimate_exported_gene_count(self) -> tuple[int | None, str | None]:
        if self._inspected_n_genes is None and not self._inspected_feature_names_by_modality:
            return None, "dataset gene count is unknown"

        selected_by_modality = self._selected_features_by_modality()
        genes = self._merge_unique(*selected_by_modality.values())
        count = len(self._merge_unique(genes))
        if count <= 0:
            return None, "add at least one feature"
        total_features = sum(len(features) for features in self._inspected_feature_names_by_modality.values())
        if total_features <= 0:
            total_features = self._inspected_n_genes or count
        return min(total_features, count), None

    def _get_section_counts_for_estimate(self, column: str) -> list[int] | None:
        if not column or column not in self._inspected_obs_cols:
            return None
        cached = self._inspected_section_counts_by_column.get(column)
        if cached is not None:
            return cached
        path = self._inspected_h5ad_path
        if path is None or not path.exists():
            return None

        adata = None
        try:
            ad_mod = _get_anndata()
            try:
                adata = ad_mod.read_h5ad(path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(path)
            counts = [int(value) for value in adata.obs[column].value_counts(dropna=False).tolist()]
            self._inspected_section_counts_by_column[column] = counts
            return counts
        except Exception:
            return None
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    def _estimate_exported_cell_count(self) -> tuple[int | None, str | None]:
        total_cells = self._inspected_n_cells
        if total_cells is None:
            return None, "dataset cell count is unknown"

        downsample_text = self.downsample_var.get().strip()
        if not downsample_text:
            return total_cells, None
        try:
            downsample = int(downsample_text)
        except ValueError:
            return None, "enter a valid downsample value"
        if downsample <= 0:
            return None, "enter a positive downsample value"

        groupby = self.section_groupby_var.get().strip()
        if not groupby:
            return min(total_cells, downsample), None
        counts = self._get_section_counts_for_estimate(groupby)
        if counts:
            return sum(min(count, downsample) for count in counts), None
        return min(total_cells, downsample), f"section counts unavailable for '{groupby}'"

    def _update_export_estimate(self) -> None:
        if not hasattr(self, "export_estimate_summary_label"):
            return

        if not self._has_inspected_input_file or self._inspected_n_cells is None or self._inspected_n_genes is None:
            self.export_estimate_summary_label.configure(
                text="Inspect a dataset to estimate exported cells, genes, elements, and HTML size.",
                text_color=self._app_palette["text"],
            )
            self.export_estimate_warning_label.configure(text="", text_color=self._app_palette["muted"])
            return

        cells, cell_note = self._estimate_exported_cell_count()
        genes, gene_note = self._estimate_exported_gene_count()
        notes = [note for note in (cell_note, gene_note) if note]
        if cells is None or genes is None:
            note = "; ".join(notes) if notes else "complete the export settings"
            self.export_estimate_summary_label.configure(
                text=f"Estimate waiting for valid settings: {note}.",
                text_color=self._app_palette["text"],
            )
            self.export_estimate_warning_label.configure(text="", text_color=self._app_palette["muted"])
            return

        elements = int(cells) * int(genes)
        predicted_mb = self._HTML_SIZE_MB_PER_ELEMENT * elements
        summary = (
            f"Export estimate: {self._format_count(cells)} cells x {self._format_count(genes)} genes = "
            f"{self._format_count(elements)} elements. Predicted HTML size: {predicted_mb:,.1f} MB. "
            "Expected if features are selected at random. HTML file size will drastically increase if top features are selected (dense matrix)."
        )
        if notes:
            summary = f"{summary} Note: {'; '.join(notes)}."
        self.export_estimate_summary_label.configure(text=summary, text_color=self._app_palette["text"])

        if elements > self._HTML_SIZE_WARNING_ELEMENTS:
            self.export_estimate_warning_label.configure(
                text=(
                    "Warning: more than 150 million exported entries. Builder will probably generate an HTML file "
                    ">500 MB, and the browser may not be able to display such a big HTML file."
                ),
                text_color=self._app_palette["danger"],
            )
        else:
            self.export_estimate_warning_label.configure(text="", text_color=self._app_palette["muted"])

    def _refresh_input_gate(self) -> None:
        valid_path = self._valid_h5ad_path()
        has_valid_input = valid_path is not None
        if has_valid_input != self._has_valid_input_file:
            self._has_valid_input_file = has_valid_input
        input_changed = (
            valid_path is not None
            and self._inspected_h5ad_path is not None
            and self._inspected_h5ad_path != valid_path
        )
        if not has_valid_input or input_changed:
            self._has_inspected_input_file = False
            if not self._loading_parameter_state:
                self._clear_inspection_metadata()
                self._spatialdata_table_choices = []
                self._spatialdata_table_source_path = None
                if not self._is_zarr_path(valid_path) or input_changed:
                    self._syncing_spatialdata_table = True
                    try:
                        self.spatialdata_table_var.set("")
                    finally:
                        self._syncing_spatialdata_table = False
        if valid_path is not None and self._spatialdata_table_source_path not in {None, valid_path}:
            self._spatialdata_table_choices = []
            self._spatialdata_table_source_path = None
            self._syncing_spatialdata_table = True
            try:
                self.spatialdata_table_var.set("")
            finally:
                self._syncing_spatialdata_table = False
        if has_valid_input and self._is_zarr_path(valid_path):
            self._autoselect_spatialdata_table(valid_path)

        busy = bool(self._export_thread and self._export_thread.is_alive())
        enabled = has_valid_input and self._has_inspected_input_file and not busy
        self._set_widgets_visible(self._inspection_gated_widgets, self._has_inspected_input_file)
        self._sync_spatialdata_table_selector(valid_path, busy=busy)

        for attr in (
            "coords_menu",
            "spatial_key_combo",
            "spatial_x_combo",
            "spatial_y_combo",
            "groupby_combo",
            "section_order_editor",
            "color_combo",
            "numba_jit_check",
            "outline_combo",
            "title_entry",
            "downsample_entry",
        ):
            self._configure_widget_state(getattr(self, attr, None), enabled)
        if enabled:
            for attr in ("spatial_key_combo", "spatial_x_combo", "spatial_y_combo", "groupby_combo", "color_combo", "outline_combo"):
                try:
                    getattr(self, attr).configure(state="readonly")
                except Exception:
                    pass

        self._sync_export_button_state(busy=busy, enabled=enabled)
        self._configure_widget_state(getattr(self, "inspect_btn", None), self._input_ready_for_inspect(valid_path) and not busy)
        self._refresh_inspect_button_state(busy=busy)
        self._configure_widget_state(getattr(self, "open_output_btn", None), self._has_inspected_input_file and not busy)
        self._configure_widget_state(getattr(self, "open_viewer_btn", None), self._has_inspected_input_file and not busy)

        self._refresh_tab_gate()
        self._refresh_right_panel_gate()
        self._update_export_estimate()

    def _sync_runtime_chip(self) -> None:
        if not hasattr(self, "runtime_chip_label"):
            return
        status = self.status_var.get().strip().lower()
        if status.startswith("export running"):
            fg = self._app_palette["accent"]
            tc = self._app_palette.get("on_accent", "#ffffff")
            self.runtime_chip_label.configure(
                fg_color=fg,
                text_color=tc,
                border_width=0,
                width=104,
            )
            self._start_runtime_chip_animation()
            return
        elif status.startswith("export complete"):
            self._stop_runtime_chip_animation()
            text = "COMPLETE"
            fg = "#1f8f5f"
            tc = "#ffffff"
        elif status.startswith("serving on"):
            self._stop_runtime_chip_animation()
            text = "SERVING"
            fg = self._app_palette["accent_strong"]
            tc = self._app_palette.get("on_accent", "#ffffff")
        elif status.startswith("export failed"):
            self._stop_runtime_chip_animation()
            text = "FAILED"
            fg = "#bf2f5e"
            tc = "#ffffff"
        else:
            self._stop_runtime_chip_animation()
            text = "READY"
            fg = self._theme["progress"]["fg_color"]
            tc = self._app_palette["text"]
        self.runtime_chip_label.configure(
            text=text,
            fg_color=fg,
            text_color=tc,
            border_width=0,
        )

    def _start_runtime_chip_animation(self) -> None:
        if self._runtime_chip_animation_after_id is None:
            self._runtime_chip_animation_step = 0
            self._animate_runtime_chip()

    def _stop_runtime_chip_animation(self) -> None:
        after_id = self._runtime_chip_animation_after_id
        self._runtime_chip_animation_after_id = None
        if after_id is None:
            return
        try:
            self.after_cancel(after_id)
        except Exception:
            pass

    def _animate_runtime_chip(self) -> None:
        if not hasattr(self, "runtime_chip_label"):
            self._runtime_chip_animation_after_id = None
            return
        if not self.status_var.get().strip().lower().startswith("export running"):
            self._runtime_chip_animation_after_id = None
            self._sync_runtime_chip()
            return
        dot_counts = (0, 1, 2, 3, 2, 1)
        count = dot_counts[self._runtime_chip_animation_step % len(dot_counts)]
        self._runtime_chip_animation_step += 1
        self.runtime_chip_label.configure(text=f"Running{'.' * count}")
        self._runtime_chip_animation_after_id = self.after(260, self._animate_runtime_chip)

    def _build_variables(self) -> None:
        self.h5ad_var = tk.StringVar()
        self.outdir_var = tk.StringVar(value=str(Path.cwd()))
        self.output_html_var = tk.StringVar(value="karospace.html")
        self.coords_var = tk.StringVar(value="auto")
        self.spatial_key_var = tk.StringVar(value="spatial")
        self.spatial_x_var = tk.StringVar()
        self.spatial_y_var = tk.StringVar()
        self.spatialdata_table_var = tk.StringVar()
        self.section_groupby_var = tk.StringVar(value="sample_id")
        self.section_order_var = tk.StringVar()
        self.metadata_value_order_var = tk.StringVar()
        self.metadata_max_columns_var = tk.StringVar()
        self.initial_color_var = tk.StringVar(value="leiden")
        self.title_var = tk.StringVar(value="KaroSpace")
        self.theme_var = tk.StringVar(value="dark")
        self.numba_jit_var = tk.BooleanVar(value=False)
        self.outline_by_var = tk.StringVar()
        self.metadata_labels_var = tk.StringVar()
        self.viewer_info_html_var = tk.StringVar()
        self.viewer_info_html_file_var = tk.StringVar()
        self.tutorial_var = tk.BooleanVar(value=False)
        self.embed_reproducibility_info_var = tk.BooleanVar(value=True)

        self.feature_encoding_var = tk.StringVar(value="auto")
        self.feature_value_encoding_var = tk.StringVar(value="uint16")
        self.feature_storage_var = tk.StringVar(value="sidecar")
        self.also_karospace_var = tk.BooleanVar(value=True)
        self.features_list_var = tk.StringVar()
        self.feature_manifest_path_var = tk.StringVar()
        self.feature_sidecar_shard_size_var = tk.StringVar(value="256")
        self.feature_sparse_zero_threshold_var = tk.StringVar(value="0.8")
        self.feature_modality_var = tk.StringVar(value="rna")
        self.min_panel_size_var = tk.StringVar(value="150")
        self.spot_size_var = tk.StringVar(value="auto")
        self.neighbor_stats_annotations_var = tk.StringVar(value="auto")
        self.neighbor_permutations_var = tk.StringVar(value="20")
        self.neighbor_stats_seed_var = tk.StringVar(value="0")
        self.statistics_contrast_categories_var = tk.StringVar()
        self.statistics_counts_layer_var = tk.StringVar(value="counts")
        self.statistics_normalization_var = tk.StringVar(value="RC")
        self.statistics_scale_factor_var = tk.StringVar(value="10000")
        self.statistics_normalized_layer_var = tk.StringVar(value="off")
        self.statistics_min_cell_counts_var = tk.StringVar(value="0")
        self.statistics_min_feature_counts_var = tk.StringVar(value="0")
        self.statistics_n_cpus_var = tk.StringVar(value="1")
        self.wilcoxon_enabled_var = tk.BooleanVar(value=True)
        self.wilcoxon_mode_var = tk.StringVar(value="auto")
        self.wilcoxon_runtime_limit_var = tk.StringVar(value="00:30:00")
        self.wilcoxon_min_cells_per_group_var = tk.StringVar(value="20")
        self.wilcoxon_min_pct_expressed_var = tk.StringVar(value="0")
        self.wilcoxon_p_adjust_method_var = tk.StringVar(value="fdr_bh")
        self.wilcoxon_padj_cutoff_var = tk.StringVar(value="0.05")
        self.wilcoxon_log2fc_cutoff_var = tk.StringVar(value="1")
        self.wilcoxon_embed_top_n_per_comparison_var = tk.StringVar(value="2")
        self.wilcoxon_top_n_per_category_var = tk.StringVar(value="300")
        self.pseudobulk_enabled_var = tk.BooleanVar(value=False)
        self.pseudobulk_replicate_annotation_var = tk.StringVar()
        self.pseudobulk_min_cells_per_pseudobulk_var = tk.StringVar(value="20")
        self.pseudobulk_min_replicates_var = tk.StringVar(value="2")
        self.pseudobulk_min_pct_expressed_var = tk.StringVar(value="0")
        self.pseudobulk_p_adjust_method_var = tk.StringVar(value="fdr_bh")
        self.pseudobulk_padj_cutoff_var = tk.StringVar(value="0.05")
        self.pseudobulk_log2fc_cutoff_var = tk.StringVar(value="1")
        self.pseudobulk_deseq2_fit_type_var = tk.StringVar(value="parametric")
        self.pseudobulk_embed_top_n_per_comparison_var = tk.StringVar(value="2")
        self.pathway_enabled_var = tk.BooleanVar(value=False)
        self.pathway_gmt_var = tk.StringVar()
        self.pathway_organism_var = tk.StringVar(value="Mouse")
        self.pathway_top_n_var = tk.StringVar(value="10")
        self.pathway_min_overlap_var = tk.StringVar(value="3")
        self.pathway_gsea_permutations_var = tk.StringVar(value="100")
        self.interaction_markers_enabled_var = tk.BooleanVar(value=True)
        self.interaction_markers_top_targets_var = tk.StringVar(value="5")
        self.interaction_markers_top_features_var = tk.StringVar(value="20")
        self.interaction_markers_min_cells_var = tk.StringVar(value="30")
        self.interaction_markers_min_neighbors_var = tk.StringVar(value="1")
        self.section_rotations_var = tk.StringVar()
        self.deconvolutions_var = tk.StringVar()
        self.feature_correlation_top_n_var = tk.StringVar(value="5")
        self.spatial_variable_features_n_var = tk.StringVar(value="20")
        self.scalebar_unit_var = tk.StringVar(value="μm")
        self.section_images_var = tk.StringVar()
        self.section_images_max_px_var = tk.StringVar(value="4096")
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
        root.grid(row=0, column=0, sticky="nsew", padx=16, pady=16)
        root.columnconfigure(0, weight=3, minsize=730, uniform="main_columns")
        root.columnconfigure(1, weight=2, minsize=430, uniform="main_columns")
        root.rowconfigure(0, weight=1)

        controls = ctk.CTkFrame(root, width=730, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", controls)
        controls.grid(row=0, column=0, sticky="nsew", padx=(0, 12))
        controls.columnconfigure(0, weight=1)
        controls.rowconfigure(1, weight=1)
        controls_outer = controls

        side = ctk.CTkFrame(root, width=430, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", side)
        side.grid(row=0, column=1, sticky="nsew")

        controls_inner = ctk.CTkFrame(controls, **self._theme["sub_frame"])
        self._register_theme_widget("sub_frame", controls_inner)
        controls_inner.grid(row=0, column=0, sticky="ew", padx=22, pady=(22, 8))
        controls = controls_inner

        side_inner = ctk.CTkFrame(side, **self._theme["sub_frame"])
        self.side_outer = side
        self.side_panel = side_inner
        self._register_theme_widget("sub_frame", side_inner)
        side_inner.pack(fill="both", expand=True, padx=22, pady=22)
        side = side_inner

        controls.columnconfigure(1, weight=1)
        side.columnconfigure(0, weight=1)
        side.rowconfigure(3, weight=1)

        hero = ctk.CTkFrame(controls, **self._theme["hero_card"])
        self._register_theme_widget("hero_card", hero)
        hero.grid(row=0, column=0, columnspan=3, sticky="ew", pady=(0, 16))
        hero.columnconfigure(0, weight=1)

        hero_inner = self._make_sub_frame(hero)
        hero_inner.grid(row=0, column=0, sticky="ew", padx=18, pady=18)
        hero_inner.columnconfigure(0, weight=1)

        hero_left = self._make_sub_frame(hero_inner)
        hero_left.grid(row=0, column=0, sticky="w")
        self._section_label(hero_left, "DESKTOP BUILDER").pack(anchor="w")
        self._hero_label(hero_left, "KaroSpaceBuilder").pack(anchor="w", pady=(4, 0))
        self._subheader_label(
            hero_left,
            "Export AnnData into a static KaroSpace viewer bundle with guided inputs and inspected field pickers.",
        ).pack(anchor="w", pady=(4, 0))

        tabs_toolbar = self._make_sub_frame(controls)
        self.tabs_toolbar = tabs_toolbar
        tabs_toolbar.grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 8))
        tabs_toolbar.columnconfigure(0, weight=1)

        self.theme_toggle_btn = self._secondary_button(tabs_toolbar, "☀", self._toggle_theme, width=42)
        self.theme_toggle_btn.grid(row=0, column=1, sticky="e", padx=(0, 8))
        self.guide_btn = self._secondary_button(tabs_toolbar, "?", self._show_guide_popup, width=42)
        self.guide_btn.grid(row=0, column=2, sticky="e")

        self._path_field(controls, 2, "Input file/store", self.h5ad_var, choose_file=True)
        self._build_parameter_state_buttons(controls, 3)
        self._build_spatialdata_table_selector(controls, 4)

        notebook = ctk.CTkTabview(controls_outer, command=self._on_tab_selected, width=1, **self._theme["tabview"])
        self.notebook = notebook
        self._last_enabled_tab = "General"
        self._register_theme_widget("tabview", notebook)
        notebook.grid(row=1, column=0, sticky="ew", padx=18)
        self._inspection_gated_widgets.append(notebook)
        controls_outer.rowconfigure(1, weight=0)

        notebook.add("General")
        notebook.add("Metadata")
        notebook.add("Features")
        notebook.add("Statistics")
        notebook.add("Neighborhoods")
        notebook.add("Connections")
        notebook.add("Viewer")
        notebook.add("Sections")
        input_tab = notebook.tab("General")
        metadata_tab = notebook.tab("Metadata")
        features_tab = notebook.tab("Features")
        statistics_tab = notebook.tab("Statistics")
        neighborhoods_tab = notebook.tab("Neighborhoods")
        connections_tab = notebook.tab("Connections")
        viewer_tab = notebook.tab("Viewer")
        overlays_tab = notebook.tab("Sections")
        for tab in (
            input_tab,
            metadata_tab,
            features_tab,
            statistics_tab,
            neighborhoods_tab,
            connections_tab,
            viewer_tab,
            overlays_tab,
        ):
            self._register_theme_widget("sub_frame", tab)

        input_tab_outer = input_tab
        metadata_tab_outer = metadata_tab
        features_tab_outer = features_tab
        statistics_tab_outer = statistics_tab
        neighborhoods_tab_outer = neighborhoods_tab
        connections_tab_outer = connections_tab
        viewer_tab_outer = viewer_tab
        overlays_tab_outer = overlays_tab
        for outer_tab in (
            input_tab_outer,
            metadata_tab_outer,
            features_tab_outer,
            statistics_tab_outer,
            neighborhoods_tab_outer,
            connections_tab_outer,
            viewer_tab_outer,
            overlays_tab_outer,
        ):
            outer_tab.columnconfigure(0, weight=1)
            outer_tab.rowconfigure(0, weight=1)

        input_tab = self._make_sub_frame(input_tab_outer)
        input_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        metadata_tab = self._make_sub_frame(metadata_tab_outer)
        metadata_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        features_tab = self._make_sub_frame(features_tab_outer)
        features_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        statistics_tab = self._make_sub_frame(statistics_tab_outer)
        statistics_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        neighborhoods_tab = self._make_sub_frame(neighborhoods_tab_outer)
        neighborhoods_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        connections_tab = self._make_sub_frame(connections_tab_outer)
        connections_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        viewer_tab = self._make_sub_frame(viewer_tab_outer)
        viewer_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)
        overlays_tab = self._make_sub_frame(overlays_tab_outer)
        overlays_tab.grid(row=0, column=0, sticky="nsew", padx=22, pady=18)

        input_tab.columnconfigure(1, weight=1)
        self._section_label(input_tab, "INPUT / OUTPUT").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 8))
        self._divider(input_tab, height=1).grid(row=1, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        row = 2
        row = self._option_row(
            input_tab,
            row,
            "Coordinates",
            widget=self._coords_dropdown(input_tab),
            hint="auto | obsm:<spatial-key> | obs:<spatial-x>,<spatial-y>. Mirrors --spatial-key/--spatial-x/--spatial-y.",
            gated=True,
        )
        spatial_row = self._make_sub_frame(input_tab)
        for column in range(3):
            spatial_row.columnconfigure(column, weight=1, uniform="spatial_coordinate")

        spatial_key_cell = self._make_sub_frame(spatial_row)
        spatial_key_cell.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        spatial_key_cell.columnconfigure(0, weight=1)
        self._body_label(spatial_key_cell, "Spatial key").grid(row=0, column=0, sticky="w", pady=(0, 3))
        self.spatial_key_combo = ctk.CTkComboBox(
            spatial_key_cell,
            variable=self.spatial_key_var,
            values=["spatial"],
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.spatial_key_combo)
        self.spatial_key_combo.grid(row=1, column=0, sticky="ew")

        spatial_x_cell = self._make_sub_frame(spatial_row)
        spatial_x_cell.grid(row=0, column=1, sticky="ew", padx=(0, 8))
        spatial_x_cell.columnconfigure(0, weight=1)
        self._body_label(spatial_x_cell, "X").grid(row=0, column=0, sticky="w", pady=(0, 3))
        self.spatial_x_combo = ctk.CTkComboBox(
            spatial_x_cell,
            variable=self.spatial_x_var,
            values=[""],
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.spatial_x_combo)
        self.spatial_x_combo.grid(row=1, column=0, sticky="ew")

        spatial_y_cell = self._make_sub_frame(spatial_row)
        spatial_y_cell.grid(row=0, column=2, sticky="ew")
        spatial_y_cell.columnconfigure(0, weight=1)
        self._body_label(spatial_y_cell, "Y").grid(row=0, column=0, sticky="w", pady=(0, 3))
        self.spatial_y_combo = ctk.CTkComboBox(
            spatial_y_cell,
            variable=self.spatial_y_var,
            values=[""],
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.spatial_y_combo)
        self.spatial_y_combo.grid(row=1, column=0, sticky="ew")
        row = self._option_row(
            input_tab,
            row,
            "Spatial coordinates",
            widget=spatial_row,
            hint="Choose an inspected obsm key, or choose both X/Y obs columns to build spatial coordinates.",
            gated=True,
        )
        self.modalities_editor = SearchableListEditor(
            input_tab,
            label="Modality",
            height=4,
            help_text="Modalities to embed/export. The Features tab modality dropdown is limited to this selection.",
            palette=self._app_palette,
            on_change=self._on_embed_modalities_changed,
        )
        self.modalities_editor.grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        self._inspection_gated_widgets.append(self.modalities_editor)
        row += 1
        row = self._option_row(
            input_tab,
            row,
            "Section key",
            widget=self._groupby_dropdown(input_tab),
            hint="obs column used to split sections. Leave empty for a single-sample H5AD.",
            gated=True,
        )
        section_order_row = self._make_sub_frame(input_tab)
        section_order_row.columnconfigure(0, weight=1)
        self.section_order_editor = SectionOrderEditor(
            section_order_row,
            variable=self.section_order_var,
            palette=self._app_palette,
            on_change=self._update_export_estimate,
        )
        self.section_order_editor.grid(row=0, column=0, sticky="ew")
        row = self._option_row(
            input_tab,
            row,
            "Section order",
            widget=section_order_row,
            gated=True,
        )
        row = self._option_row(
            input_tab,
            row,
            "Main cell annotation",
            widget=self._color_dropdown(input_tab),
            gated=True,
        )
        row = self._option_row(
            input_tab,
            row,
            "Outline by",
            widget=self._outline_dropdown(input_tab),
            gated=True,
        )
        self._section_label(input_tab, "OUTPUT").grid(row=row, column=0, columnspan=3, sticky="w", pady=(10, 8))
        row += 1
        self._divider(input_tab, height=1).grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        row += 1
        row = self._path_field(input_tab, row, "Output directory", self.outdir_var, choose_file=False, gated=True)
        output_html_row = self._make_sub_frame(input_tab)
        output_html_entry = ctk.CTkEntry(output_html_row, textvariable=self.output_html_var, **self._theme["entry"])
        self._register_entry_widget(output_html_entry)
        output_html_entry.pack(side="left", fill="x", expand=True)
        output_html_button = self._secondary_button(
            output_html_row,
            "Browse",
            self._choose_output_file,
            width=96,
        )
        output_html_button.pack(side="left", padx=(8, 0))
        row = self._option_row(
            input_tab,
            row,
            "Output path",
            widget=output_html_row,
            gated=True,
        )
        self.title_entry = self._entry(input_tab, self.title_var)
        row = self._option_row(
            input_tab,
            row,
            "Viewer title",
            widget=self.title_entry,
            hint="Title shown in the exported HTML.",
            gated=True,
        )
        self._section_label(input_tab, "RUNTIME").grid(row=row, column=0, columnspan=3, sticky="w", pady=(10, 8))
        row += 1
        self._divider(input_tab, height=1).grid(row=row, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        row += 1
        runtime_mode_row = self._make_sub_frame(input_tab)
        self.runtime_mode_row = runtime_mode_row
        self.numba_jit_check = ctk.CTkCheckBox(
            runtime_mode_row,
            text="Performance mode (enable numba JIT)",
            variable=self.numba_jit_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.numba_jit_check)
        self.numba_jit_check.pack(side="left")
        row = self._option_row(
            input_tab,
            row,
            "Runtime mode",
            widget=runtime_mode_row,
            hint="Faster on some datasets. If unstable in the desktop app, turn this off.",
            gated=True,
        )
        viewer_tab.columnconfigure(1, weight=1)
        self._section_label(viewer_tab, "VIEWER LAYOUT").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        self._subheader_label(
            viewer_tab,
            "Configure exported viewer layout and cell-level display volume.",
        ).grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 8))
        self._divider(viewer_tab, height=1).grid(row=2, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        downsample_container = self._make_sub_frame(viewer_tab)
        downsample_entry = ctk.CTkEntry(downsample_container, textvariable=self.downsample_var, width=100, **self._theme["entry"])
        self.downsample_entry = downsample_entry
        self._register_entry_widget(downsample_entry)
        downsample_entry.pack(side="left")
        self._body_label(downsample_container, "cells per section (blank = all)").pack(side="left", padx=(8, 0))
        self._option_row(
            viewer_tab,
            3,
            "Downsample",
            widget=downsample_container,
            hint="Maps to export_to_html(downsample=...).",
            gated=True,
        )

        metadata_tab.columnconfigure(0, weight=0)
        metadata_tab.columnconfigure(1, weight=1)
        self._section_label(metadata_tab, "CELL ANNOTATIONS AND SECTION METADATA").grid(row=0, column=0, columnspan=3, sticky="w", pady=(0, 6))
        self._subheader_label(
            metadata_tab,
            "Pick cell annotations and section metadata from inspected obs fields.",
        ).grid(row=1, column=0, columnspan=3, sticky="w", pady=(0, 8))
        self._divider(metadata_tab, height=1).grid(row=2, column=0, columnspan=3, sticky="ew", pady=(0, 10))
        self.additional_colors_editor = SearchableListEditor(
            metadata_tab,
            label="cell_annotations",
            height=6,
            help_text="Additional cell/sample annotation obs columns selectable in the KaroSpace viewer.",
            palette=self._app_palette,
        )
        self.additional_colors_editor.grid(row=3, column=0, columnspan=3, sticky="ew", pady=(0, 10))

        self.section_metadata_editor = SearchableListEditor(
            metadata_tab,
            label="section_metadata (filter chips)",
            height=5,
            help_text="Section-level obs fields shown in the visual params bar/filter chips.",
            palette=self._app_palette,
        )
        self.section_metadata_editor.grid(row=4, column=0, columnspan=3, sticky="ew", pady=(0, 10))

        self.section_metadata_extra_editor = SearchableListEditor(
            metadata_tab,
            label="section_metadata_extra (stored only)",
            height=5,
            help_text="Section-level obs fields stored in the payload without filter chips.",
            palette=self._app_palette,
        )
        self.section_metadata_extra_editor.grid(row=5, column=0, columnspan=3, sticky="ew", pady=(0, 12))

        features_tab.columnconfigure(0, weight=1)
        self._section_label(features_tab, "FEATURE CONTENT AND STORAGE").grid(row=0, column=0, sticky="w", pady=(0, 6))
        features_intro = self._subheader_label(
            features_tab,
            "Select expression features and configure embedded or sidecar feature storage.",
        )
        features_intro.configure(wraplength=520, justify="left")
        features_intro.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self._divider(features_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))

        modality_row = self._make_sub_frame(features_tab)
        modality_row.grid(row=3, column=0, sticky="ew", pady=(0, 10))
        modality_row.columnconfigure(1, weight=1)
        self._body_label(modality_row, "Modality").grid(row=0, column=0, sticky="w", padx=(0, 8))
        self.feature_modality_combo = ctk.CTkComboBox(
            modality_row,
            variable=self.feature_modality_var,
            values=["rna"],
            state="readonly",
            command=lambda _choice: self._on_feature_modality_changed(),
            **self._theme["combo"],
        )
        self._register_combo_widget(self.feature_modality_combo)
        self.feature_modality_combo.grid(row=0, column=1, sticky="ew")

        genes_card = ctk.CTkFrame(features_tab, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", genes_card)
        genes_card.grid(row=4, column=0, sticky="ew")
        genes_card.columnconfigure(0, weight=1)
        genes_inner = self._make_sub_frame(genes_card)
        genes_inner.grid(row=0, column=0, sticky="ew", padx=12, pady=12)
        genes_inner.columnconfigure(0, weight=1)
        self._field_label(genes_inner, "Feature Selection").grid(row=0, column=0, sticky="w", pady=(0, 6))
        feature_selection_hint = self._subheader_label(
            genes_inner,
            "Selections are saved separately for each Modality and exported together.",
        )
        feature_selection_hint.configure(wraplength=500, justify="left")
        feature_selection_hint.grid(row=1, column=0, sticky="ew", pady=(0, 8))

        features_list_row = self._make_sub_frame(genes_inner)
        features_list_row.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        features_list_row.columnconfigure(0, weight=1)
        features_list_entry = ctk.CTkEntry(
            features_list_row,
            textvariable=self.features_list_var,
            placeholder_text="Text file with one feature per line",
            **self._theme["entry"],
        )
        self._register_entry_widget(features_list_entry)
        features_list_entry.grid(row=0, column=0, sticky="ew")
        features_list_button = self._secondary_button(
            features_list_row,
            "Features list",
            lambda: self._choose_file(self.features_list_var, optional=True),
            width=112,
        )
        features_list_button.grid(row=0, column=1, sticky="e", padx=(8, 0))

        self.manual_genes_editor = SearchableListEditor(
            genes_inner,
            label="features",
            height=8,
            palette=self._app_palette,
            on_change=self._update_export_estimate,
            stack_controls=True,
        )
        self.manual_genes_editor.grid(row=3, column=0, sticky="ew", pady=(0, 12))

        storage_card = ctk.CTkFrame(features_tab, **self._theme["card_frame"])
        self._register_theme_widget("card_frame", storage_card)
        storage_card.grid(row=5, column=0, sticky="ew", pady=(12, 0))
        storage_card.columnconfigure(0, weight=1)
        storage_inner = self._make_sub_frame(storage_card)
        storage_inner.grid(row=0, column=0, sticky="ew", padx=12, pady=12)
        storage_inner.columnconfigure(0, weight=1)
        self._field_label(storage_inner, "Feature Storage").grid(row=0, column=0, sticky="w", pady=(0, 6))
        storage_hint = self._subheader_label(
            storage_inner,
            "Storage settings apply to all selected modalities.",
        )
        storage_hint.configure(wraplength=500, justify="left")
        storage_hint.grid(row=1, column=0, sticky="ew", pady=(0, 8))

        feature_storage_row = self._make_sub_frame(storage_inner)
        feature_storage_row.grid(row=2, column=0, sticky="ew", pady=(0, 8))
        for column in range(3):
            feature_storage_row.columnconfigure(column, weight=1, uniform="feature_storage")

        storage_cell = self._make_sub_frame(feature_storage_row)
        storage_cell.grid(row=0, column=0, sticky="ew", padx=(0, 8))
        storage_cell.columnconfigure(0, weight=1)
        self._body_label(storage_cell, "Storage").grid(row=0, column=0, sticky="w", pady=(0, 3))
        self.feature_storage_combo = ctk.CTkComboBox(
            storage_cell,
            variable=self.feature_storage_var,
            values=["embedded", "sidecar"],
            width=120,
            state="readonly",
            command=lambda _choice: self._on_feature_storage_changed(),
            **self._theme["combo"],
        )
        self._register_combo_widget(self.feature_storage_combo)
        self.feature_storage_combo.grid(row=1, column=0, sticky="ew")

        encoding_cell = self._make_sub_frame(feature_storage_row)
        encoding_cell.grid(row=0, column=1, sticky="ew", padx=(0, 8))
        encoding_cell.columnconfigure(0, weight=1)
        self._body_label(encoding_cell, "Encoding").grid(row=0, column=0, sticky="w", pady=(0, 3))
        self.feature_encoding_combo = ctk.CTkComboBox(
            encoding_cell,
            variable=self.feature_encoding_var,
            values=["auto", "dense", "sparse"],
            width=110,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.feature_encoding_combo)
        self.feature_encoding_combo.grid(row=1, column=0, sticky="ew")

        value_cell = self._make_sub_frame(feature_storage_row)
        value_cell.grid(row=0, column=2, sticky="ew")
        value_cell.columnconfigure(0, weight=1)
        self._body_label(value_cell, "Value").grid(row=0, column=0, sticky="w", pady=(0, 3))
        self.feature_value_encoding_combo = ctk.CTkComboBox(
            value_cell,
            variable=self.feature_value_encoding_var,
            values=["uint16", "uint8"],
            width=96,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.feature_value_encoding_combo)
        self.feature_value_encoding_combo.grid(row=1, column=0, sticky="ew")

        feature_options_row = self._make_sub_frame(storage_inner)
        feature_options_row.grid(row=3, column=0, sticky="ew", pady=(0, 8))
        for column in range(4):
            feature_options_row.columnconfigure(column, weight=1, uniform="feature_options")
        self._body_label(feature_options_row, "Shard").grid(row=0, column=0, sticky="w", padx=(0, 4))
        self.feature_sidecar_shard_entry = ctk.CTkEntry(
            feature_options_row,
            textvariable=self.feature_sidecar_shard_size_var,
            width=74,
            **self._theme["entry"],
        )
        self._register_entry_widget(self.feature_sidecar_shard_entry)
        self.feature_sidecar_shard_entry.grid(row=0, column=1, sticky="ew", padx=(0, 12))
        self._body_label(feature_options_row, "Sparse >=").grid(row=0, column=2, sticky="w", padx=(0, 4))
        self.feature_sparse_zero_threshold_entry = ctk.CTkEntry(
            feature_options_row,
            textvariable=self.feature_sparse_zero_threshold_var,
            width=74,
            **self._theme["entry"],
        )
        self._register_entry_widget(self.feature_sparse_zero_threshold_entry)
        self.feature_sparse_zero_threshold_entry.grid(row=0, column=3, sticky="ew")

        karospace_row = self._make_sub_frame(storage_inner)
        karospace_row.grid(row=4, column=0, sticky="ew", pady=(4, 0))
        karospace_row.columnconfigure(0, weight=1)
        self.also_karospace_check = ctk.CTkCheckBox(
            karospace_row,
            text="Also create .karospace package (bundles all features)",
            variable=self.also_karospace_var,
            command=self._on_also_karospace_toggled,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.also_karospace_check)
        self.also_karospace_check.grid(row=0, column=0, sticky="w")
        self.also_karospace_hint = self._subheader_label(
            karospace_row,
            "Requires sidecar storage. Produces a portable .karospace next to the HTML.",
        )
        self.also_karospace_hint.configure(wraplength=500, justify="left")
        self.also_karospace_hint.grid(row=1, column=0, sticky="ew", pady=(2, 0))

        connections_tab.columnconfigure(0, weight=1)
        self._section_label(connections_tab, "CONNECTIONS").grid(row=0, column=0, sticky="w", pady=(0, 6))
        connections_intro = self._subheader_label(
            connections_tab,
            "Configure optional KaroSpace discovery summaries.",
        )
        connections_intro.configure(wraplength=520, justify="left")
        connections_intro.grid(row=1, column=0, sticky="ew", pady=(0, 8))
        self._divider(connections_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        discovery_inner = self._make_sub_frame(connections_tab)
        discovery_inner.grid(row=3, column=0, sticky="ew")
        discovery_inner.columnconfigure(0, weight=1)
        self._field_label(discovery_inner, "Discovery panels").grid(row=0, column=0, sticky="w", pady=(0, 6))
        discovery_grid = self._make_sub_frame(discovery_inner)
        discovery_grid.grid(row=1, column=0, sticky="ew")
        for column in range(2):
            discovery_grid.columnconfigure(column, weight=1, uniform="discovery")
        for column, (label, variable) in enumerate(
            (
                ("Feature correlations", self.feature_correlation_top_n_var),
                ("Spatial variable", self.spatial_variable_features_n_var),
            )
        ):
            cell = self._make_sub_frame(discovery_grid)
            cell.grid(row=0, column=column, sticky="ew", padx=(0, 8 if column < 2 else 0))
            cell.columnconfigure(0, weight=1)
            self._body_label(cell, label).grid(row=0, column=0, sticky="w", pady=(0, 3))
            entry = ctk.CTkEntry(cell, textvariable=variable, width=74, **self._theme["entry"])
            self._register_entry_widget(entry)
            entry.grid(row=1, column=0, sticky="ew")
        discovery_hint = self._subheader_label(
            discovery_inner,
            "Set numeric values to 0 to disable optional feature discovery summaries.",
        )
        discovery_hint.configure(wraplength=500, justify="left")
        discovery_hint.grid(row=2, column=0, sticky="ew", pady=(6, 0))

        statistics_tab.columnconfigure(0, weight=1)
        self._section_label(statistics_tab, "STATISTICS AND DIFFERENTIAL FEATURES").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._subheader_label(
            statistics_tab,
            "Configure Wilcoxon marker statistics, optional pseudobulk DE, and downstream pathway enrichment.",
        ).grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._divider(statistics_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        self.statistics_content = self._make_sub_frame(statistics_tab)
        self.statistics_content.grid(row=3, column=0, sticky="ew")
        self.statistics_content.columnconfigure(1, weight=1)
        self._wilcoxon_detail_rows: list[list[tk.Widget]] = []
        self._pseudobulk_detail_rows: list[list[tk.Widget]] = []
        self._pathway_detail_rows: list[list[tk.Widget]] = []

        stats_row = 0
        self.statistics_annotations_editor = SearchableListEditor(
            self.statistics_content,
            label="statistics_additional_annotations",
            height=6,
            help_text="Additional annotation columns analyzed with Wilcoxon statistics, pseudobulk DE, and interaction markers.",
            palette=self._app_palette,
        )
        self.statistics_annotations_editor.grid(row=stats_row, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        stats_row += 1

        self.statistics_modalities_editor = SearchableListEditor(
            self.statistics_content,
            label="statistics_modalities",
            height=4,
            help_text="Modalities used for Wilcoxon, pseudobulk, and interaction statistics. Empty uses the dataset default modality.",
            palette=self._app_palette,
        )
        self.statistics_modalities_editor.grid(row=stats_row, column=0, columnspan=3, sticky="ew", pady=(0, 12))
        stats_row += 1

        statistics_distribution_row = self._make_sub_frame(self.statistics_content)
        for column in range(4):
            statistics_distribution_row.columnconfigure(column, weight=1, uniform="statistics_distribution")
        for column, (label, variable, values) in enumerate(
            (
                ("Counts layer", self.statistics_counts_layer_var, ["none"]),
                ("Normalization", self.statistics_normalization_var, ["RC", "LogNormalize"]),
                ("Scale factor", self.statistics_scale_factor_var, None),
                ("Normalized layer", self.statistics_normalized_layer_var, ["off"]),
            )
        ):
            cell = self._make_sub_frame(statistics_distribution_row)
            cell.grid(row=0, column=column, sticky="ew", padx=(0, 8 if column < 3 else 0))
            cell.columnconfigure(0, weight=1)
            self._body_label(cell, label).grid(row=0, column=0, sticky="w", pady=(0, 3))
            if values is None:
                entry = ctk.CTkEntry(cell, textvariable=variable, width=100, **self._theme["entry"])
                self._register_entry_widget(entry)
                entry.grid(row=1, column=0, sticky="ew")
            else:
                combo = ctk.CTkComboBox(cell, variable=variable, values=values, state="readonly", **self._theme["combo"])
                self._register_combo_widget(combo)
                combo.grid(row=1, column=0, sticky="ew")
                if variable is self.statistics_counts_layer_var:
                    self.statistics_counts_layer_combo = combo
                elif variable is self.statistics_normalized_layer_var:
                    self.statistics_normalized_layer_combo = combo
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Distribution values",
            widget=statistics_distribution_row,
            hint="Controls Statistics feature distributions. A normalized layer overrides counts, normalization, and scale factor.",
        )

        statistics_filter_row = self._make_sub_frame(self.statistics_content)
        for label, variable, width in (
            ("Min cell counts", self.statistics_min_cell_counts_var, 76),
            ("Min feature counts", self.statistics_min_feature_counts_var, 76),
            ("CPUs", self.statistics_n_cpus_var, 58),
        ):
            self._body_label(statistics_filter_row, label).pack(side="left", padx=(0, 4))
            entry = ctk.CTkEntry(statistics_filter_row, textvariable=variable, width=width, **self._theme["entry"])
            self._register_entry_widget(entry)
            entry.pack(side="left", padx=(0, 10))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Statistics filters",
            widget=statistics_filter_row,
            hint="Applies to Wilcoxon and pseudobulk statistics before analysis.",
        )

        statistics_categories_row = self._make_sub_frame(self.statistics_content)
        statistics_categories_entry = self._entry(
            statistics_categories_row,
            self.statistics_contrast_categories_var,
            placeholder='{"cell_type":["Astrocyte","B cell"]}',
        )
        statistics_categories_entry.pack(side="left", fill="x", expand=True)
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Contrast categories",
            widget=statistics_categories_row,
        )

        wilcoxon_enable_row = self._make_sub_frame(self.statistics_content)
        self.wilcoxon_enabled_check = ctk.CTkCheckBox(
            wilcoxon_enable_row,
            text="Run Wilcoxon statistics",
            variable=self.wilcoxon_enabled_var,
            command=self._sync_analysis_controls,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.wilcoxon_enabled_check)
        self.wilcoxon_enabled_check.pack(side="left")
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Wilcoxon",
            widget=wilcoxon_enable_row,
        )

        wilcoxon_mode_row = self._make_sub_frame(self.statistics_content)
        self._body_label(wilcoxon_mode_row, "Mode").pack(side="left")
        self.wilcoxon_mode_combo = ctk.CTkComboBox(
            wilcoxon_mode_row,
            variable=self.wilcoxon_mode_var,
            values=["auto", "force"],
            width=100,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.wilcoxon_mode_combo)
        self.wilcoxon_mode_combo.pack(side="left", padx=(6, 12))
        self._body_label(wilcoxon_mode_row, "Runtime limit").pack(side="left")
        wilcoxon_runtime_entry = ctk.CTkEntry(
            wilcoxon_mode_row,
            textvariable=self.wilcoxon_runtime_limit_var,
            width=100,
            **self._theme["entry"],
        )
        self._register_entry_widget(wilcoxon_runtime_entry)
        wilcoxon_runtime_entry.pack(side="left", padx=(6, 0))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Wilcoxon mode",
            widget=wilcoxon_mode_row,
            row_collector=self._wilcoxon_detail_rows,
        )

        wilcoxon_threshold_row = self._make_sub_frame(self.statistics_content)
        for label, variable, width in (
            ("Min cells/group", self.wilcoxon_min_cells_per_group_var, 76),
            ("Min pct", self.wilcoxon_min_pct_expressed_var, 62),
            ("Padj", self.wilcoxon_padj_cutoff_var, 62),
            ("Log2FC", self.wilcoxon_log2fc_cutoff_var, 62),
        ):
            self._body_label(wilcoxon_threshold_row, label).pack(side="left", padx=(0, 4))
            entry = ctk.CTkEntry(wilcoxon_threshold_row, textvariable=variable, width=width, **self._theme["entry"])
            self._register_entry_widget(entry)
            entry.pack(side="left", padx=(0, 10))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Wilcoxon thresholds",
            widget=wilcoxon_threshold_row,
            row_collector=self._wilcoxon_detail_rows,
        )

        wilcoxon_display_row = self._make_sub_frame(self.statistics_content)
        self._body_label(wilcoxon_display_row, "P adjust").pack(side="left")
        self.wilcoxon_p_adjust_combo = ctk.CTkComboBox(
            wilcoxon_display_row,
            variable=self.wilcoxon_p_adjust_method_var,
            values=["fdr_bh", "bonferroni", "holm", "none"],
            width=120,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.wilcoxon_p_adjust_combo)
        self.wilcoxon_p_adjust_combo.pack(side="left", padx=(6, 12))
        self._body_label(wilcoxon_display_row, "Embed top").pack(side="left")
        wilcoxon_embed_entry = ctk.CTkEntry(
            wilcoxon_display_row,
            textvariable=self.wilcoxon_embed_top_n_per_comparison_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(wilcoxon_embed_entry)
        wilcoxon_embed_entry.pack(side="left", padx=(4, 10))
        self._body_label(wilcoxon_display_row, "Rows/category").pack(side="left")
        wilcoxon_rows_entry = ctk.CTkEntry(
            wilcoxon_display_row,
            textvariable=self.wilcoxon_top_n_per_category_var,
            width=80,
            **self._theme["entry"],
        )
        self._register_entry_widget(wilcoxon_rows_entry)
        wilcoxon_rows_entry.pack(side="left", padx=(4, 0))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Wilcoxon results",
            widget=wilcoxon_display_row,
            row_collector=self._wilcoxon_detail_rows,
        )

        neighborhoods_tab.columnconfigure(0, weight=1)
        self._section_label(neighborhoods_tab, "NEIGHBORHOODS AND INTERACTIONS").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._subheader_label(
            neighborhoods_tab,
            "Configure neighbor composition statistics and contact-conditioned interaction markers.",
        ).grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._divider(neighborhoods_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        self.neighborhoods_content = self._make_sub_frame(neighborhoods_tab)
        self.neighborhoods_content.grid(row=3, column=0, sticky="ew")
        self.neighborhoods_content.columnconfigure(1, weight=1)

        overlays_tab.columnconfigure(0, weight=1)
        self._section_label(overlays_tab, "SECTIONS").grid(row=0, column=0, sticky="w", pady=(0, 6))
        self._subheader_label(
            overlays_tab,
            "Configure section layout, scale, image overlays, rotations, and preview serving.",
        ).grid(row=1, column=0, sticky="w", pady=(0, 8))
        self._divider(overlays_tab, height=1).grid(row=2, column=0, sticky="ew", pady=(0, 10))
        self.overlays_content = self._make_sub_frame(overlays_tab)
        self.overlays_content.grid(row=3, column=0, sticky="ew")
        self.overlays_content.columnconfigure(1, weight=1)

        min_panel_row = self._make_sub_frame(self.overlays_content)
        min_panel_entry = ctk.CTkEntry(min_panel_row, textvariable=self.min_panel_size_var, width=90, **self._theme["entry"])
        self._register_entry_widget(min_panel_entry)
        min_panel_entry.pack(side="left")
        self._body_label(min_panel_row, "px").pack(side="left", padx=(8, 0))
        self._option_row(
            self.overlays_content,
            0,
            "Min panel size",
            widget=min_panel_row,
            hint="Minimum section panel width in exported KaroSpace HTML.",
        )

        spot_row = self._make_sub_frame(self.overlays_content)
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
            self.overlays_content,
            2,
            "Spot size",
            widget=spot_row,
            hint="Use auto/adaptive/density or a positive number.",
        )

        pseudobulk_enable_row = self._make_sub_frame(self.statistics_content)
        self.pseudobulk_enabled_check = ctk.CTkCheckBox(
            pseudobulk_enable_row,
            text="Run pseudobulk analysis",
            variable=self.pseudobulk_enabled_var,
            command=self._sync_analysis_controls,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.pseudobulk_enabled_check)
        self.pseudobulk_enabled_check.pack(side="left")
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Pseudobulk analysis",
            widget=pseudobulk_enable_row,
        )

        neighbor_row = self._make_sub_frame(self.neighborhoods_content)
        self._body_label(neighbor_row, "Annotations").pack(side="left")
        self.neighbor_stats_annotations_entry = ctk.CTkEntry(
            neighbor_row,
            textvariable=self.neighbor_stats_annotations_var,
            width=140,
            **self._theme["entry"],
        )
        self._register_entry_widget(self.neighbor_stats_annotations_entry)
        self.neighbor_stats_annotations_entry.pack(side="left", padx=(6, 12))
        self._body_label(neighbor_row, "Permutations").pack(side="left", padx=(12, 4))
        neighbor_perm_entry = ctk.CTkEntry(neighbor_row, textvariable=self.neighbor_permutations_var, width=76, **self._theme["entry"])
        self._register_entry_widget(neighbor_perm_entry)
        neighbor_perm_entry.pack(side="left")
        self._body_label(neighbor_row, "Seed").pack(side="left", padx=(12, 4))
        neighbor_seed_entry = ctk.CTkEntry(neighbor_row, textvariable=self.neighbor_stats_seed_var, width=76, **self._theme["entry"])
        self._register_entry_widget(neighbor_seed_entry)
        neighbor_seed_entry.pack(side="left")
        self._option_row(
            self.neighborhoods_content,
            0,
            "Neighbor stats",
            widget=neighbor_row,
            hint="--neighbor-stats-annotations accepts auto, blank to disable, or comma-separated obs columns.",
        )

        pseudobulk_row_1 = self._make_sub_frame(self.statistics_content)
        self.pseudobulk_replicate_combo = ctk.CTkComboBox(
            pseudobulk_row_1,
            variable=self.pseudobulk_replicate_annotation_var,
            values=[""],
            width=220,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.pseudobulk_replicate_combo)
        self.pseudobulk_replicate_combo.pack(side="left", fill="x", expand=True)
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Replicate annotation",
            widget=pseudobulk_row_1,
            row_collector=self._pseudobulk_detail_rows,
        )

        pseudobulk_row_2 = self._make_sub_frame(self.statistics_content)
        for label, variable, width in (
            ("Min cells/sample", self.pseudobulk_min_cells_per_pseudobulk_var, 70),
            ("Min reps", self.pseudobulk_min_replicates_var, 58),
            ("Min pct", self.pseudobulk_min_pct_expressed_var, 62),
            ("Padj", self.pseudobulk_padj_cutoff_var, 62),
            ("Log2FC", self.pseudobulk_log2fc_cutoff_var, 62),
        ):
            self._body_label(pseudobulk_row_2, label).pack(side="left", padx=(0, 4))
            entry = ctk.CTkEntry(pseudobulk_row_2, textvariable=variable, width=width, **self._theme["entry"])
            self._register_entry_widget(entry)
            entry.pack(side="left", padx=(0, 8))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Pseudobulk thresholds",
            widget=pseudobulk_row_2,
            row_collector=self._pseudobulk_detail_rows,
        )

        pseudobulk_row_3 = self._make_sub_frame(self.statistics_content)
        self._body_label(pseudobulk_row_3, "P adjust").pack(side="left")
        self.pseudobulk_p_adjust_combo = ctk.CTkComboBox(
            pseudobulk_row_3,
            variable=self.pseudobulk_p_adjust_method_var,
            values=["fdr_bh", "bonferroni", "holm", "none"],
            width=120,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.pseudobulk_p_adjust_combo)
        self.pseudobulk_p_adjust_combo.pack(side="left", padx=(6, 12))
        self._body_label(pseudobulk_row_3, "Fit").pack(side="left")
        self.pseudobulk_fit_type_combo = ctk.CTkComboBox(
            pseudobulk_row_3,
            variable=self.pseudobulk_deseq2_fit_type_var,
            values=["parametric", "mean"],
            width=120,
            state="readonly",
            **self._theme["combo"],
        )
        self._register_combo_widget(self.pseudobulk_fit_type_combo)
        self.pseudobulk_fit_type_combo.pack(side="left", padx=(6, 12))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Pseudobulk fit",
            widget=pseudobulk_row_3,
            row_collector=self._pseudobulk_detail_rows,
        )

        marker_row = self._make_sub_frame(self.statistics_content)
        marker_top_n_entry = ctk.CTkEntry(
            marker_row,
            textvariable=self.pseudobulk_embed_top_n_per_comparison_var,
            width=90,
            **self._theme["entry"],
        )
        self._register_entry_widget(marker_top_n_entry)
        marker_top_n_entry.pack(side="left")
        self._body_label(marker_row, "significant DE features per comparison").pack(side="left", padx=(8, 0))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Auto-embedded DE genes",
            widget=marker_row,
            hint="Maps to pseudobulk_embed_top_n_per_comparison in embedded feature mode.",
            row_collector=self._pseudobulk_detail_rows,
        )

        pathway_enable_row = self._make_sub_frame(self.statistics_content)
        self.pathway_enabled_check = ctk.CTkCheckBox(
            pathway_enable_row,
            text="Run pathway enrichment",
            variable=self.pathway_enabled_var,
            command=self._sync_analysis_controls,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.pathway_enabled_check)
        self.pathway_enabled_check.pack(side="left")
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Pathway enrichment",
            widget=pathway_enable_row,
        )

        pathway_row = self._make_sub_frame(self.statistics_content)
        self._body_label(pathway_row, "Organism").pack(side="left")
        pathway_organism_entry = ctk.CTkEntry(pathway_row, textvariable=self.pathway_organism_var, width=100, **self._theme["entry"])
        self._register_entry_widget(pathway_organism_entry)
        pathway_organism_entry.pack(side="left", padx=(4, 10))
        self._body_label(pathway_row, "Top N").pack(side="left")
        pathway_top_entry = ctk.CTkEntry(pathway_row, textvariable=self.pathway_top_n_var, width=58, **self._theme["entry"])
        self._register_entry_widget(pathway_top_entry)
        pathway_top_entry.pack(side="left", padx=(4, 10))
        self._body_label(pathway_row, "Min overlap").pack(side="left")
        pathway_overlap_entry = ctk.CTkEntry(pathway_row, textvariable=self.pathway_min_overlap_var, width=58, **self._theme["entry"])
        self._register_entry_widget(pathway_overlap_entry)
        pathway_overlap_entry.pack(side="left", padx=(4, 10))
        self._body_label(pathway_row, "GSEA perms").pack(side="left")
        pathway_perm_entry = ctk.CTkEntry(pathway_row, textvariable=self.pathway_gsea_permutations_var, width=70, **self._theme["entry"])
        self._register_entry_widget(pathway_perm_entry)
        pathway_perm_entry.pack(side="left", padx=(4, 0))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Pathways",
            widget=pathway_row,
            hint="GMT files can be supplied below as comma-separated paths; blank uses KaroSpace defaults.",
            row_collector=self._pathway_detail_rows,
        )

        pathway_gmt_row = self._make_sub_frame(self.statistics_content)
        pathway_gmt_entry = ctk.CTkEntry(pathway_gmt_row, textvariable=self.pathway_gmt_var, **self._theme["entry"])
        self._register_entry_widget(pathway_gmt_entry)
        pathway_gmt_entry.pack(side="left", fill="x", expand=True)
        pathway_gmt_button = self._secondary_button(
            pathway_gmt_row,
            "Browse",
            self._choose_pathway_gmt_files,
            width=96,
        )
        self._secondary_action_buttons.append(pathway_gmt_button)
        self._style_secondary_action_button(pathway_gmt_button)
        pathway_gmt_button.pack(side="left", padx=(8, 0))
        stats_row = self._option_row(
            self.statistics_content,
            stats_row,
            "Pathway GMT",
            widget=pathway_gmt_row,
            row_collector=self._pathway_detail_rows,
        )

        self._interaction_detail_rows: list[list[tk.Widget]] = []
        interaction_row_1 = self._make_sub_frame(self.neighborhoods_content)
        self.interaction_markers_enabled_check = ctk.CTkCheckBox(
            interaction_row_1,
            text="Run contact-conditioned interaction markers",
            variable=self.interaction_markers_enabled_var,
            command=self._sync_analysis_controls,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.interaction_markers_enabled_check)
        self.interaction_markers_enabled_check.pack(side="left")
        self._option_row(
            self.neighborhoods_content,
            2,
            "Interaction markers",
            widget=interaction_row_1,
            hint="Uses the main annotation plus statistics_additional_annotations.",
        )

        interaction_row_2 = self._make_sub_frame(self.neighborhoods_content)
        self._body_label(interaction_row_2, "Top targets").pack(side="left")
        interaction_top_targets_entry = ctk.CTkEntry(
            interaction_row_2,
            textvariable=self.interaction_markers_top_targets_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(interaction_top_targets_entry)
        interaction_top_targets_entry.pack(side="left", padx=(4, 10))
        self._body_label(interaction_row_2, "Top features").pack(side="left")
        interaction_top_features_entry = ctk.CTkEntry(
            interaction_row_2,
            textvariable=self.interaction_markers_top_features_var,
            width=70,
            **self._theme["entry"],
        )
        self._register_entry_widget(interaction_top_features_entry)
        interaction_top_features_entry.pack(side="left", padx=(4, 10))
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
            self.neighborhoods_content,
            3,
            "Interaction limits",
            widget=interaction_row_2,
            hint="Maps to interaction_markers_top_targets/top_features/min_cells/min_neighbors.",
            row_collector=self._interaction_detail_rows,
        )

        metadata_json_row = self._make_sub_frame(metadata_tab)
        metadata_json_row.columnconfigure(0, weight=1)
        metadata_value_order_entry = self._entry(
            metadata_json_row,
            self.metadata_value_order_var,
            placeholder='{"condition":["control","treated"]}',
        )
        metadata_value_order_entry.pack(side="left", fill="x", expand=True)
        self._option_row(
            metadata_tab,
            8,
            "Metadata value order JSON",
            widget=metadata_json_row,
        )

        metadata_labels_row = self._make_sub_frame(metadata_tab)
        metadata_labels_entry = self._entry(
            metadata_labels_row,
            self.metadata_labels_var,
            placeholder='{"sample_id":"Sample"}',
        )
        metadata_labels_entry.pack(side="left", fill="x", expand=True)
        self._option_row(metadata_tab, 9, "Metadata labels JSON", widget=metadata_labels_row)

        metadata_max_row = self._make_sub_frame(metadata_tab)
        metadata_max_entry = ctk.CTkEntry(metadata_max_row, textvariable=self.metadata_max_columns_var, width=90, **self._theme["entry"])
        self._register_entry_widget(metadata_max_entry)
        metadata_max_entry.pack(side="left")
        self._option_row(
            metadata_tab,
            10,
            "Metadata max columns",
            widget=metadata_max_row,
        )

        scalebar_row = self._make_sub_frame(self.overlays_content)
        scalebar_entry = ctk.CTkEntry(scalebar_row, textvariable=self.scalebar_unit_var, width=90, **self._theme["entry"])
        self._register_entry_widget(scalebar_entry)
        scalebar_entry.pack(side="left")
        self._option_row(self.overlays_content, 4, "Scalebar unit", widget=scalebar_row)

        viewer_info_direct_row = self._make_sub_frame(viewer_tab)
        viewer_info_direct_entry = ctk.CTkEntry(
            viewer_info_direct_row,
            textvariable=self.viewer_info_html_var,
            placeholder_text="<p>Viewer information</p>",
            **self._theme["entry"],
        )
        self._register_entry_widget(viewer_info_direct_entry)
        viewer_info_direct_entry.pack(side="left", fill="x", expand=True)
        self._option_row(
            viewer_tab,
            9,
            "Viewer info HTML",
            widget=viewer_info_direct_row,
            hint="Maps to --viewer-info-html. A file below overrides this value, matching the CLI.",
        )

        overlay_json_row = self._make_sub_frame(self.overlays_content)
        self._body_label(overlay_json_row, "Rotations").pack(side="left")
        section_rotations_entry = ctk.CTkEntry(overlay_json_row, textvariable=self.section_rotations_var, width=180, **self._theme["entry"])
        self._register_entry_widget(section_rotations_entry)
        section_rotations_entry.pack(side="left", padx=(4, 10))
        self._body_label(overlay_json_row, "Images JSON").pack(side="left")
        section_images_entry = ctk.CTkEntry(overlay_json_row, textvariable=self.section_images_var, width=180, **self._theme["entry"])
        self._register_entry_widget(section_images_entry)
        section_images_entry.pack(side="left", padx=(4, 0))
        self._option_row(
            self.overlays_content,
            5,
            "Section overlays",
            widget=overlay_json_row,
            hint="Rotations accept section:angle CSV or JSON. Images expects the KaroSpace section_images JSON object.",
        )

        image_size_row = self._make_sub_frame(self.overlays_content)
        self._body_label(image_size_row, "Images max px").pack(side="left")
        section_images_max_entry = ctk.CTkEntry(image_size_row, textvariable=self.section_images_max_px_var, width=80, **self._theme["entry"])
        self._register_entry_widget(section_images_max_entry)
        section_images_max_entry.pack(side="left", padx=(4, 0))
        self._option_row(self.overlays_content, 7, "Section image size", widget=image_size_row)

        deconv_row = self._make_sub_frame(self.overlays_content)
        deconv_entry = ctk.CTkEntry(deconv_row, textvariable=self.deconvolutions_var, width=220, **self._theme["entry"])
        self._register_entry_widget(deconv_entry)
        deconv_entry.pack(side="left", fill="x", expand=True)
        self._option_row(
            self.overlays_content,
            9,
            "Deconvolutions JSON",
            widget=deconv_row,
        )

        viewer_info_row = self._make_sub_frame(viewer_tab)
        viewer_info_entry = ctk.CTkEntry(viewer_info_row, textvariable=self.viewer_info_html_file_var, **self._theme["entry"])
        self._register_entry_widget(viewer_info_entry)
        viewer_info_entry.pack(side="left", fill="x", expand=True)
        viewer_info_button = self._secondary_button(
            viewer_info_row,
            "Info HTML",
            lambda: self._choose_file(self.viewer_info_html_file_var, optional=True),
            width=96,
        )
        viewer_info_button.pack(side="left", padx=(8, 0))
        self._option_row(
            viewer_tab,
            13,
            "Viewer info HTML file",
            widget=viewer_info_row,
        )

        tutorial_row = self._make_sub_frame(viewer_tab)
        self.tutorial_check = ctk.CTkCheckBox(tutorial_row, text="Enable tutorial", variable=self.tutorial_var, **self._theme["checkbox"])
        self._register_checkbox_widget(self.tutorial_check)
        self.tutorial_check.pack(side="left")
        self._option_row(viewer_tab, 15, "Tutorial", widget=tutorial_row)

        reproducibility_row = self._make_sub_frame(viewer_tab)
        self.reproducibility_check = ctk.CTkCheckBox(
            reproducibility_row,
            text="Embed export arguments and resolved settings",
            variable=self.embed_reproducibility_info_var,
            **self._theme["checkbox"],
        )
        self._register_checkbox_widget(self.reproducibility_check)
        self.reproducibility_check.pack(side="left")
        self._option_row(
            viewer_tab,
            16,
            "Reproducibility info",
            widget=reproducibility_row,
            hint="Uncheck to match CLI --no-reproducibility-info.",
        )

        serve_row = self._make_sub_frame(self.overlays_content)
        self.serve_check = ctk.CTkCheckBox(serve_row, text="Serve after export", variable=self.serve_var, **self._theme["checkbox"])
        self._register_checkbox_widget(self.serve_check)
        self.serve_check.pack(side="left")
        self._body_label(serve_row, "Port").pack(side="left", padx=(12, 4))
        serve_port_entry = ctk.CTkEntry(serve_row, textvariable=self.port_var, width=90, **self._theme["entry"])
        self._register_entry_widget(serve_port_entry)
        serve_port_entry.pack(side="left")
        self._option_row(
            self.overlays_content,
            11,
            "Preview server",
            widget=serve_row,
            hint="Optional local server to open the configured output.",
        )
        button_row = self._make_sub_frame(controls_outer)
        button_row.grid(row=2, column=0, sticky="ew", padx=18, pady=(16, 22))

        self.inspect_btn = self._secondary_button(button_row, "Inspect Dataset", self._inspect_h5ad, width=145)
        self._style_secondary_action_button(self.inspect_btn)
        self.inspect_btn.pack(side="left")

        self.inspect_loading_label = ctk.CTkLabel(button_row, text="⟳ Inspecting...", **self._theme["section_label"])
        self._register_theme_widget("section_label", self.inspect_loading_label)

        self.export_btn = self._secondary_button(button_row, "Build Viewer", self._on_export_button, width=140)
        self.export_btn.pack(side="right")
        self._inspection_gated_widgets.append(self.export_btn)

        runtime_top = self._make_sub_frame(side)
        runtime_top.grid(row=0, column=0, sticky="ew")
        runtime_top.columnconfigure(0, weight=1)
        self.runtime_title_label = self._section_label(runtime_top, "RUNTIME")
        self.runtime_title_label.grid(row=0, column=0, sticky="w")
        self.runtime_chip_label = self._pill_label(runtime_top, text="READY", muted=True)
        self.runtime_chip_label.grid(row=0, column=1, sticky="e")
        self._right_panel_widgets.append(self.runtime_chip_label)

        self.progress = ctk.CTkProgressBar(side, mode="determinate", **self._theme["progress"])
        self._register_theme_widget("progress", self.progress)
        self.progress.grid(row=1, column=0, sticky="ew", pady=(8, 0))
        self.progress.set(0.0)
        self._right_panel_widgets.append(self.progress)

        launch_row = self._make_sub_frame(side)
        launch_row.grid(row=2, column=0, sticky="ew", pady=(14, 12))
        self.open_output_btn = self._secondary_button(launch_row, "Open Output Folder", self._open_output_folder, width=165)
        self.open_output_btn.pack(side="left")
        self.open_viewer_btn = self._secondary_button(launch_row, "Open Viewer", self._open_viewer, width=125)
        self.open_viewer_btn.pack(side="left", padx=(10, 0))
        self._right_panel_widgets.extend([self.open_output_btn, self.open_viewer_btn])

        log_wrap = self._make_sub_frame(side)
        log_wrap.grid(row=3, column=0, sticky="nsew")
        log_wrap.columnconfigure(0, weight=1)
        log_wrap.rowconfigure(1, weight=1)

        self.event_log_label = self._section_label(log_wrap, "EVENT LOG")
        self.event_log_label.grid(row=0, column=0, sticky="w", pady=(0, 6))
        self.log_text = ctk.CTkTextbox(log_wrap, wrap="word", **self._theme["textbox"])
        self.log_text.grid(row=1, column=0, sticky="nsew")
        self.log_text.configure(state="disabled")

        self.downsample_var.trace_add("write", lambda *_: self._update_export_estimate())
        self.section_groupby_var.trace_add("write", lambda *_: self._on_section_groupby_changed())
        self.initial_color_var.trace_add("write", lambda *_: self._sync_required_cell_annotations())
        self.outline_by_var.trace_add("write", lambda *_: self._sync_required_section_metadata())
        self.wilcoxon_enabled_var.trace_add("write", lambda *_: self._sync_analysis_controls())
        self.pseudobulk_enabled_var.trace_add("write", lambda *_: self._sync_analysis_controls())
        self.pathway_enabled_var.trace_add("write", lambda *_: self._sync_analysis_controls())
        self.interaction_markers_enabled_var.trace_add("write", lambda *_: self._sync_analysis_controls())
        self.status_var.trace_add("write", lambda *_: self._sync_runtime_chip())
        self.h5ad_var.trace_add("write", lambda *_: self._refresh_input_gate())
        self.theme_var.trace_add("write", lambda *_: self._sync_theme_toggle())
        self._apply_preset("default", log=False)
        self._sync_app_theme_to_viewer_setting()
        self._sync_runtime_chip()
        self._sync_analysis_controls()
        self._refresh_input_gate()
        self._update_export_estimate()

    def _path_field(
        self,
        parent: tk.Widget,
        row: int,
        label: str,
        variable: tk.StringVar,
        choose_file: bool,
        optional: bool = False,
        gated: bool = False,
    ) -> int:
        row_widgets: list[tk.Widget] = []
        label_widget = self._field_label(parent, label)
        label_widget.grid(row=row, column=0, sticky="nw", pady=(0, 8), padx=(0, 14))
        row_widgets.append(label_widget)
        placeholder = "/absolute/path/to/input.h5ad or input.zarr" if variable is self.h5ad_var else None
        entry = self._entry(parent, variable, placeholder=placeholder)
        entry.grid(row=row, column=1, sticky="ew", pady=(0, 8))
        row_widgets.append(entry)

        if choose_file and variable is self.h5ad_var:
            button_wrap = self._make_sub_frame(parent)
            file_button = self._secondary_button(
                button_wrap,
                "H5AD",
                lambda: self._choose_file(variable, optional=optional),
                width=74,
            )
            self._secondary_action_buttons.append(file_button)
            self._style_secondary_action_button(file_button)
            file_button.pack(side="left")
            zarr_button = self._secondary_button(
                button_wrap,
                "Zarr Dir",
                self._choose_input_directory,
                width=82,
            )
            self._secondary_action_buttons.append(zarr_button)
            self._style_secondary_action_button(zarr_button)
            zarr_button.pack(side="left", padx=(6, 0))
            button = button_wrap
        elif choose_file:
            button = self._secondary_button(parent, "Browse", lambda: self._choose_file(variable, optional=optional), width=96)
            self._secondary_action_buttons.append(button)
            self._style_secondary_action_button(button)
        else:
            button = self._secondary_button(parent, "Browse", lambda: self._choose_dir(variable), width=96)
            self._secondary_action_buttons.append(button)
            self._style_secondary_action_button(button)
        button.grid(row=row, column=2, sticky="e", pady=(0, 8), padx=(10, 0))
        row_widgets.append(button)
        if gated:
            self._inspection_gated_widgets.extend(row_widgets)
        return row + 1

    def _build_parameter_state_buttons(self, parent: tk.Widget, row: int) -> int:
        label = self._field_label(parent, "Parameter JSON")
        label.grid(row=row, column=0, sticky="nw", pady=(0, 8), padx=(0, 14))
        state_row = self._make_sub_frame(parent)
        state_row.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(0, 8))
        save_button = self._secondary_button(state_row, "Save JSON", self._save_parameter_state_json, width=104)
        import_button = self._secondary_button(state_row, "Import JSON", self._import_parameter_state_json, width=112)
        for button in (save_button, import_button):
            self._secondary_action_buttons.append(button)
            self._style_secondary_action_button(button)
        save_button.pack(side="left")
        import_button.pack(side="left", padx=(8, 0))
        self.save_parameters_btn = save_button
        self.import_parameters_btn = import_button
        return row + 1

    def _build_spatialdata_table_selector(self, parent: tk.Widget, row: int) -> None:
        label = self._field_label(parent, "SpatialData table")
        label.grid(row=row, column=0, sticky="nw", pady=(0, 8), padx=(0, 14))
        self.spatialdata_table_combo = ctk.CTkComboBox(
            parent,
            variable=self.spatialdata_table_var,
            values=[],
            state="disabled",
            command=self._on_spatialdata_table_selected,
            **self._theme["combo"],
        )
        self._register_combo_widget(self.spatialdata_table_combo)
        self.spatialdata_table_combo.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(0, 4))
        hint = self._subheader_label(
            parent,
            "Shown for SpatialData .zarr inputs. Inspect selects the first available table and loads all table choices.",
        )
        hint.grid(row=row + 1, column=1, columnspan=2, sticky="w", pady=(0, 10))
        self._spatialdata_table_widgets = [label, self.spatialdata_table_combo, hint]
        self._set_widgets_visible(self._spatialdata_table_widgets, False)

    def _option_row(
        self,
        parent: tk.Widget,
        row: int,
        label: str,
        widget: tk.Widget,
        hint: str | None = None,
        gated: bool = False,
        row_collector: list[list[tk.Widget]] | None = None,
    ) -> int:
        row_widgets: list[tk.Widget] = []
        label_widget = self._field_label(parent, label)
        label_widget.grid(row=row, column=0, sticky="nw", pady=(0, 4), padx=(0, 14))
        row_widgets.append(label_widget)
        widget.grid(row=row, column=1, columnspan=2, sticky="ew", pady=(0, 4))
        row_widgets.append(widget)
        row += 1
        if hint:
            hint_widget = self._subheader_label(parent, hint)
            hint_widget.grid(row=row, column=1, columnspan=2, sticky="w", pady=(0, 10))
            row_widgets.append(hint_widget)
            row += 1
        if gated:
            self._inspection_gated_widgets.extend(row_widgets)
        if row_collector is not None:
            row_collector.append(row_widgets)
        return row

    def _entry(self, parent: tk.Widget, variable: tk.StringVar, *, placeholder: str | None = None) -> ctk.CTkEntry:
        kwargs = dict(self._theme["entry"])
        if not placeholder:
            entry = ctk.CTkEntry(parent, textvariable=variable, **kwargs)
            self._register_entry_widget(entry)
            return entry

        entry = ctk.CTkEntry(parent, **kwargs)
        self._register_entry_widget(entry)
        placeholder_state = {"active": False, "syncing": False}

        def show_placeholder() -> None:
            if variable.get().strip():
                return
            placeholder_state["syncing"] = True
            placeholder_state["active"] = True
            entry.configure(text_color=self._app_palette["muted"])
            entry.delete(0, "end")
            entry.insert(0, placeholder)
            placeholder_state["syncing"] = False

        def show_value(value: str) -> None:
            placeholder_state["syncing"] = True
            placeholder_state["active"] = False
            entry.configure(text_color=self._app_palette["text"])
            entry.delete(0, "end")
            entry.insert(0, value)
            placeholder_state["syncing"] = False

        def sync_from_variable(*_args: object) -> None:
            if placeholder_state["syncing"]:
                return
            value = variable.get()
            if value.strip():
                show_value(value)
            elif entry.focus_get() is entry:
                placeholder_state["active"] = False
                entry.configure(text_color=self._app_palette["text"])
                entry.delete(0, "end")
            else:
                show_placeholder()

        def on_focus_in(_event: object) -> None:
            if placeholder_state["active"]:
                placeholder_state["active"] = False
                entry.configure(text_color=self._app_palette["text"])
                entry.delete(0, "end")

        def on_focus_out(_event: object) -> None:
            if not entry.get().strip():
                variable.set("")
                show_placeholder()

        def on_key_release(_event: object) -> None:
            if placeholder_state["active"] or placeholder_state["syncing"]:
                return
            placeholder_state["syncing"] = True
            variable.set(entry.get())
            placeholder_state["syncing"] = False

        variable.trace_add("write", sync_from_variable)
        entry.bind("<FocusIn>", on_focus_in)
        entry.bind("<FocusOut>", on_focus_out)
        entry.bind("<KeyRelease>", on_key_release)
        self._placeholder_refreshers.append(sync_from_variable)
        show_placeholder()
        return entry

    def _coords_dropdown(self, parent: tk.Widget) -> ctk.CTkOptionMenu:
        self.coords_menu = ctk.CTkOptionMenu(
            parent,
            variable=self.coords_var,
            values=["auto", "obsm:spatial", "obs:centroid_x_y"],
            **self._theme["combo"],
        )
        self._register_combo_widget(self.coords_menu)
        return self.coords_menu

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

    def _sync_analysis_controls(self) -> None:
        def enabled(var_name: str, default: bool) -> bool:
            var = getattr(self, var_name, None)
            if var is None:
                return default
            return bool(var.get())

        row_groups = (
            ("_wilcoxon_detail_rows", enabled("wilcoxon_enabled_var", True)),
            ("_pseudobulk_detail_rows", enabled("pseudobulk_enabled_var", False)),
            ("_pathway_detail_rows", enabled("pathway_enabled_var", False)),
            ("_interaction_detail_rows", enabled("interaction_markers_enabled_var", True)),
        )
        for attr, visible in row_groups:
            rows = getattr(self, attr, None)
            if rows is None:
                continue
            widgets = [widget for row_widgets in rows for widget in row_widgets]
            self._set_widgets_visible(widgets, visible)

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
        if not self.outdir_var.get().strip():
            self.outdir_var.set(str(Path.cwd().resolve()))
        if not self.output_html_var.get().strip():
            self.output_html_var.set("karospace.html")
        self.coords_var.set("auto")
        self.spatial_key_var.set("spatial")
        self.spatial_x_var.set("")
        self.spatial_y_var.set("")
        self.spatialdata_table_var.set("")
        self.section_order_var.set("")
        self.metadata_value_order_var.set("")
        self.metadata_max_columns_var.set("")
        self.serve_var.set(False)
        self.port_var.set("8000")
        self.downsample_var.set("")
        self.section_groupby_var.set("sample_id")
        self.initial_color_var.set("leiden")
        self.title_var.set("KaroSpace")
        if not self.theme_var.get().strip():
            self.theme_var.set("dark")
        self.outline_by_var.set("")
        self.metadata_labels_var.set("")
        self.viewer_info_html_var.set("")
        self.viewer_info_html_file_var.set("")
        self.tutorial_var.set(False)
        self.embed_reproducibility_info_var.set(True)
        self.min_panel_size_var.set("150")
        self.spot_size_var.set("auto")
        self.numba_jit_var.set(False)
        self.feature_encoding_var.set("auto")
        self.feature_value_encoding_var.set("uint16")
        self.feature_storage_var.set("sidecar")
        self.also_karospace_var.set(True)
        self.features_list_var.set("")
        self.feature_manifest_path_var.set("")
        self.feature_sidecar_shard_size_var.set("256")
        self.feature_sparse_zero_threshold_var.set("0.8")
        self._feature_items_by_modality = {}
        self._active_feature_modality = "rna"
        self.feature_modality_var.set("rna")
        if hasattr(self, "modalities_editor"):
            choices = list(self._inspected_feature_names_by_modality) or ["rna"]
            self.modalities_editor.set_choices(choices)
            self.modalities_editor.set_items(["rna"] if "rna" in choices else choices[:1])
        self._sync_feature_modality_choices()
        self.neighbor_stats_annotations_var.set("auto")
        self.neighbor_permutations_var.set("20")
        self.neighbor_stats_seed_var.set("0")
        self.statistics_contrast_categories_var.set("")
        self.statistics_counts_layer_var.set("counts")
        self.statistics_normalization_var.set("RC")
        self.statistics_scale_factor_var.set("10000")
        self.statistics_normalized_layer_var.set("off")
        self.statistics_min_cell_counts_var.set("0")
        self.statistics_min_feature_counts_var.set("0")
        self.statistics_n_cpus_var.set("1")
        self.wilcoxon_enabled_var.set(True)
        self.wilcoxon_mode_var.set("auto")
        self.wilcoxon_runtime_limit_var.set("00:30:00")
        self.wilcoxon_min_cells_per_group_var.set("20")
        self.wilcoxon_min_pct_expressed_var.set("0")
        self.wilcoxon_p_adjust_method_var.set("fdr_bh")
        self.wilcoxon_padj_cutoff_var.set("0.05")
        self.wilcoxon_log2fc_cutoff_var.set("1")
        self.wilcoxon_embed_top_n_per_comparison_var.set("2")
        self.wilcoxon_top_n_per_category_var.set("300")
        self.pseudobulk_enabled_var.set(False)
        self.pseudobulk_replicate_annotation_var.set("")
        if hasattr(self, "pseudobulk_replicate_combo"):
            self.pseudobulk_replicate_combo.configure(values=[""])
        self.pseudobulk_min_cells_per_pseudobulk_var.set("20")
        self.pseudobulk_min_replicates_var.set("2")
        self.pseudobulk_min_pct_expressed_var.set("0")
        self.pseudobulk_p_adjust_method_var.set("fdr_bh")
        self.pseudobulk_padj_cutoff_var.set("0.05")
        self.pseudobulk_log2fc_cutoff_var.set("1")
        self.pseudobulk_deseq2_fit_type_var.set("parametric")
        self.pseudobulk_embed_top_n_per_comparison_var.set("2")
        self.pathway_enabled_var.set(False)
        self.pathway_gmt_var.set("")
        self.pathway_organism_var.set("Mouse")
        self.pathway_top_n_var.set("10")
        self.pathway_min_overlap_var.set("3")
        self.pathway_gsea_permutations_var.set("100")
        self.interaction_markers_enabled_var.set(True)
        self.interaction_markers_top_targets_var.set("5")
        self.interaction_markers_top_features_var.set("20")
        self.interaction_markers_min_cells_var.set("30")
        self.interaction_markers_min_neighbors_var.set("1")
        self.section_rotations_var.set("")
        self.deconvolutions_var.set("")
        self.feature_correlation_top_n_var.set("5")
        self.spatial_variable_features_n_var.set("20")
        self.scalebar_unit_var.set("μm")
        self.section_images_var.set("")
        self.section_images_max_px_var.set("4096")
        self.section_groupby_var.set("sample_id")
        self.initial_color_var.set("leiden")
        self.outline_by_var.set("")
        self.min_panel_size_var.set("150")
        self.additional_colors_editor.set_items([])
        self.statistics_annotations_editor.set_items([])
        self.section_metadata_editor.set_items([])
        self.section_metadata_extra_editor.set_items([])
        self.manual_genes_editor.set_items([])
        self.statistics_modalities_editor.set_items([])
        self.pseudobulk_embed_top_n_per_comparison_var.set("2")
        self.neighbor_stats_annotations_var.set("auto")
        self.neighbor_permutations_var.set("20")
        self.interaction_markers_enabled_var.set(True)
        self.status_var.set("Ready")
        self._sync_required_cell_annotations()
        self._sync_statistics_layer_choices()

        if log:
            self._log("Applied default input values.")

    def _choose_file(self, variable: tk.StringVar, optional: bool = False) -> None:
        initial_dir = str(Path(variable.get()).expanduser().parent) if variable.get() else str(Path.home() / "Downloads")
        if variable is self.h5ad_var:
            path = filedialog.askopenfilename(
                initialdir=initial_dir,
                filetypes=[("AnnData files", "*.h5ad"), ("All files", "*.*")],
            )
        else:
            path = filedialog.askopenfilename(initialdir=initial_dir)
        if path:
            if variable is self.h5ad_var:
                selected_path = Path(path).expanduser()
                if selected_path.is_dir() and not self._is_zarr_path(selected_path):
                    messagebox.showerror("Invalid input", "Select a .h5ad file or a Zarr/SpatialData directory.")
                    return
                self._inspect_locked_after_press = False
                self._has_inspected_input_file = False
                self._clear_inspection_metadata()
                self._spatialdata_table_choices = []
                self._spatialdata_table_source_path = None
            variable.set(path)
        elif not optional and not variable.get():
            self._log("File selection canceled.")

    def _choose_pathway_gmt_files(self) -> None:
        existing_paths = self._parse_csv_list(self.pathway_gmt_var.get()) or []
        initial_dir = str(Path(existing_paths[0]).expanduser().parent) if existing_paths else str(Path.home() / "Downloads")
        paths = filedialog.askopenfilenames(
            initialdir=initial_dir,
            filetypes=[("GMT files", "*.gmt"), ("All files", "*.*")],
        )
        if paths:
            self.pathway_gmt_var.set(",".join(paths))

    def _choose_input_directory(self) -> None:
        current = self.h5ad_var.get().strip()
        if current:
            current_path = Path(current).expanduser()
            initial_dir = str(current_path if current_path.is_dir() else current_path.parent)
        else:
            initial_dir = str(Path.home() / "Downloads")
        path = filedialog.askdirectory(
            initialdir=initial_dir,
            title="Select SpatialData .zarr directory",
        )
        if not path:
            self._log("Directory selection canceled.")
            return
        selected = Path(path).expanduser()
        if not self._is_zarr_path(selected):
            messagebox.showerror("Invalid input", "Select a Zarr/SpatialData directory.")
            return
        self._inspect_locked_after_press = False
        self._has_inspected_input_file = False
        self._clear_inspection_metadata()
        self._spatialdata_table_choices = []
        self._spatialdata_table_source_path = None
        self.h5ad_var.set(str(selected))

    def _choose_dir(self, variable: tk.StringVar) -> None:
        initial_dir = variable.get() or str(Path.home() / "Downloads")
        path = filedialog.askdirectory(initialdir=initial_dir)
        if path:
            variable.set(path)

    def _choose_output_file(self) -> None:
        base_dir = Path(self.outdir_var.get().strip() or Path.cwd()).expanduser()
        current = Path(self.output_html_var.get().strip() or "karospace.html").expanduser()
        initial_dir = current.parent if current.is_absolute() else base_dir
        initial_file = current.name if current.name else "karospace.html"
        path = filedialog.asksaveasfilename(
            initialdir=str(initial_dir),
            initialfile=initial_file,
            defaultextension=".html",
            filetypes=[
                ("KaroSpace outputs", "*.html *.karospace"),
                ("HTML files", "*.html"),
                ("KaroSpace packages", "*.karospace"),
                ("All files", "*.*"),
            ],
        )
        if path:
            selected = Path(path).expanduser()
            try:
                self.output_html_var.set(str(selected.relative_to(base_dir)))
            except ValueError:
                self.output_html_var.set(str(selected))

    def _parameter_state_variable_names(self) -> list[str]:
        names: list[str] = []
        for name, value in self.__dict__.items():
            if name == "status_var" or not name.endswith("_var"):
                continue
            if isinstance(value, tk.Variable):
                names.append(name)
        return sorted(names)

    @staticmethod
    def _coerce_bool_state_value(value: object) -> bool:
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)):
            return bool(value)
        if isinstance(value, str):
            return value.strip().lower() in {"1", "true", "yes", "on", "checked"}
        return bool(value)

    @staticmethod
    def _json_string_list(value: object) -> list[str]:
        if not isinstance(value, list):
            return []
        items: list[str] = []
        seen: set[str] = set()
        for raw in value:
            text = str(raw).strip()
            if not text or text in seen:
                continue
            seen.add(text)
            items.append(text)
        return items

    def _collect_parameter_state(self) -> dict[str, object]:
        self._save_active_feature_selection()
        variables: dict[str, object] = {}
        for name in self._parameter_state_variable_names():
            variable = getattr(self, name)
            variables[name] = variable.get()

        lists: dict[str, list[str]] = {}
        for name in self._PARAMETER_LIST_EDITORS:
            editor = getattr(self, name, None)
            if editor is not None and hasattr(editor, "get_items"):
                lists[name] = editor.get_items()

        feature_items_by_modality = {
            str(modality): self._json_string_list(items)
            for modality, items in self._feature_items_by_modality.items()
        }
        return {
            "schema": self._PARAMETER_STATE_SCHEMA,
            "version": 1,
            "saved_at": datetime.now().isoformat(timespec="seconds"),
            "variables": variables,
            "lists": lists,
            "feature_items_by_modality": feature_items_by_modality,
            "active_feature_modality": self._active_feature_modality,
        }

    def _validate_parameter_state_format(self, state: object) -> dict[str, object]:
        if not isinstance(state, dict):
            raise ValueError("Parameter JSON must contain an object.")
        schema = state.get("schema")
        if schema != self._PARAMETER_STATE_SCHEMA:
            raise ValueError(f"Unsupported parameter JSON schema: {schema!r}.")
        version = state.get("version")
        if version != 1:
            raise ValueError(f"Unsupported parameter JSON version: {version!r}.")
        variables = state.get("variables")
        if not isinstance(variables, dict):
            raise ValueError("Parameter JSON is missing a variables object.")
        lists = state.get("lists", {})
        if lists is not None and not isinstance(lists, dict):
            raise ValueError("Parameter JSON lists must be an object.")
        feature_items = state.get("feature_items_by_modality", {})
        if feature_items is not None and not isinstance(feature_items, dict):
            raise ValueError("Parameter JSON feature_items_by_modality must be an object.")
        active_modality = state.get("active_feature_modality", "rna")
        if active_modality is not None and not isinstance(active_modality, str):
            raise ValueError("Parameter JSON active_feature_modality must be a string.")
        return state

    def _save_parameter_state_json(self) -> None:
        current_output = Path(self.output_html_var.get().strip() or "karospace.html").expanduser()
        initial_file = f"{current_output.stem or 'karospace'}-parameters.json"
        initial_dir = Path(self.outdir_var.get().strip() or Path.cwd()).expanduser()
        path = filedialog.asksaveasfilename(
            initialdir=str(initial_dir),
            initialfile=initial_file,
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Save KaroSpaceBuilder parameters",
        )
        if not path:
            self._log("Parameter JSON save canceled.")
            return
        state = self._collect_parameter_state()
        target = Path(path).expanduser()
        try:
            target.write_text(json.dumps(state, indent=2, sort_keys=True), encoding="utf-8")
        except Exception as exc:
            messagebox.showerror("Save failed", f"Could not save parameter JSON:\n{exc}")
            return
        self._log(f"Saved parameter JSON: {target}")

    def _apply_parameter_state(self, state: dict[str, object]) -> None:
        state = self._validate_parameter_state_format(state)
        variables = state["variables"]

        incoming_input = str(variables.get("h5ad_var") or self.h5ad_var.get() or "").strip()
        try:
            incoming_path = Path(incoming_input).expanduser() if incoming_input else None
        except Exception:
            incoming_path = None
        preserve_inspection = (
            self._has_inspected_input_file
            and incoming_path is not None
            and self._inspected_h5ad_path is not None
            and incoming_path == self._inspected_h5ad_path
        )

        self._loading_parameter_state = True
        try:
            self._apply_preset("default", log=False)
            if not preserve_inspection:
                self._has_inspected_input_file = False
                self._clear_inspection_metadata()
            for name, value in variables.items():
                variable = getattr(self, str(name), None)
                if not isinstance(variable, tk.Variable) or name == "status_var":
                    continue
                if isinstance(variable, tk.BooleanVar):
                    variable.set(self._coerce_bool_state_value(value))
                else:
                    variable.set("" if value is None else str(value))

            lists = state.get("lists", {})
            if isinstance(lists, dict):
                for name in self._PARAMETER_LIST_EDITORS:
                    if name == "manual_genes_editor":
                        continue
                    editor = getattr(self, name, None)
                    if editor is not None and hasattr(editor, "set_items"):
                        editor.set_items(self._json_string_list(lists.get(name)))

            feature_items_raw = state.get("feature_items_by_modality", {})
            feature_items_by_modality: dict[str, list[str]] = {}
            if isinstance(feature_items_raw, dict):
                for modality, items in feature_items_raw.items():
                    modality_name = str(modality).strip()
                    if modality_name:
                        feature_items_by_modality[modality_name] = self._json_string_list(items)
            self._feature_items_by_modality = feature_items_by_modality

            active_modality = str(state.get("active_feature_modality") or self.feature_modality_var.get() or "rna").strip() or "rna"
            self._active_feature_modality = active_modality
            self.feature_modality_var.set(active_modality)
            manual_features = feature_items_by_modality.get(active_modality)
            if manual_features is None and isinstance(lists, dict):
                manual_features = self._json_string_list(lists.get("manual_genes_editor"))
            if manual_features is None:
                manual_features = []
            if hasattr(self, "manual_genes_editor"):
                self.manual_genes_editor.set_items(manual_features)
            if manual_features and active_modality not in self._feature_items_by_modality:
                self._feature_items_by_modality[active_modality] = manual_features

            self._sync_required_cell_annotations()
            self._sync_required_section_metadata()
            self._sync_statistics_layer_choices(preserve_current=True)
            self._sync_feature_modality_choices()
            self._sync_analysis_controls()
            self._sync_app_theme_to_viewer_setting()
            self._refresh_input_gate()
            self._update_export_estimate()
        finally:
            self._loading_parameter_state = False

    def _import_parameter_state_json(self) -> None:
        initial_dir = Path(self.outdir_var.get().strip() or Path.cwd()).expanduser()
        path = filedialog.askopenfilename(
            initialdir=str(initial_dir),
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Import KaroSpaceBuilder parameters",
        )
        if not path:
            self._log("Parameter JSON import canceled.")
            return
        source = Path(path).expanduser()
        try:
            state = json.loads(source.read_text(encoding="utf-8"))
            self._apply_parameter_state(state)
        except Exception as exc:
            messagebox.showerror("Import failed", f"Could not import parameter JSON:\n{exc}")
            return
        self._log(f"Imported parameter JSON: {source}")
        self._inspect_imported_parameter_input()

    def _inspect_imported_parameter_input(self) -> None:
        valid_path = self._valid_h5ad_path()
        if valid_path is None:
            messagebox.showerror(
                "Input unavailable",
                "Parameter JSON was imported, but the input file/store path is missing or unavailable.",
            )
            return
        self._inspect_locked_after_press = False
        self._inspect_h5ad()

    def _log(self, message: str) -> None:
        stamp = datetime.now().strftime("%H:%M:%S")
        line = f"[{stamp}] {message}\n"
        self.log_text.configure(state="normal")
        self.log_text.insert("end", line)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    @staticmethod
    def _inspect_loaded_dataset(
        adata: object,
        *,
        path: Path,
        spatialdata_table: str | None,
    ) -> dict[str, object]:
        import importlib

        inspect_input_file = getattr(importlib.import_module("karospace"), "inspect_input_file")
        summary = dict(inspect_input_file(adata))
        summary["path"] = str(path)
        summary["spatialdata_table"] = spatialdata_table
        return summary

    @staticmethod
    def _format_dataset_inspection_output(summary: dict[str, object]) -> str:
        return json.dumps(summary, indent=2, default=str)

    def _inspect_h5ad(self) -> None:
        path_text = self.h5ad_var.get().strip()
        if not path_text:
            messagebox.showerror("Missing input", "Pick an input .h5ad file or SpatialData .zarr directory.")
            return

        path = Path(path_text).expanduser().resolve()
        if not path.exists():
            messagebox.showerror("Missing file", f"Input file not found:\n{path}")
            return

        self._inspect_locked_after_press = True
        self._refresh_input_gate()
        self._set_inspect_loading(True)
        self._log(f"Inspecting {path}")
        adata = None
        requested_table = None
        try:
            if self._is_zarr_path(path):
                import importlib

                self._import_karospace_api(enable_numba_jit=False)
                table_choices = self._discover_spatialdata_tables(path)
                self._spatialdata_table_choices = table_choices
                self._spatialdata_table_source_path = path
                requested_table = self._parse_optional_text(self.spatialdata_table_var.get())
                if table_choices:
                    if requested_table not in table_choices:
                        requested_table = table_choices[0]
                    self._syncing_spatialdata_table = True
                    try:
                        self.spatialdata_table_var.set(requested_table)
                    finally:
                        self._syncing_spatialdata_table = False
                    self._sync_spatialdata_table_selector(path)
                    table_note = (
                        f"Selected SpatialData table '{requested_table}' "
                        f"({len(table_choices)} available: {', '.join(table_choices)})."
                    )
                    self._log(table_note)
                else:
                    requested_table = None
                    self._syncing_spatialdata_table = True
                    try:
                        self.spatialdata_table_var.set("")
                    finally:
                        self._syncing_spatialdata_table = False
                    self._sync_spatialdata_table_selector(path)
                    self._log("No SpatialData tables directory detected; reading .zarr as AnnData or default SpatialData input.")

                data_loader = importlib.import_module("karospace.data_loader")
                coerce = getattr(data_loader, "_coerce_input_to_anndata")
                adata, _source_label, _table_key = coerce(
                    str(path),
                    requested_table,
                )
            else:
                self._spatialdata_table_choices = []
                self._sync_spatialdata_table_selector(path)
                ad_mod = _get_anndata()
                try:
                    adata = ad_mod.read_h5ad(path, backed="r")
                except Exception:
                    adata = ad_mod.read_h5ad(path)

            obs_cols = [str(c) for c in adata.obs.columns]
            obs_col_set = set(obs_cols)
            section_key_cols = self._eligible_section_key_columns(adata, obs_cols)
            section_key_col_set = set(section_key_cols)
            obsm_keys = self._obsm_keys_from_adata(adata)
            spatial_obs_cols = self._spatial_obs_columns_from_adata(adata, obs_cols)
            layer_keys = self._layer_keys_from_adata(adata)
            total_var_count = int(adata.n_vars)
            max_gene_choices = 120000
            feature_names_by_modality = self._feature_names_by_modality_from_adata(
                adata,
                max_features=max_gene_choices,
            )
            if total_var_count > max_gene_choices:
                self._log(
                    f"Large gene table detected ({total_var_count}). "
                    f"Loaded first {max_gene_choices} genes into pickers for responsiveness."
                )

            self.additional_colors_editor.set_choices(section_key_cols)
            self.statistics_annotations_editor.set_choices(obs_cols)
            self.section_metadata_editor.set_choices(section_key_cols)
            self.section_metadata_extra_editor.set_choices(section_key_cols)
            self._set_feature_modalities(feature_names_by_modality)
            modality_choices = list(feature_names_by_modality)
            if hasattr(self, "groupby_combo"):
                self.groupby_combo.configure(values=[""] + section_key_cols)
            if hasattr(self, "color_combo"):
                self.color_combo.configure(values=section_key_cols)
            if hasattr(self, "outline_combo"):
                self.outline_combo.configure(values=[""] + section_key_cols)
            if hasattr(self, "spatial_key_combo"):
                self.spatial_key_combo.configure(values=obsm_keys)
            if hasattr(self, "spatial_x_combo"):
                self.spatial_x_combo.configure(values=[""] + spatial_obs_cols)
            if hasattr(self, "spatial_y_combo"):
                self.spatial_y_combo.configure(values=[""] + spatial_obs_cols)
            if hasattr(self, "pseudobulk_replicate_combo"):
                self.pseudobulk_replicate_combo.configure(values=[""] + obs_cols)

            current_section_groupby = self.section_groupby_var.get().strip()
            if current_section_groupby and current_section_groupby not in section_key_col_set:
                if "sample_id" in section_key_col_set:
                    self.section_groupby_var.set("sample_id")
                elif section_key_cols:
                    self.section_groupby_var.set(section_key_cols[0])
                else:
                    self.section_groupby_var.set("")
            section_groupby = self.section_groupby_var.get().strip()
            section_counts = (
                [int(value) for value in adata.obs[section_groupby].value_counts(dropna=False).tolist()]
                if section_groupby in section_key_col_set
                else []
            )
            section_values = self._section_values_from_adata(adata, section_groupby) if section_groupby in section_key_col_set else []
            initial_color = self.initial_color_var.get().strip()
            initial_in_obs = initial_color in section_key_col_set
            if not initial_in_obs:
                clustering_color = next(
                    (
                        column
                        for column in section_key_cols
                        if column.lower().startswith(("louvain", "leiden"))
                    ),
                    None,
                )
                if clustering_color:
                    self.initial_color_var.set(clustering_color)
                elif section_key_cols:
                    self.initial_color_var.set(section_key_cols[0])
                else:
                    self.initial_color_var.set("")
            if self.outline_by_var.get().strip() and self.outline_by_var.get().strip() not in section_key_col_set:
                self.outline_by_var.set("")
            if self.pseudobulk_replicate_annotation_var.get().strip() not in obs_col_set:
                self.pseudobulk_replicate_annotation_var.set("")

            existing_additional = [name for name in self.additional_colors_editor.get_items() if name in section_key_col_set]
            self.additional_colors_editor.set_items(existing_additional)

            existing_statistics_annotations = [
                name for name in self.statistics_annotations_editor.get_items() if name in obs_col_set
            ]
            self.statistics_annotations_editor.set_items(existing_statistics_annotations)

            selected_modality_choices = set(self._feature_modality_choices())
            existing_statistics_modalities = [
                name for name in self.statistics_modalities_editor.get_items() if name in selected_modality_choices
            ]
            self.statistics_modalities_editor.set_items(existing_statistics_modalities)
            self._sync_statistics_modality_choices()

            existing_section_metadata = [
                name for name in self.section_metadata_editor.get_items() if name in section_key_col_set
            ]
            self.section_metadata_editor.set_items(existing_section_metadata)

            existing_section_metadata_extra = [
                name for name in self.section_metadata_extra_editor.get_items() if name in section_key_col_set
            ]
            self.section_metadata_extra_editor.set_items(existing_section_metadata_extra)

            if obs_cols:
                self._log(f"Loaded {len(obs_cols)} obs columns into analysis and metadata pickers.")
            if section_key_cols:
                self._log(
                    f"Loaded {len(section_key_cols)} low-cardinality choices with fewer than 500 unique values."
                )
            else:
                self._log("No low-cardinality choices found with fewer than 500 unique values.")
            if feature_names_by_modality:
                self._log(
                    "Loaded feature pickers for modalities: "
                    + ", ".join(
                        f"{name} ({len(features)})"
                        for name, features in feature_names_by_modality.items()
                    )
                    + "."
                )

            if self.spatial_key_var.get().strip() not in set(obsm_keys):
                self.spatial_key_var.set("spatial" if "spatial" in obsm_keys else obsm_keys[0])
            spatial_key = self.spatial_key_var.get().strip() or "spatial"
            has_spatial = spatial_key in adata.obsm
            has_centroid = {"centroid_x", "centroid_y"}.issubset(set(obs_cols))
            if has_spatial:
                self.coords_var.set("obsm:spatial")
                self.spatial_x_var.set("")
                self.spatial_y_var.set("")
                inspected_coords_mode = "obsm:spatial"
            elif has_centroid:
                self.coords_var.set("obs:centroid_x_y")
                self.spatial_x_var.set("centroid_x")
                self.spatial_y_var.set("centroid_y")
                inspected_coords_mode = "obs:centroid_x_y"
            else:
                self.coords_var.set("auto")
                if self.spatial_x_var.get().strip() not in set(spatial_obs_cols):
                    self.spatial_x_var.set("")
                if self.spatial_y_var.get().strip() not in set(spatial_obs_cols):
                    self.spatial_y_var.set("")
                inspected_coords_mode = None

            self._inspected_h5ad_path = path
            self._inspected_coords_mode = inspected_coords_mode
            self._inspected_n_cells = int(adata.n_obs)
            self._inspected_n_genes = int(adata.n_vars)
            self._inspected_obs_cols = set(obs_cols)
            self._inspected_section_key_cols = section_key_cols
            self._inspected_feature_names_by_modality = feature_names_by_modality
            self._inspected_obsm_keys = obsm_keys
            self._inspected_spatial_obs_cols = spatial_obs_cols
            self._inspected_layer_keys = layer_keys
            self._inspected_section_counts_by_column = {}
            self._inspected_section_values_by_column = {}
            if section_groupby and section_counts:
                self._inspected_section_counts_by_column[section_groupby] = section_counts
            if section_groupby and section_values:
                self._inspected_section_values_by_column[section_groupby] = section_values
            self._has_inspected_input_file = True
            self._sync_required_cell_annotations()
            self._sync_statistics_layer_choices()
            self._refresh_section_order_values(preserve_current=True)

            try:
                inspection_output = self._inspect_loaded_dataset(
                    adata,
                    path=path,
                    spatialdata_table=requested_table,
                )
            except Exception as exc:
                self._log(f"Dataset inspect output unavailable: {exc}")
            else:
                self._log(
                    "Dataset inspect output:\n"
                    f"{self._format_dataset_inspection_output(inspection_output)}"
                )

            self._log(
                f"obs columns: {len(obs_cols)} | cells: {adata.n_obs} | genes: {adata.n_vars} | "
                f"coords: {f'obsm:{spatial_key}' if has_spatial else 'obs centroids' if has_centroid else 'not detected'}"
            )
            self.status_var.set("Inspection complete")
            self._refresh_input_gate()
            self._update_export_estimate()
        except Exception as exc:
            self._has_inspected_input_file = False
            self._clear_inspection_metadata()
            self._refresh_input_gate()
            messagebox.showerror("Inspect failed", str(exc))
            self._log(f"Inspect failed: {exc}")
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()
            self._set_inspect_loading(False)

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
    def _parse_non_negative_float(label: str, raw: str) -> float:
        text = str(raw).strip()
        try:
            value = float(text)
        except ValueError as exc:
            raise ValueError(f"{label} must be a number.") from exc
        if value < 0:
            raise ValueError(f"{label} must be >= 0.")
        return value

    @staticmethod
    def _parse_probability(label: str, raw: str) -> float:
        value = ExportApp._parse_non_negative_float(label, raw)
        if value > 1:
            raise ValueError(f"{label} must be between 0 and 1.")
        return value

    @staticmethod
    def _parse_csv_list(raw: str) -> list[str] | None:
        values = [item.strip() for item in str(raw or "").split(",") if item.strip()]
        return values or None

    @classmethod
    def _parse_modalities_option(cls, raw: str, label: str) -> list[str] | str | None:
        text = str(raw or "").strip()
        if not text:
            return None
        lowered = text.lower()
        if lowered == "all":
            return "all"
        if lowered in {"none", "null"}:
            return None
        values = cls._parse_csv_list(text)
        if values is None:
            return None
        return values

    @staticmethod
    def _parse_optional_text(raw: str) -> str | None:
        text = str(raw or "").strip()
        if not text or text.lower() in {"none", "null", "off"}:
            return None
        return text

    @staticmethod
    def _parse_auto_or_none(raw: str, label: str) -> str | None:
        text = str(raw or "").strip().lower()
        if not text or text == "auto":
            return "auto"
        if text in {"none", "null"}:
            return None
        raise ValueError(f"{label} must be 'auto' or 'None'.")

    @classmethod
    def _parse_neighbor_stats_annotations(cls, raw: str, main_cell_annotation: str) -> list[str] | None:
        text = str(raw or "").strip()
        if text.lower() == "auto":
            return [main_cell_annotation]
        if not text:
            return None
        return cls._parse_csv_list(text)

    @staticmethod
    def _clean_statistics_category_list(value: object, label: str) -> list[str] | None:
        if value is None:
            return None
        if isinstance(value, str):
            categories = [item.strip() for item in value.split(",") if item.strip()]
            return categories or None
        if isinstance(value, (list, tuple, set)):
            categories = [str(item).strip() for item in value if str(item).strip()]
            return categories or None
        raise ValueError(f"{label} values must be category strings or lists of category strings.")

    @classmethod
    def _normalize_statistics_contrast_categories(
        cls,
        raw: str,
        annotation_columns: list[str],
        *,
        label: str = "Statistics contrast categories",
    ) -> dict[str, list[str] | None] | None:
        columns = [str(column).strip() for column in annotation_columns if str(column).strip()]
        text = str(raw or "").strip()
        if not columns or not text:
            return None

        value: object = text
        if text[0] in "[{":
            try:
                value = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"{label} must be valid JSON: {exc}") from exc
        elif len(columns) == 1:
            return {columns[0]: cls._clean_statistics_category_list(text, label)}
        else:
            raise ValueError(
                f"{label} is ambiguous with multiple statistics annotations. "
                "Use a JSON object keyed by annotation name or a nested JSON list in order: "
                + ", ".join(columns)
                + "."
            )

        if isinstance(value, dict):
            unknown = sorted(str(key) for key in value.keys() if str(key) not in set(columns))
            if unknown:
                raise ValueError(
                    f"{label} contains annotation column(s) not requested for statistics: "
                    + ", ".join(unknown)
                )
            return {
                column: cls._clean_statistics_category_list(value.get(column), label)
                for column in columns
            }

        if isinstance(value, list):
            nested = any(isinstance(item, (dict, list, tuple, set)) for item in value)
            if not nested:
                if len(columns) > 1:
                    raise ValueError(
                        f"{label} is ambiguous with multiple statistics annotations. "
                        "Use a JSON object keyed by annotation name or a nested JSON list in order: "
                        + ", ".join(columns)
                        + "."
                    )
                return {columns[0]: cls._clean_statistics_category_list(value, label)}
            if len(value) != len(columns):
                raise ValueError(
                    f"{label} nested list must contain one category list per statistics annotation "
                    f"({len(columns)} expected: {', '.join(columns)})."
                )
            normalized: dict[str, list[str] | None] = {}
            for column, item in zip(columns, value):
                if isinstance(item, dict):
                    if "categories" in item:
                        item = item.get("categories")
                    elif column in item:
                        item = item.get(column)
                    else:
                        raise ValueError(
                            f"{label} nested mapping for '{column}' must contain 'categories' or the annotation name."
                        )
                normalized[column] = cls._clean_statistics_category_list(item, label)
            return normalized

        raise ValueError(
            f"{label} must be a category list, a JSON object keyed by annotation, "
            "or a nested list matching the statistics annotation order."
        )

    @staticmethod
    def _parse_json_mapping(raw: str, label: str) -> dict | None:
        text = str(raw or "").strip()
        if not text:
            return None
        try:
            parsed = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"{label} must be valid JSON: {exc}") from exc
        if not isinstance(parsed, dict):
            raise ValueError(f"{label} must be a JSON object.")
        return parsed

    @classmethod
    def _parse_metadata_value_order(cls, raw: str) -> dict[str, list[str]] | None:
        parsed = cls._parse_json_mapping(raw, "Metadata value order JSON")
        if parsed is None:
            return None
        out: dict[str, list[str]] = {}
        for key, values in parsed.items():
            if not isinstance(values, list):
                raise ValueError("Metadata value order JSON values must be lists.")
            out[str(key)] = [str(value) for value in values]
        return out or None

    @classmethod
    def _parse_metadata_labels(cls, raw: str) -> dict[str, str] | None:
        parsed = cls._parse_json_mapping(raw, "Metadata labels JSON")
        if parsed is None:
            return None
        out = {str(key): str(value) for key, value in parsed.items() if str(key).strip() and value is not None}
        return out or None

    @classmethod
    def _parse_string_mapping(cls, raw: str, label: str) -> dict[str, str] | None:
        parsed = cls._parse_json_mapping(raw, label)
        if parsed is None:
            return None
        out = {str(key): str(value) for key, value in parsed.items() if str(key).strip() and value is not None}
        return out or None

    @classmethod
    def _parse_section_images(cls, raw: str) -> dict[str, object] | None:
        parsed = cls._parse_json_mapping(raw, "Section images JSON")
        return parsed or None

    @staticmethod
    def _parse_section_rotations(raw: str) -> dict[str, float] | None:
        text = str(raw or "").strip()
        if not text:
            return None
        if text.startswith("{"):
            try:
                parsed = json.loads(text)
            except json.JSONDecodeError as exc:
                raise ValueError(f"Section rotations JSON must be valid JSON: {exc}") from exc
            if not isinstance(parsed, dict):
                raise ValueError("Section rotations JSON must be an object.")
            source = parsed.items()
        else:
            items = []
            for token in text.split(","):
                token = token.strip()
                if not token:
                    continue
                if ":" not in token:
                    raise ValueError("Section rotations must use section_id:angle CSV or JSON object syntax.")
                key, value = token.split(":", 1)
                items.append((key.strip(), value.strip()))
            source = items

        out: dict[str, float] = {}
        for key, value in source:
            section_id = str(key).strip()
            if not section_id:
                continue
            try:
                out[section_id] = float(value)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Invalid rotation angle for section {section_id!r}.") from exc
        return out or None

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

    def _read_adata_for_feature_ops(self, path: Path):
        if self._is_zarr_path(path):
            import importlib

            self._import_karospace_api(enable_numba_jit=False)
            spatialdata_table = self._parse_optional_text(self.spatialdata_table_var.get())
            if spatialdata_table is None:
                if not self._spatialdata_table_choices:
                    self._spatialdata_table_choices = self._discover_spatialdata_tables(path)
                    self._spatialdata_table_source_path = path
                if self._spatialdata_table_choices:
                    spatialdata_table = self._spatialdata_table_choices[0]
                    self._syncing_spatialdata_table = True
                    try:
                        self.spatialdata_table_var.set(spatialdata_table)
                    finally:
                        self._syncing_spatialdata_table = False
            data_loader = importlib.import_module("karospace.data_loader")
            coerce = getattr(data_loader, "_coerce_input_to_anndata")
            adata, _source_label, _table_key = coerce(
                str(path),
                spatialdata_table,
            )
            return adata

        ad_mod = _get_anndata()
        try:
            return ad_mod.read_h5ad(path, backed="r")
        except Exception:
            return ad_mod.read_h5ad(path)

    def _load_feature_names_by_modality(self, h5ad_path: Path) -> dict[str, set[str]]:
        if self._matches_inspected_h5ad(h5ad_path) and self._inspected_feature_names_by_modality:
            return {
                modality: set(features)
                for modality, features in self._inspected_feature_names_by_modality.items()
            }

        adata = None
        try:
            adata = self._read_adata_for_feature_ops(h5ad_path)
            return {
                modality: set(features)
                for modality, features in self._feature_names_by_modality_from_adata(
                    adata,
                    max_features=None,
                ).items()
            }
        finally:
            if adata is not None and getattr(adata, "isbacked", False):
                file_obj = getattr(adata, "file", None)
                if file_obj is not None:
                    file_obj.close()

    def _resolve_features(
        self,
        h5ad_path: Path,
        *,
        features_by_modality: dict[str, list[str]] | None = None,
        require_features: bool = True,
    ) -> list[str] | None:
        selected = features_by_modality if features_by_modality is not None else self._selected_features_by_modality()
        genes = self._merge_unique(*selected.values())
        if not genes:
            if require_features:
                raise ValueError("Add at least one feature or choose a features_list file.")
            return None

        feature_names_by_modality = self._load_feature_names_by_modality(h5ad_path)
        missing_by_modality: dict[str, list[str]] = {}
        for modality, features in selected.items():
            allowed = feature_names_by_modality.get(modality)
            if allowed is None:
                missing_by_modality[modality] = features
                continue
            missing = [feature for feature in features if feature not in allowed]
            if missing:
                missing_by_modality[modality] = missing
        if missing_by_modality:
            modality, missing = next(iter(missing_by_modality.items()))
            preview = ", ".join(missing[:10])
            raise ValueError(f"{len(missing)} features are missing in {modality} feature names: {preview}")
        return genes

    def _parse_config(self) -> BuilderConfig:
        h5ad_text = self.h5ad_var.get().strip()
        outdir_text = self.outdir_var.get().strip()
        output_text = self.output_html_var.get().strip()

        if not h5ad_text:
            raise ValueError("Input file is required.")
        if not outdir_text:
            raise ValueError("Output directory is required.")
        if not output_text:
            raise ValueError("Output path is required.")

        h5ad_path = Path(h5ad_text).expanduser().resolve()
        outdir = Path(outdir_text).expanduser().resolve()
        output_candidate = Path(output_text).expanduser()
        output_html_path = output_candidate if output_candidate.is_absolute() else outdir / output_candidate
        if output_html_path.suffix.lower() not in {".html", ".karospace"}:
            output_html_path = output_html_path.with_suffix(".html")
        if not h5ad_path.exists():
            raise ValueError(f"Input file not found: {h5ad_path}")

        output_html_path.parent.mkdir(parents=True, exist_ok=True)
        outdir = output_html_path.parent

        coords_raw = self.coords_var.get().strip().lower() or "auto"
        if coords_raw not in {"auto", "obsm:spatial", "obs:centroid_x_y"}:
            raise ValueError("Coordinates must be auto, obsm:spatial, or obs:centroid_x_y.")
        coords_mode = None if coords_raw == "auto" else coords_raw
        spatial_key = self.spatial_key_var.get().strip() or "spatial"
        spatial_x = self._parse_optional_text(self.spatial_x_var.get())
        spatial_y = self._parse_optional_text(self.spatial_y_var.get())
        if bool(spatial_x) != bool(spatial_y):
            raise ValueError("Spatial X and Spatial Y must be provided together.")
        spatial_columns = (spatial_x, spatial_y) if spatial_x and spatial_y else None
        if spatial_columns is not None and coords_mode == "obs:centroid_x_y":
            coords_mode = None

        spatialdata_table = self._parse_optional_text(self.spatialdata_table_var.get())
        if self._is_zarr_path(h5ad_path) and spatialdata_table is None and self._spatialdata_table_choices:
            spatialdata_table = self._spatialdata_table_choices[0]
            self._syncing_spatialdata_table = True
            try:
                self.spatialdata_table_var.set(spatialdata_table)
            finally:
                self._syncing_spatialdata_table = False
        section_groupby = self.section_groupby_var.get().strip()
        initial_color = self.initial_color_var.get().strip()
        if not initial_color:
            raise ValueError("Main cell annotation is required.")
        title = self.title_var.get().strip() or "KaroSpace"
        outline_by = self._parse_optional_text(self.outline_by_var.get())
        eligible_annotation_cols = set(self._inspected_section_key_cols)
        if eligible_annotation_cols:
            if section_groupby and section_groupby not in eligible_annotation_cols:
                raise ValueError("Section key must have fewer than 500 unique values.")
            if initial_color not in eligible_annotation_cols:
                raise ValueError("Main cell annotation must have fewer than 500 unique values.")
            if outline_by is not None and outline_by not in eligible_annotation_cols:
                raise ValueError("Outline by must have fewer than 500 unique values.")
        min_panel_size = self._parse_positive_int("Min panel size", self.min_panel_size_var.get())
        spot_size = self._parse_spot_size(self.spot_size_var.get())

        downsample_text = self.downsample_var.get().strip()
        downsample = None
        if downsample_text:
            downsample = self._parse_positive_int("Downsample", downsample_text)

        additional_colors = self._merge_unique([initial_color], self.additional_colors_editor.get_items())
        statistics_annotation_lists = self._merge_unique(self.statistics_annotations_editor.get_items())
        section_metadata = self._merge_unique(
            self._required_section_metadata_items(),
            self.section_metadata_editor.get_items(),
        )
        section_metadata_extra = self._merge_unique(self.section_metadata_extra_editor.get_items())
        if eligible_annotation_cols:
            for label, values in (
                ("cell_annotations", additional_colors),
                ("section_metadata", section_metadata),
                ("section_metadata_extra", section_metadata_extra),
            ):
                invalid = [value for value in values if value not in eligible_annotation_cols]
                if invalid:
                    raise ValueError(
                        f"{label} choices must have fewer than 500 unique values: "
                        + ", ".join(invalid[:10])
                    )
        features_list = self._parse_optional_text(self.features_list_var.get())
        if features_list is not None and not Path(features_list).expanduser().is_file():
            raise ValueError(f"Features list file not found: {features_list}")
        embed_modalities = self._selected_embed_modalities()
        if self._inspected_feature_names_by_modality and embed_modalities:
            available_modalities = set(self._inspected_feature_names_by_modality)
            invalid_modalities = [modality for modality in embed_modalities if modality not in available_modalities]
            if invalid_modalities:
                raise ValueError("Modality choices are not available in the inspected dataset: " + ", ".join(invalid_modalities[:10]))
        selected_features_by_modality = self._selected_features_by_modality()

        genes = self._resolve_features(
            h5ad_path,
            features_by_modality=selected_features_by_modality,
            require_features=features_list is None,
        )

        feature_encoding = self.feature_encoding_var.get().strip().lower() or "auto"
        if feature_encoding not in {"auto", "dense", "sparse"}:
            raise ValueError("Feature encoding must be auto, dense, or sparse.")
        feature_value_encoding = self.feature_value_encoding_var.get().strip().lower() or "uint16"
        if feature_value_encoding not in {"uint16", "uint8"}:
            raise ValueError("Feature value encoding must be uint16 or uint8.")
        feature_storage = self.feature_storage_var.get().strip().lower() or "embedded"
        if feature_storage not in {"embedded", "sidecar"}:
            raise ValueError("Feature storage must be embedded or sidecar.")
        if output_html_path.suffix.lower() == ".karospace" and feature_storage != "sidecar":
            raise ValueError(".karospace output requires Feature storage to be sidecar.")
        # Only package a companion .karospace when the export produces a sidecar
        # HTML bundle (embedded storage has nothing to package, and a .karospace
        # output path is already the package itself).
        also_export_karospace = (
            bool(self.also_karospace_var.get())
            and feature_storage == "sidecar"
            and output_html_path.suffix.lower() != ".karospace"
        )
        feature_manifest_path = self._parse_optional_text(self.feature_manifest_path_var.get())
        feature_sidecar_shard_size = self._parse_positive_int(
            "Feature sidecar shard size", self.feature_sidecar_shard_size_var.get()
        )
        feature_sparse_zero_threshold = self._parse_probability(
            "Feature sparse zero threshold", self.feature_sparse_zero_threshold_var.get()
        )
        if embed_modalities:
            modalities = embed_modalities
        elif selected_features_by_modality:
            modalities = list(selected_features_by_modality.keys())
        elif features_list is not None:
            modalities = [self.feature_modality_var.get().strip() or "rna"]
        else:
            modalities = None

        neighbor_stats_groupby = self._parse_neighbor_stats_annotations(
            self.neighbor_stats_annotations_var.get(),
            initial_color,
        )
        neighbor_permutations = self._parse_neighbor_permutations(self.neighbor_permutations_var.get())
        neighbor_seed = self._parse_non_negative_int("Neighbor stats seed", self.neighbor_stats_seed_var.get() or "0")
        statistics_additional_annotations = self._merge_unique(statistics_annotation_lists) or None
        statistics_annotation_columns = self._merge_unique([initial_color], statistics_additional_annotations or [])
        statistics_modalities = self._merge_unique(self.statistics_modalities_editor.get_items()) or None
        statistics_contrast_categories = self._normalize_statistics_contrast_categories(
            self.statistics_contrast_categories_var.get(),
            statistics_annotation_columns,
        )
        statistics_counts_layer = self._parse_optional_text(self.statistics_counts_layer_var.get())
        statistics_normalization = self.statistics_normalization_var.get().strip() or "RC"
        if statistics_normalization not in {"RC", "LogNormalize"}:
            raise ValueError("Statistics normalization must be RC or LogNormalize.")
        statistics_scale_factor = self._parse_non_negative_float(
            "Statistics scale factor",
            self.statistics_scale_factor_var.get(),
        )
        if statistics_scale_factor <= 0:
            raise ValueError("Statistics scale factor must be > 0.")
        statistics_normalized_layer = self._parse_optional_text(self.statistics_normalized_layer_var.get())
        statistics_min_cell_counts = self._parse_non_negative_int(
            "Statistics min cell counts", self.statistics_min_cell_counts_var.get()
        )
        statistics_min_feature_counts = self._parse_non_negative_int(
            "Statistics min feature counts", self.statistics_min_feature_counts_var.get()
        )
        statistics_n_cpus = self._parse_positive_int("Statistics CPUs", self.statistics_n_cpus_var.get())
        wilcoxon_enabled = bool(self.wilcoxon_enabled_var.get())
        wilcoxon = self.wilcoxon_mode_var.get().strip().lower() or "auto"
        if not wilcoxon_enabled:
            wilcoxon = "off"
        if wilcoxon not in {"auto", "off", "force"}:
            raise ValueError("Wilcoxon mode must be auto, off, or force.")
        if wilcoxon_enabled:
            wilcoxon_runtime_limit = self.wilcoxon_runtime_limit_var.get().strip() or "00:30:00"
            wilcoxon_min_cells_per_group = self._parse_positive_int(
                "Wilcoxon min cells/group", self.wilcoxon_min_cells_per_group_var.get()
            )
            wilcoxon_min_pct_expressed = self._parse_non_negative_float(
                "Wilcoxon min pct expressed", self.wilcoxon_min_pct_expressed_var.get()
            )
            wilcoxon_p_adjust_method = self.wilcoxon_p_adjust_method_var.get().strip().lower() or "fdr_bh"
            if wilcoxon_p_adjust_method not in {"fdr_bh", "bonferroni", "holm", "none"}:
                raise ValueError("Wilcoxon p adjust method must be fdr_bh, bonferroni, holm, or none.")
            wilcoxon_padj_cutoff = self._parse_probability("Wilcoxon padj cutoff", self.wilcoxon_padj_cutoff_var.get())
            wilcoxon_log2fc_cutoff = self._parse_non_negative_float(
                "Wilcoxon log2FC cutoff", self.wilcoxon_log2fc_cutoff_var.get()
            )
            wilcoxon_embed_top_n_per_comparison = self._parse_non_negative_int(
                "Wilcoxon embed top/comparison", self.wilcoxon_embed_top_n_per_comparison_var.get()
            )
            wilcoxon_top_n_per_category = self._parse_positive_int(
                "Wilcoxon rows/category", self.wilcoxon_top_n_per_category_var.get()
            )
        else:
            wilcoxon_runtime_limit = "00:30:00"
            wilcoxon_min_cells_per_group = 20
            wilcoxon_min_pct_expressed = 0.0
            wilcoxon_p_adjust_method = "fdr_bh"
            wilcoxon_padj_cutoff = 0.05
            wilcoxon_log2fc_cutoff = 1.0
            wilcoxon_embed_top_n_per_comparison = 2
            wilcoxon_top_n_per_category = 300
        pseudobulk_enabled = bool(self.pseudobulk_enabled_var.get())
        pathway_enabled = bool(self.pathway_enabled_var.get())
        interaction_enabled = bool(self.interaction_markers_enabled_var.get())
        pseudobulk = "auto" if pseudobulk_enabled else None
        pathway = "auto" if pathway_enabled else None
        interaction_markers = "auto" if interaction_enabled else None
        if pseudobulk_enabled:
            pseudobulk_replicate_annotation = self._parse_optional_text(self.pseudobulk_replicate_annotation_var.get())
            if pseudobulk_replicate_annotation and self._inspected_obs_cols and pseudobulk_replicate_annotation not in self._inspected_obs_cols:
                raise ValueError("Pseudobulk replicate annotation must be an inspected obs column.")
            pseudobulk_min_cells_per_pseudobulk = self._parse_positive_int(
                "Pseudobulk min cells per sample", self.pseudobulk_min_cells_per_pseudobulk_var.get()
            )
            pseudobulk_min_replicates = self._parse_positive_int(
                "Pseudobulk min replicates", self.pseudobulk_min_replicates_var.get()
            )
            pseudobulk_min_pct_expressed = self._parse_non_negative_float(
                "Pseudobulk min pct expressed", self.pseudobulk_min_pct_expressed_var.get()
            )
            pseudobulk_p_adjust_method = self.pseudobulk_p_adjust_method_var.get().strip().lower() or "fdr_bh"
            if pseudobulk_p_adjust_method not in {"fdr_bh", "bonferroni", "holm", "none"}:
                raise ValueError("Pseudobulk p adjust method must be fdr_bh, bonferroni, holm, or none.")
            pseudobulk_padj_cutoff = self._parse_probability("Pseudobulk padj cutoff", self.pseudobulk_padj_cutoff_var.get())
            pseudobulk_log2fc_cutoff = self._parse_non_negative_float(
                "Pseudobulk log2FC cutoff", self.pseudobulk_log2fc_cutoff_var.get()
            )
            pseudobulk_deseq2_fit_type = self.pseudobulk_deseq2_fit_type_var.get().strip().lower() or "parametric"
            if pseudobulk_deseq2_fit_type not in {"parametric", "mean"}:
                raise ValueError("Pseudobulk DESeq2 fit type must be parametric or mean.")
            pseudobulk_embed_top_n_per_comparison = self._parse_non_negative_int(
                "Auto-embedded DE genes", self.pseudobulk_embed_top_n_per_comparison_var.get()
            )
        else:
            pseudobulk_replicate_annotation = None
            pseudobulk_min_cells_per_pseudobulk = 20
            pseudobulk_min_replicates = 2
            pseudobulk_min_pct_expressed = 0.0
            pseudobulk_p_adjust_method = "fdr_bh"
            pseudobulk_padj_cutoff = 0.05
            pseudobulk_log2fc_cutoff = 1.0
            pseudobulk_deseq2_fit_type = "parametric"
            pseudobulk_embed_top_n_per_comparison = 2
        if pathway_enabled:
            pathway_gmt = self._parse_csv_list(self.pathway_gmt_var.get())
            pathway_organism = self.pathway_organism_var.get().strip() or "Mouse"
            pathway_top_n = self._parse_positive_int("Pathway top N", self.pathway_top_n_var.get())
            pathway_min_overlap = self._parse_positive_int("Pathway min overlap", self.pathway_min_overlap_var.get())
            pathway_gsea_permutations = self._parse_non_negative_int(
                "Pathway GSEA permutations", self.pathway_gsea_permutations_var.get()
            )
        else:
            pathway_gmt = None
            pathway_organism = "Mouse"
            pathway_top_n = 10
            pathway_min_overlap = 3
            pathway_gsea_permutations = 100
        if interaction_enabled:
            interaction_top_targets = self._parse_positive_int(
                "Interaction top targets", self.interaction_markers_top_targets_var.get()
            )
            interaction_top_features = self._parse_positive_int(
                "Interaction top features", self.interaction_markers_top_features_var.get()
            )
            interaction_min_cells = self._parse_positive_int(
                "Interaction min cells", self.interaction_markers_min_cells_var.get()
            )
            interaction_min_neighbors = self._parse_positive_int(
                "Interaction min neighbors", self.interaction_markers_min_neighbors_var.get()
            )
        else:
            interaction_top_targets = 5
            interaction_top_features = 20
            interaction_min_cells = 10
            interaction_min_neighbors = 10
        metadata_value_order = self._parse_metadata_value_order(self.metadata_value_order_var.get())
        metadata_labels = self._parse_metadata_labels(self.metadata_labels_var.get())
        metadata_max_columns = None
        metadata_max_columns_text = self.metadata_max_columns_var.get().strip()
        if metadata_max_columns_text:
            metadata_max_columns = self._parse_positive_int("Metadata max columns", metadata_max_columns_text)
        section_order = self._parse_csv_list(self.section_order_var.get())
        viewer_info_html = self._parse_optional_text(self.viewer_info_html_var.get())
        viewer_info_path = self._parse_optional_text(self.viewer_info_html_file_var.get())
        if viewer_info_path is not None:
            path = Path(viewer_info_path).expanduser()
            if not path.exists():
                raise ValueError(f"Viewer info HTML file not found: {path}")
            viewer_info_html = path.read_text(encoding="utf-8")
        source_input_path = str(h5ad_path)
        section_rotations = self._parse_section_rotations(self.section_rotations_var.get())
        deconvolutions = self._parse_string_mapping(self.deconvolutions_var.get(), "Deconvolutions JSON")
        feature_correlation_top_n = self._parse_non_negative_int("Feature correlations", self.feature_correlation_top_n_var.get())
        spatial_variable_features_n = self._parse_non_negative_int(
            "Spatial variable features", self.spatial_variable_features_n_var.get()
        )
        scalebar_unit = self.scalebar_unit_var.get().strip() or "μm"
        section_images = self._parse_section_images(self.section_images_var.get())
        section_images_max_px = self._parse_positive_int("Section images max px", self.section_images_max_px_var.get())

        return BuilderConfig(
            h5ad_path=h5ad_path,
            outdir=outdir,
            output_html_path=output_html_path,
            coords_mode=coords_mode,
            spatial_key=spatial_key,
            spatial_columns=spatial_columns,
            spatialdata_table=spatialdata_table,
            section_groupby=section_groupby,
            section_order=section_order,
            section_metadata=section_metadata or None,
            section_metadata_extra=section_metadata_extra or None,
            metadata_value_order=metadata_value_order,
            metadata_max_columns=metadata_max_columns,
            initial_color=initial_color,
            title=title,
            outline_by=outline_by,
            metadata_labels=metadata_labels,
            viewer_info_html=viewer_info_html,
            tutorial=bool(self.tutorial_var.get()),
            embed_reproducibility_info=bool(self.embed_reproducibility_info_var.get()),
            source_input_path=source_input_path,
            min_panel_size=min_panel_size,
            spot_size=spot_size,
            enable_numba_jit=bool(self.numba_jit_var.get()),
            downsample=downsample,
            additional_colors=additional_colors or None,
            genes=genes,
            features_list=features_list,
            feature_encoding=feature_encoding,
            feature_value_encoding=feature_value_encoding,
            feature_storage=feature_storage,
            also_export_karospace=also_export_karospace,
            feature_manifest_path=feature_manifest_path,
            feature_sidecar_shard_size=feature_sidecar_shard_size,
            feature_sparse_zero_threshold=feature_sparse_zero_threshold,
            modalities=modalities,
            neighbor_stats_groupby=neighbor_stats_groupby,
            neighbor_stats_permutations=neighbor_permutations,
            neighbor_stats_seed=neighbor_seed,
            statistics_additional_annotations=statistics_additional_annotations,
            statistics_modalities=statistics_modalities,
            statistics_contrast_categories=statistics_contrast_categories,
            statistics_counts_layer=statistics_counts_layer,
            statistics_normalization=statistics_normalization,
            statistics_scale_factor=statistics_scale_factor,
            statistics_normalized_layer=statistics_normalized_layer,
            statistics_min_cell_counts=statistics_min_cell_counts,
            statistics_min_feature_counts=statistics_min_feature_counts,
            statistics_n_cpus=statistics_n_cpus,
            wilcoxon=wilcoxon,
            wilcoxon_runtime_limit=wilcoxon_runtime_limit,
            wilcoxon_min_cells_per_group=wilcoxon_min_cells_per_group,
            wilcoxon_min_pct_expressed=wilcoxon_min_pct_expressed,
            wilcoxon_p_adjust_method=wilcoxon_p_adjust_method,
            wilcoxon_padj_cutoff=wilcoxon_padj_cutoff,
            wilcoxon_log2fc_cutoff=wilcoxon_log2fc_cutoff,
            wilcoxon_embed_top_n_per_comparison=wilcoxon_embed_top_n_per_comparison,
            wilcoxon_top_n_per_category=wilcoxon_top_n_per_category,
            pseudobulk=pseudobulk,
            pseudobulk_replicate_annotation=pseudobulk_replicate_annotation,
            pseudobulk_min_cells_per_pseudobulk=pseudobulk_min_cells_per_pseudobulk,
            pseudobulk_min_replicates=pseudobulk_min_replicates,
            pseudobulk_min_pct_expressed=pseudobulk_min_pct_expressed,
            pseudobulk_p_adjust_method=pseudobulk_p_adjust_method,
            pseudobulk_padj_cutoff=pseudobulk_padj_cutoff,
            pseudobulk_log2fc_cutoff=pseudobulk_log2fc_cutoff,
            pseudobulk_deseq2_fit_type=pseudobulk_deseq2_fit_type,
            pseudobulk_embed_top_n_per_comparison=pseudobulk_embed_top_n_per_comparison,
            pathway_gmt=pathway_gmt,
            pathway=pathway,
            pathway_organism=pathway_organism,
            pathway_top_n=pathway_top_n,
            pathway_min_overlap=pathway_min_overlap,
            pathway_gsea_permutations=pathway_gsea_permutations,
            interaction_markers=interaction_markers,
            interaction_markers_top_targets=interaction_top_targets,
            interaction_markers_top_features=interaction_top_features,
            interaction_markers_min_cells=interaction_min_cells,
            interaction_markers_min_neighbors=interaction_min_neighbors,
            section_rotations=section_rotations,
            deconvolutions=deconvolutions,
            feature_correlation_top_n=feature_correlation_top_n,
            spatial_variable_features_n=spatial_variable_features_n,
            scalebar_unit=scalebar_unit,
            section_images=section_images,
            section_images_max_px=section_images_max_px,
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

        def _is_new_export_api(export_func) -> bool:
            try:
                params = inspect.signature(export_func).parameters
            except Exception:
                return False
            required_params = {
                "main_cell_annotation",
                "statistics_additional_annotations",
                "statistics_modalities",
                "wilcoxon",
                "pathway",
                "interaction_markers_top_features",
                "feature_correlation_top_n",
                "spatial_variable_features_n",
            }
            return required_params.issubset(params)

        def _snapshot_karospace_modules() -> dict[str, object]:
            return {name: module for name, module in sys.modules.items() if name == "karospace" or name.startswith("karospace.")}

        def _clear_karospace_modules() -> None:
            for name in list(sys.modules):
                if name == "karospace" or name.startswith("karospace."):
                    del sys.modules[name]

        def _restore_karospace_modules(snapshot: dict[str, object]) -> None:
            _clear_karospace_modules()
            sys.modules.update(snapshot)

        def _candidate_paths() -> list[Path]:
            root = Path(__file__).resolve().parents[2]
            return [
                (root / "KaroSpace").resolve(),
                (root / "spatial-viewer").resolve(),
                (Path.cwd().parent / "KaroSpace").resolve(),
                (Path.cwd().parent / "spatial-viewer").resolve(),
            ]

        def _try_candidate(candidate: Path, previous_modules: dict[str, object]):
            package_dir = candidate / "karospace"
            if not package_dir.exists():
                return None
            candidate_str = str(candidate)
            if candidate_str not in sys.path:
                sys.path.insert(0, candidate_str)
            try:
                _clear_karospace_modules()
                module = importlib.import_module("karospace")
                load_func = module.load_spatial_data
                export_func = module.export_to_html
                if _is_new_export_api(export_func):
                    return load_func, export_func
            except Exception:
                return None
            finally:
                if "karospace" not in sys.modules or not _is_new_export_api(getattr(sys.modules.get("karospace"), "export_to_html", None)):
                    _restore_karospace_modules(previous_modules)
            return None

        try:
            module = importlib.import_module("karospace")
            load_spatial_data = module.load_spatial_data
            export_to_html = module.export_to_html
            if _is_new_export_api(export_to_html):
                return load_spatial_data, export_to_html
            installed_modules = _snapshot_karospace_modules()
            for candidate in _candidate_paths():
                resolved = _try_candidate(candidate, installed_modules)
                if resolved is not None:
                    return resolved
            _restore_karospace_modules(installed_modules)
            return load_spatial_data, export_to_html
        except Exception as exc:
            root_exc = exc
            empty_modules: dict[str, object] = {}
            for candidate in _candidate_paths():
                package_dir = candidate / "karospace"
                if not package_dir.exists():
                    continue
                candidate_str = str(candidate)
                if candidate_str not in sys.path:
                    sys.path.insert(0, candidate_str)
                try:
                    _clear_karospace_modules()
                    module = importlib.import_module("karospace")
                    return module.load_spatial_data, module.export_to_html
                except Exception as inner_exc:
                    _restore_karospace_modules(empty_modules)
                    _raise_dependency_error(inner_exc)
                    continue
            _raise_dependency_error(root_exc)
            raise RuntimeError(
                "Could not import 'karospace'. Install it in this environment before exporting "
                "(for example: pip install -e /path/to/spatial-viewer)."
            ) from root_exc

    @staticmethod
    def _detect_coords_mode(h5ad_path: Path, spatial_key: str = "spatial") -> str:
        if ExportApp._is_zarr_path(h5ad_path):
            return "obsm:spatial"
        ad_mod = _get_anndata()
        adata = None
        try:
            try:
                adata = ad_mod.read_h5ad(h5ad_path, backed="r")
            except Exception:
                adata = ad_mod.read_h5ad(h5ad_path)
            if spatial_key in adata.obsm:
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
    def _build_obs_spatial_h5ad(
        h5ad_path: Path,
        *,
        x_column: str = "centroid_x",
        y_column: str = "centroid_y",
        spatial_key: str = "spatial",
    ) -> Path:
        ad_mod = _get_anndata()
        np_mod = _get_numpy()
        adata = ad_mod.read_h5ad(h5ad_path)
        if x_column not in adata.obs.columns or y_column not in adata.obs.columns:
            raise ValueError(f"Spatial coordinate columns are missing from obs: {x_column}, {y_column}.")
        coords = adata.obs[[x_column, y_column]].to_numpy(dtype=np_mod.float32)
        adata.obsm[spatial_key] = coords
        with tempfile.NamedTemporaryFile(suffix=".h5ad", prefix="karospace_builder_coords_", delete=False) as handle:
            temp_path = Path(handle.name)
        adata.write_h5ad(temp_path)
        return temp_path

    def _resolve_export_input(self, config: BuilderConfig) -> tuple[Path, str, tuple[str, str] | None, Path | None]:
        if config.spatial_columns is not None:
            return config.h5ad_path, config.spatial_key, config.spatial_columns, None
        mode = config.coords_mode
        if mode is None and self._matches_inspected_h5ad(config.h5ad_path) and self._inspected_coords_mode:
            mode = self._inspected_coords_mode
        if mode is None:
            mode = self._detect_coords_mode(config.h5ad_path, config.spatial_key)
        if mode == "obsm:spatial":
            return config.h5ad_path, config.spatial_key, None, None
        if mode == "obs:centroid_x_y":
            temp_h5ad = self._build_obs_spatial_h5ad(
                config.h5ad_path,
                x_column="centroid_x",
                y_column="centroid_y",
                spatial_key=config.spatial_key,
            )
            return temp_h5ad, config.spatial_key, None, temp_h5ad
        raise ValueError(f"Unsupported coordinates mode: {mode}")

    def _set_busy(self, busy: bool) -> None:
        widgets = [
            self.inspect_btn,
            self.numba_jit_check,
            self.wilcoxon_enabled_check,
            self.pseudobulk_enabled_check,
            self.pathway_enabled_check,
            self.interaction_markers_enabled_check,
            self.feature_modality_combo,
            self.save_parameters_btn,
            self.import_parameters_btn,
            self.theme_toggle_btn,
        ]
        for widget in widgets:
            self._configure_widget_state(widget, not busy)
        self.additional_colors_editor.set_enabled(not busy)
        self.modalities_editor.set_enabled(not busy)
        self.statistics_annotations_editor.set_enabled(not busy)
        self.statistics_modalities_editor.set_enabled(not busy)
        self.section_metadata_editor.set_enabled(not busy)
        self.section_metadata_extra_editor.set_enabled(not busy)
        self.manual_genes_editor.set_enabled(not busy)

        if busy:
            self._set_progress(0, "Queued")
        else:
            if self.status_var.get().startswith("Export running..."):
                self.status_var.set("Ready")
        self._refresh_input_gate()

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

        self._cancel_requested.clear()
        self._set_busy(True)
        self._log(f"Starting export: {config.h5ad_path} -> {config.output_html_path}")
        self._log(
            "Export options: "
            f"coords={config.coords_mode or 'auto'}, "
            f"spatial_key={config.spatial_key}, "
            f"spatial_columns={config.spatial_columns or 'obsm'}, "
            f"section_key={config.section_groupby or '(single section)'}, "
            f"main_cell_annotation={config.initial_color}, "
            f"feature_storage={config.feature_storage}, "
            f"downsample={config.downsample if config.downsample is not None else 'all'}."
        )
        self._log(
            "Feature settings: "
            f"features={len(config.genes or [])}, "
            f"feature_encoding={config.feature_encoding}."
        )
        self._log(
            "Analytics settings: "
            f"pseudobulk={config.pseudobulk or 'None'}, "
            f"statistics_annotations={len(config.statistics_additional_annotations or [])}, "
            f"statistics_modalities={config.statistics_modalities or 'default'}, "
            f"wilcoxon={config.wilcoxon}, "
            f"pathway={config.pathway or 'None'}, "
            f"neighbor_annotations={len(config.neighbor_stats_groupby or []) if config.neighbor_stats_groupby is not None else 'default'}, "
            f"neighbor_permutations={config.neighbor_stats_permutations if config.neighbor_stats_permutations is not None else 'auto'}, "
            f"interaction_markers={config.interaction_markers or 'None'}."
        )
        self._log(
            "Runtime mode: "
            f"numba_jit={'on' if config.enable_numba_jit else 'off (safe mode)'}."
        )
        serve_after_export = bool(self.serve_var.get())
        thread = threading.Thread(target=self._run_export, args=(config, serve_after_export), daemon=True)
        self._export_thread = thread
        thread.start()

    def _raise_if_cancelled(self) -> None:
        if self._cancel_requested.is_set():
            raise _ExportCancelled()

    def _run_export(self, config: BuilderConfig, serve_after_export: bool) -> None:
        stdout_tee = _EventLogTee(sys.stdout, self._queue)
        stderr_tee = _EventLogTee(sys.stderr, self._queue)
        with contextlib.redirect_stdout(stdout_tee), contextlib.redirect_stderr(stderr_tee):
            try:
                self._run_export_body(config, serve_after_export)
            finally:
                stdout_tee.flush()
                stderr_tee.flush()

    def _run_export_body(self, config: BuilderConfig, serve_after_export: bool) -> None:
        temp_h5ad: Path | None = None
        total_started = time.perf_counter()

        def emit_progress(percent: int, stage: str, detail: str | None = None) -> None:
            self._queue.put(("progress", (percent, stage, detail)))

        try:
            self._raise_if_cancelled()
            emit_progress(5, "Importing API", "Resolving karospace export functions.")
            load_spatial_data, export_to_html = self._import_karospace_api(enable_numba_jit=config.enable_numba_jit)

            self._raise_if_cancelled()
            emit_progress(12, "Preparing input", "Resolving coordinates mode and source data.")
            input_path, spatial_key, spatial_columns, temp_h5ad = self._resolve_export_input(config)
            if temp_h5ad is not None:
                self._queue.put(("log", f"Converted obs coordinate columns to temporary obsm['{spatial_key}']."))

            self._raise_if_cancelled()
            emit_progress(
                25,
                "Loading spatial data",
                f"Reading {input_path} with section_key='{config.section_groupby}' and spatial_key='{spatial_key}'.",
            )
            load_started = time.perf_counter()
            import inspect

            load_params = inspect.signature(load_spatial_data).parameters
            if "section_key" in load_params:
                load_kwargs: dict[str, object] = {
                    "section_key": config.section_groupby,
                    "spatial_key": spatial_key,
                }
                if config.spatialdata_table:
                    load_kwargs["spatialdata_table"] = config.spatialdata_table
                if spatial_columns is not None:
                    load_kwargs["spatial_columns"] = spatial_columns
                if config.section_order is not None:
                    load_kwargs["section_order"] = config.section_order
                if config.section_metadata is not None:
                    load_kwargs["section_metadata"] = config.section_metadata
                if config.section_metadata_extra is not None:
                    load_kwargs["section_metadata_extra"] = config.section_metadata_extra
                if config.metadata_value_order is not None:
                    load_kwargs["metadata_value_order"] = config.metadata_value_order
                if config.metadata_max_columns is not None:
                    load_kwargs["metadata_max_columns"] = config.metadata_max_columns
                dataset = load_spatial_data(str(input_path), **load_kwargs)
            else:
                old_input_path = input_path
                if spatial_columns is not None:
                    old_input_path = self._build_obs_spatial_h5ad(
                        config.h5ad_path,
                        x_column=spatial_columns[0],
                        y_column=spatial_columns[1],
                        spatial_key=spatial_key,
                    )
                    temp_h5ad = old_input_path
                dataset = load_spatial_data(
                    str(old_input_path),
                    groupby=config.section_groupby,
                    spatial_key=spatial_key,
                )
            load_elapsed = time.perf_counter() - load_started
            self._raise_if_cancelled()
            emit_progress(
                55,
                "Validating analytics",
                f"Dataset loaded: sections={int(dataset.n_sections)}, cells={int(dataset.n_cells)} in {load_elapsed:.1f}s.",
            )

            neighbor_groupby = config.neighbor_stats_groupby
            interaction_groupby = config.statistics_additional_annotations if config.interaction_markers else None
            neighbor_permutations = config.neighbor_stats_permutations
            emit_progress(
                68,
                "Preparing output",
                "Resolved analytics arguments: "
                f"statistics_additional={len(config.statistics_additional_annotations or [])}, "
                f"neighbor={len(neighbor_groupby or []) if neighbor_groupby is not None else 'auto/default'}, "
                f"interaction={len(interaction_groupby or []) if config.interaction_markers else 0}, "
                f"neighbor_permutations={neighbor_permutations if neighbor_permutations is not None else 'auto'}.",
            )

            output_html = config.output_html_path
            self._raise_if_cancelled()
            emit_progress(76, "Writing viewer", f"Exporting HTML viewer to {output_html}.")
            export_started = time.perf_counter()
            export_params = inspect.signature(export_to_html).parameters
            if "main_cell_annotation" in export_params:
                export_kwargs: dict[str, object] = {
                    "output_path": str(output_html),
                    "main_cell_annotation": config.initial_color,
                    "title": config.title,
                    "min_panel_size": config.min_panel_size,
                    "spot_size": config.spot_size,
                    "downsample": config.downsample,
                    "outline_by": config.outline_by,
                    "metadata_labels": config.metadata_labels,
                    "viewer_info_html": config.viewer_info_html,
                    "tutorial": config.tutorial,
                    "embed_reproducibility_info": config.embed_reproducibility_info,
                    "source_input_path": config.source_input_path,
                    "cell_annotations": config.additional_colors,
                    "features": config.genes,
                    "features_list": config.features_list,
                    "feature_encoding": config.feature_encoding,
                    "feature_value_encoding": config.feature_value_encoding,
                    "feature_storage": config.feature_storage,
                    "feature_manifest_path": config.feature_manifest_path,
                    "feature_sidecar_shard_size": config.feature_sidecar_shard_size,
                    "feature_sparse_zero_threshold": config.feature_sparse_zero_threshold,
                    "statistics_additional_annotations": config.statistics_additional_annotations,
                    "statistics_modalities": config.statistics_modalities,
                    "statistics_contrast_categories": config.statistics_contrast_categories,
                    "statistics_counts_layer": config.statistics_counts_layer,
                    "statistics_normalization": config.statistics_normalization,
                    "statistics_scale_factor": config.statistics_scale_factor,
                    "statistics_normalized_layer": config.statistics_normalized_layer,
                    "statistics_min_cell_counts": config.statistics_min_cell_counts,
                    "statistics_min_feature_counts": config.statistics_min_feature_counts,
                    "statistics_n_cpus": config.statistics_n_cpus,
                    "wilcoxon": config.wilcoxon,
                    "wilcoxon_runtime_limit": config.wilcoxon_runtime_limit,
                    "wilcoxon_min_cells_per_group": config.wilcoxon_min_cells_per_group,
                    "wilcoxon_min_pct_expressed": config.wilcoxon_min_pct_expressed,
                    "wilcoxon_p_adjust_method": config.wilcoxon_p_adjust_method,
                    "wilcoxon_padj_cutoff": config.wilcoxon_padj_cutoff,
                    "wilcoxon_log2fc_cutoff": config.wilcoxon_log2fc_cutoff,
                    "wilcoxon_embed_top_n_per_comparison": config.wilcoxon_embed_top_n_per_comparison,
                    "wilcoxon_top_n_per_category": config.wilcoxon_top_n_per_category,
                    "pseudobulk": config.pseudobulk,
                    "pseudobulk_replicate_annotation": config.pseudobulk_replicate_annotation,
                    "pseudobulk_min_cells_per_pseudobulk": config.pseudobulk_min_cells_per_pseudobulk,
                    "pseudobulk_min_replicates": config.pseudobulk_min_replicates,
                    "pseudobulk_min_pct_expressed": config.pseudobulk_min_pct_expressed,
                    "pseudobulk_p_adjust_method": config.pseudobulk_p_adjust_method,
                    "pseudobulk_padj_cutoff": config.pseudobulk_padj_cutoff,
                    "pseudobulk_log2fc_cutoff": config.pseudobulk_log2fc_cutoff,
                    "pseudobulk_deseq2_fit_type": config.pseudobulk_deseq2_fit_type,
                    "pseudobulk_embed_top_n_per_comparison": config.pseudobulk_embed_top_n_per_comparison,
                    "pathway": config.pathway,
                    "pathway_gmt": config.pathway_gmt,
                    "pathway_organism": config.pathway_organism,
                    "pathway_top_n": config.pathway_top_n,
                    "pathway_min_overlap": config.pathway_min_overlap,
                    "pathway_gsea_permutations": config.pathway_gsea_permutations,
                    "neighbor_stats_annotations": neighbor_groupby,
                    "neighbor_stats_permutations": neighbor_permutations,
                    "neighbor_stats_seed": config.neighbor_stats_seed,
                    "interaction_markers": config.interaction_markers,
                    "interaction_markers_top_targets": config.interaction_markers_top_targets,
                    "interaction_markers_top_features": config.interaction_markers_top_features,
                    "interaction_markers_min_cells": config.interaction_markers_min_cells,
                    "interaction_markers_min_neighbors": config.interaction_markers_min_neighbors,
                    "section_rotations": config.section_rotations,
                    "deconvolutions": config.deconvolutions,
                    "feature_correlation_top_n": config.feature_correlation_top_n,
                    "spatial_variable_features_n": config.spatial_variable_features_n,
                    "scalebar_unit": config.scalebar_unit,
                    "modalities": config.modalities,
                    "section_images": config.section_images,
                    "section_images_max_px": config.section_images_max_px,
                }
                export_kwargs = {key: value for key, value in export_kwargs.items() if key in export_params}
                output_value = export_to_html(dataset, **export_kwargs)
            else:
                raise RuntimeError(
                    "KaroSpaceBuilder requires the current Karospace export_to_html API "
                    "(expected main_cell_annotation/features/pseudobulk/pathway arguments)."
                )
            output_html_path = Path(output_value).expanduser()
            export_elapsed = time.perf_counter() - export_started
            self._raise_if_cancelled()
            emit_progress(95, "Finalizing", f"Viewer bundle created in {export_elapsed:.1f}s: {output_html_path}")

            if getattr(config, "also_export_karospace", False) and output_html_path.suffix.lower() == ".html":
                self._raise_if_cancelled()
                emit_progress(97, "Packaging", "Creating .karospace package from the sidecar bundle...")
                try:
                    import karospace as _karospace

                    package_str = _karospace.package_sidecar_viewer(str(output_html_path))
                    package_path = Path(package_str)
                    emit_progress(98, "Packaging", f"Created .karospace package: {package_path}")
                    self._queue.put(("log", f"Loader for the package: {package_path.with_suffix('.loader.html')}"))
                except Exception as exc:  # noqa: BLE001
                    self._queue.put(
                        ("log", f"Warning: HTML viewer was exported but the .karospace package failed: {exc}")
                    )

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
                self._raise_if_cancelled()
                preview_path = output_html_path
                if output_html_path.suffix.lower() == ".karospace":
                    loader_path = output_html_path.with_suffix(".loader.html")
                    if loader_path.exists():
                        preview_path = loader_path
                self._queue.put(("log", "Starting preview server for the exported viewer."))
                self._queue.put(("start_server", preview_path))
        except _ExportCancelled:
            self._queue.put(("canceled", None))
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
                self._cancel_requested.clear()
                result = payload
                assert isinstance(result, AppResult)
                self._last_outdir = result.outdir
                self._last_output_html = result.output_html
                self._log(
                    f"Export complete. sections={result.n_sections}, cells={result.n_cells}, html={result.output_html}"
                )
                self.status_var.set("Export complete")
            elif kind == "canceled":
                self._set_busy(False)
                self._cancel_requested.clear()
                self._log("Export canceled.")
                self.status_var.set("Canceled")
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
                self._cancel_requested.clear()
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
        output_text = self.output_html_var.get().strip() if hasattr(self, "output_html_var") else ""
        if output_text:
            base = Path(self.outdir_var.get().strip()).expanduser() if self.outdir_var.get().strip() else Path.cwd()
            configured = Path(output_text).expanduser()
            configured = configured if configured.is_absolute() else base / configured
            if configured.exists():
                return configured
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
        if self._last_outdir is not None:
            path = self._last_outdir
        elif self.output_html_var.get().strip():
            base = Path(self.outdir_var.get().strip()).expanduser() if self.outdir_var.get().strip() else Path.cwd()
            configured = Path(self.output_html_var.get().strip()).expanduser()
            path = (configured if configured.is_absolute() else base / configured).parent
        else:
            path = Path(self.outdir_var.get().strip()).expanduser() if self.outdir_var.get().strip() else None
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
        self._stop_runtime_chip_animation()
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
