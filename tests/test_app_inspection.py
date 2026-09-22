from __future__ import annotations

from dataclasses import fields
import inspect
import json

import numpy as np
import pandas as pd

from karospace_export.app import BuilderConfig, ExportApp, SearchableListEditor, SectionOrderEditor
from tests.test_smoke_export import _make_test_adata


def test_inspect_loaded_dataset_formats_karospace_inspect_output(tmp_path):
    h5ad_path = tmp_path / "toy.h5ad"
    adata = _make_test_adata()

    summary = ExportApp._inspect_loaded_dataset(
        adata,
        path=h5ad_path,
        spatialdata_table="table_a",
    )
    output = ExportApp._format_dataset_inspection_output(summary)
    parsed = json.loads(output)

    assert parsed["path"] == str(h5ad_path)
    assert parsed["spatialdata_table"] == "table_a"
    assert parsed["n_cells"] == 80
    assert parsed["n_genes"] == 24
    assert any(column["name"] == "cell_type" for column in parsed["metadata"])
    assert any(modality["name"] == "rna" for modality in parsed["feature_modalities"])


def test_builder_exposes_statistics_modalities_export_argument():
    config_fields = {field.name for field in fields(BuilderConfig)}
    run_export_source = inspect.getsource(ExportApp._run_export_body)
    layout_source = inspect.getsource(ExportApp._build_layout)
    inspect_source = inspect.getsource(ExportApp._inspect_h5ad)
    parse_source = inspect.getsource(ExportApp._parse_config)

    assert "statistics_modalities" in config_fields
    assert '"statistics_modalities": config.statistics_modalities' in run_export_source
    assert "self.statistics_modalities_editor = SearchableListEditor" in layout_source
    assert "self.pseudobulk_modalities_editor" not in layout_source
    assert "self._sync_statistics_modality_choices()" in inspect_source
    assert "self.statistics_modalities_editor.get_items()" in parse_source


def test_parse_modalities_option():
    assert ExportApp._parse_modalities_option("", "Modalities") is None
    assert ExportApp._parse_modalities_option("None", "Modalities") is None
    assert ExportApp._parse_modalities_option("all", "Modalities") == "all"
    assert ExportApp._parse_modalities_option("rna, protein", "Modalities") == ["rna", "protein"]


def test_connections_tab_hosts_discovery_layout_and_features_cleanup():
    source = inspect.getsource(ExportApp._build_layout)

    assert source.count('"Discovery panels"') == 1
    assert "notebook.add(\"Statistics\")" in source
    assert "notebook.add(\"Pseudobulk\")" not in source
    assert "notebook.add(\"Neighborhoods\")" in source
    assert "notebook.add(\"Connections\")" in source
    assert "notebook.add(\"Viewer\")" in source
    assert "notebook.add(\"Sections\")" in source
    assert "notebook.add(\"Overlays\")" not in source
    assert source.index('notebook.add("Neighborhoods")') < source.index('notebook.add("Connections")')
    assert source.index('notebook.add("Connections")') < source.index('notebook.add("Viewer")')
    assert "discovery_inner = self._make_sub_frame(connections_tab)" in source
    assert "discovery_inner = self._make_sub_frame(genes_inner)" not in source
    assert '"Category means"' not in source
    assert "self.modalities_editor = SearchableListEditor" in source
    assert 'label="Modality"' in source
    assert "on_change=self._on_embed_modalities_changed" in source
    assert "self.feature_modality_combo = ctk.CTkComboBox" in source
    assert 'self._body_label(modality_row, "Modality")' in source
    assert "self.modalities_entry = ctk.CTkEntry" not in source
    assert "self._build_parameter_state_buttons(controls, 3)" in source
    assert "self._build_spatialdata_table_selector(controls, 4)" in source
    assert "selection_mode_check" not in source
    assert "feature_manifest_entry" not in source
    assert "stack_controls=True" in source
    assert 'self._option_row(connections_tab, 5, "Deconvolutions JSON", widget=deconv_row)' not in source
    assert '"Feature Storage"' in source
    assert '"Storage settings apply to all selected modalities."' in source
    assert "self.statistics_counts_layer_combo = combo" in source
    assert "self.statistics_normalized_layer_combo = combo" in source
    assert 'self._option_row(viewer_tab, 11, "Deconvolutions JSON", widget=deconv_row)' not in source
    assert 'self._option_row(\n            self.overlays_content,\n            9,\n            "Deconvolutions JSON"' in source
    assert 'self._option_row(self.overlays_content, 4, "Scalebar unit", widget=scalebar_row)' in source
    assert "Checked maps to --pseudobulk auto." not in source
    assert "--pseudobulk auto analyzes the main annotation" not in source
    assert "Also supports cell/gene count filters below" not in source
    assert "Maps to --pseudobulk-simple-constrast-categories" not in source
    assert "self._choose_pathway_gmt_files" in source
    assert "self.wilcoxon_enabled_check = ctk.CTkCheckBox" in source
    assert "self.pseudobulk_enabled_check = ctk.CTkCheckBox" in source
    assert "self.pseudobulk_replicate_combo = ctk.CTkComboBox" in source
    assert '"Replicate annotation"' in source
    assert "self.pseudobulk_replicate_entry" not in source
    assert source.index('"Pseudobulk fit"') < source.index('"Auto-embedded DE genes"')
    assert "self.pathway_enabled_check = ctk.CTkCheckBox" in source
    assert "self.interaction_markers_enabled_check = ctk.CTkCheckBox" in source
    assert "self._sync_analysis_controls" in source
    assert "row_collector=self._wilcoxon_detail_rows" in source
    assert "row_collector=self._pseudobulk_detail_rows" in source
    assert "row_collector=self._pathway_detail_rows" in source
    assert "row_collector=self._interaction_detail_rows" in source
    assert "Values from the selected Section key" not in source
    assert "Values come from the selected Section key" not in source
    assert "Maps to --main-cell-annotation" not in source
    assert "Maps to --outlineby" not in source
    assert "Maps to CLI --output" not in source
    assert "Search var_names and build the exported feature" not in source


def test_parameter_state_json_helpers_are_present_and_stable():
    source = inspect.getsource(ExportApp)

    assert '_PARAMETER_STATE_SCHEMA = "karospace_builder_parameters"' in source
    assert "def _build_parameter_state_buttons" in source
    assert "Save JSON" in source
    assert "Import JSON" in source
    assert "def _collect_parameter_state" in source
    assert "def _validate_parameter_state_format" in source
    assert "def _apply_parameter_state" in source
    assert "def _save_parameter_state_json" in source
    assert "def _import_parameter_state_json" in source
    assert "self._inspect_imported_parameter_input()" in source
    assert "def _inspect_imported_parameter_input" in source
    assert "self._save_active_feature_selection()" in source
    assert "feature_items_by_modality" in source
    assert ExportApp._coerce_bool_state_value("true") is True
    assert ExportApp._coerce_bool_state_value("off") is False
    assert ExportApp._json_string_list(["A", "A", "", "B"]) == ["A", "B"]
    validator = ExportApp.__new__(ExportApp)
    valid_state = {
        "schema": "karospace_builder_parameters",
        "version": 1,
        "variables": {},
        "lists": {},
        "feature_items_by_modality": {},
        "active_feature_modality": "rna",
    }
    assert ExportApp._validate_parameter_state_format(validator, valid_state) is valid_state
    for invalid_state in (
        [],
        {"version": 1, "variables": {}},
        {"schema": "wrong", "version": 1, "variables": {}},
        {"schema": "karospace_builder_parameters", "version": 2, "variables": {}},
        {"schema": "karospace_builder_parameters", "version": 1, "variables": []},
    ):
        try:
            ExportApp._validate_parameter_state_format(validator, invalid_state)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Invalid parameter state accepted: {invalid_state!r}")


def test_runtime_chip_running_state_animates_ellipsis():
    sync_source = inspect.getsource(ExportApp._sync_runtime_chip)
    animate_source = inspect.getsource(ExportApp._animate_runtime_chip)
    close_source = inspect.getsource(ExportApp._on_close)

    assert "self._start_runtime_chip_animation()" in sync_source
    assert "Running{'.' * count}" in animate_source
    assert "dot_counts = (0, 1, 2, 3, 2, 1)" in animate_source
    assert "self.after(260, self._animate_runtime_chip)" in animate_source
    assert "self._stop_runtime_chip_animation()" in close_source


def test_embed_modality_selection_controls_feature_dropdown_and_export_modalities():
    source = inspect.getsource(ExportApp)
    inspect_source = inspect.getsource(ExportApp._inspect_h5ad)
    parse_source = inspect.getsource(ExportApp._parse_config)
    selected_source = inspect.getsource(ExportApp._selected_features_by_modality)

    assert "def _selected_embed_modalities" in source
    assert "def _feature_modality_choices" in source
    assert "def _sync_feature_modality_choices" in source
    assert "def _sync_statistics_modality_choices" in source
    assert "self.modalities_editor.set_choices(choices)" in source
    assert "self.modalities_editor.set_items(existing)" in source
    assert "self._sync_statistics_modality_choices()" in source
    assert "selected_modality_choices = set(self._feature_modality_choices())" in inspect_source
    assert "for modality in self._feature_modality_choices()" in selected_source
    assert "embed_modalities = self._selected_embed_modalities()" in parse_source
    assert "modalities = embed_modalities" in parse_source
    assert 'self.pseudobulk_replicate_combo.configure(values=[""] + obs_cols)' in inspect_source


def test_feature_names_are_detected_by_modality():
    adata = _make_test_adata()
    adata.obsm["protein"] = np.ones((adata.n_obs, 3))
    adata.uns["protein_var"] = pd.DataFrame({"protein": ["CD3", "CD4", "CD8"]})

    names_by_modality = ExportApp._feature_names_by_modality_from_adata(adata)

    assert names_by_modality["rna"][:2] == ["Gene000", "Gene001"]
    assert names_by_modality["protein"] == ["CD3", "CD4", "CD8"]
    assert "spatial" not in names_by_modality


def test_layer_names_are_detected_for_statistics_dropdowns():
    adata = _make_test_adata()
    adata.layers["counts"] = adata.X.copy()
    adata.layers["normalized"] = adata.X.copy()

    assert ExportApp._layer_keys_from_adata(adata) == ["counts", "normalized"]


def test_section_values_are_extracted_in_observed_order():
    adata = type("FakeAdata", (), {})()
    adata.obs = pd.DataFrame({"sample_id": ["S2", "S1", "S2", "S3", None, "S1"]})

    assert ExportApp._section_values_from_adata(adata, "sample_id") == ["S2", "S1", "S3"]


def test_section_key_columns_require_fewer_than_500_unique_values():
    adata = type("FakeAdata", (), {})()
    adata.obs = pd.DataFrame(
        {
            "sample_id": [f"S{i % 3}" for i in range(500)],
            "exactly_500": [f"v{i}" for i in range(500)],
            "under_500": [f"u{i % 499}" for i in range(500)],
            "all_missing": [None for _ in range(500)],
        }
    )

    assert ExportApp._eligible_section_key_columns(
        adata,
        ["sample_id", "exactly_500", "under_500", "all_missing"],
    ) == ["sample_id", "under_500"]


def test_spatial_dropdown_choices_are_derived_from_inspected_data():
    adata = _make_test_adata()
    adata.obs["centroid_x"] = [float(i) for i in range(adata.n_obs)]
    adata.obs["centroid_y"] = [float(i * 2) for i in range(adata.n_obs)]
    adata.obs["label"] = "A"
    adata.obsm["protein"] = adata.obsm["spatial"]

    assert ExportApp._obsm_keys_from_adata(adata) == ["spatial", "protein"]
    del adata.obsm["spatial"]
    assert ExportApp._obsm_keys_from_adata(adata) == ["protein"]
    del adata.obsm["protein"]
    assert ExportApp._obsm_keys_from_adata(adata) == ["spatial"]
    assert ExportApp._spatial_obs_columns_from_adata(
        adata,
        ["cell_type", "centroid_x", "centroid_y", "label"],
    ) == ["centroid_x", "centroid_y"]


def test_section_order_uses_reorderable_editor():
    layout_source = inspect.getsource(ExportApp._build_layout)
    editor_source = inspect.getsource(SectionOrderEditor)

    assert "self.section_order_editor = SectionOrderEditor" in layout_source
    assert "section_order_entry = ctk.CTkEntry" not in layout_source
    assert "<B1-Motion>" in editor_source


def test_spatial_coordinates_use_dropdown_controls():
    source = inspect.getsource(ExportApp._build_layout)
    inspect_source = inspect.getsource(ExportApp._inspect_h5ad)

    assert "self.spatial_key_combo = ctk.CTkComboBox" in source
    assert "self.spatial_x_combo = ctk.CTkComboBox" in source
    assert "self.spatial_y_combo = ctk.CTkComboBox" in source
    assert "self.spatial_key_entry = ctk.CTkEntry" not in source
    assert "self.spatial_x_entry = ctk.CTkEntry" not in source
    assert "self.spatial_y_entry = ctk.CTkEntry" not in source
    assert "self.spatial_key_combo.configure(values=obsm_keys)" in inspect_source
    assert 'self.spatial_x_combo.configure(values=[""] + spatial_obs_cols)' in inspect_source
    assert 'self.spatial_y_combo.configure(values=[""] + spatial_obs_cols)' in inspect_source


def test_inspection_uses_filtered_section_and_annotation_choices():
    source = inspect.getsource(ExportApp._inspect_h5ad)
    parse_source = inspect.getsource(ExportApp._parse_config)

    assert "section_key_cols = self._eligible_section_key_columns(adata, obs_cols)" in source
    assert 'self.groupby_combo.configure(values=[""] + section_key_cols)' in source
    assert "if current_section_groupby and current_section_groupby not in section_key_col_set" in source
    assert "self.color_combo.configure(values=section_key_cols)" in source
    assert 'self.outline_combo.configure(values=[""] + section_key_cols)' in source
    assert "self.additional_colors_editor.set_choices(section_key_cols)" in source
    assert "self.section_metadata_editor.set_choices(section_key_cols)" in source
    assert "self.section_metadata_extra_editor.set_choices(section_key_cols)" in source
    assert "Section key must have fewer than 500 unique values." in parse_source
    assert "Section key is required." not in parse_source
    assert "if section_groupby and section_groupby not in eligible_annotation_cols" in parse_source
    assert "Main cell annotation must have fewer than 500 unique values." in parse_source
    assert "Outline by must have fewer than 500 unique values." in parse_source
    assert "cell_annotations" in parse_source
    assert "section_metadata" in parse_source
    assert "section_metadata_extra" in parse_source


def test_required_metadata_keeps_section_key_out_of_section_metadata():
    editor_source = inspect.getsource(SearchableListEditor)
    layout_source = inspect.getsource(ExportApp._build_layout)
    required_source = inspect.getsource(ExportApp._required_section_metadata_items)
    parse_source = inspect.getsource(ExportApp._parse_config)

    assert "def set_required_items" in editor_source
    assert "if str(self.listbox.get(idx)) in required" in editor_source
    assert "self.outline_by_var.trace_add" in layout_source
    assert "self._required_section_metadata_items()" in parse_source
    assert "section_groupby_var" not in required_source
    assert "statistics_additional_annotations" in parse_source
    assert "def _sync_required_cell_annotations" in inspect.getsource(ExportApp)
    assert "self.initial_color_var.trace_add" in layout_source
    assert "additional_colors = self._merge_unique([initial_color]" in parse_source


def test_analysis_checkboxes_disable_hidden_argument_validation():
    parse_source = inspect.getsource(ExportApp._parse_config)
    variables_source = inspect.getsource(ExportApp._build_variables)
    preset_source = inspect.getsource(ExportApp._apply_preset)
    sync_source = inspect.getsource(ExportApp._sync_analysis_controls)

    assert "if wilcoxon_enabled:" in parse_source
    assert "else:\n            wilcoxon_runtime_limit = \"00:30:00\"" in parse_source
    assert "if pseudobulk_enabled:" in parse_source
    assert "else:\n            pseudobulk_replicate_annotation = None" in parse_source
    assert "if pathway_enabled:" in parse_source
    assert "else:\n            pathway_gmt = None" in parse_source
    assert "if interaction_enabled:" in parse_source
    assert "else:\n            interaction_top_targets = 5" in parse_source
    assert 'self.pathway_enabled_var = tk.BooleanVar(value=False)' in variables_source
    assert 'self.pathway_enabled_var.set(False)' in preset_source
    assert '("_pathway_detail_rows", enabled("pathway_enabled_var", False))' in sync_source
    assert 'self.neighbor_permutations_var = tk.StringVar(value="20")' in variables_source
    assert 'self.neighbor_permutations_var.set("20")' in preset_source


def test_export_api_probe_requires_current_karospace_arguments():
    source = inspect.getsource(ExportApp._import_karospace_api)

    assert '"statistics_additional_annotations"' in source
    assert '"statistics_modalities"' in source
    assert '"wilcoxon"' in source
    assert '"pathway"' in source
    assert '"interaction_markers_top_features"' in source
    assert '"feature_correlation_top_n"' in source
    assert '"spatial_variable_features_n"' in source
    assert "required_params.issubset(params)" in source
