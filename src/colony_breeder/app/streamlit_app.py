"""Streamlit upload page for loading colony structure only."""

import hashlib
from typing import Dict, Optional, Sequence, Tuple

import streamlit as st

from colony_breeder.app.streamlit_logic import (
    colony_sex_warnings,
    default_column_selections,
    default_sex_value_mapping,
    drop_stale_genes,
    format_load_result,
    gene_column_options,
    input_signature,
    preview_rows,
    sex_mapping_check,
    sex_value_label,
)
from colony_breeder.app.ui_helpers import (
    NONE_COLUMN,
    UploadLoadResult,
    load_from_upload,
    summarize_distinct_values,
)
from colony_breeder.app.upload import SheetChoiceRequired, UploadResult, read_upload

_OPTIONAL_COLUMNS = {"cage", "mother_id", "father_id"}
_MAPPING_WIDGETS = ("animal_id", "sex", "cage", "mother_id", "father_id")
_SEX_OPTIONS = ("M", "F", "unknown")
_STRUCTURE_NOTICE = (
    "Genotypes have NOT been loaded or validated. Ranking is unavailable "
    "until step 7c."
)


def _reset_upload_state() -> None:
    for key in (
        "mapping_selection",
        "sex_value_mapping",
        "gene_selection",
        "loaded_result",
        "loaded_signature",
        "header_context",
    ):
        st.session_state.pop(key, None)
    for key in _MAPPING_WIDGETS:
        st.session_state.pop(f"mapping_widget_{key}", None)
    st.session_state.pop("gene_selection_widget", None)
    for key in list(st.session_state):
        if key.startswith("sex_widget_"):
            st.session_state.pop(key, None)


def _mapping_selectboxes(
    fieldnames: Sequence[str],
    previous: Optional[Dict[str, Optional[str]]],
) -> Dict[str, Optional[str]]:
    defaults = default_column_selections(fieldnames, previous)
    mapping: Dict[str, Optional[str]] = {}
    st.subheader("Column mapping")
    for key in _MAPPING_WIDGETS:
        options = list(fieldnames)
        if key in _OPTIONAL_COLUMNS:
            options.append(NONE_COLUMN)
        label = key.replace("_", " ").title()
        widget_key = f"mapping_widget_{key}"
        default = defaults.get(key)
        if default not in options:
            default = NONE_COLUMN if key in _OPTIONAL_COLUMNS else None
        if (
            widget_key not in st.session_state
            or st.session_state[widget_key] not in options
        ):
            st.session_state[widget_key] = default
        mapping[key] = st.selectbox(
            label,
            options=options,
            placeholder=f"Choose {label.lower()}",
            key=widget_key,
        )
    st.session_state["mapping_selection"] = mapping
    return mapping


def _sex_mapping_table(
    upload: UploadResult,
    sex_column: Optional[str],
) -> Tuple[Dict[str, str], Optional[str]]:
    if sex_column is None or sex_column == NONE_COLUMN:
        return {}, "Load disabled: choose a sex column to map observed values."

    check = sex_mapping_check(upload.rows, sex_column)
    if not check.displayable:
        return {}, "Load disabled: " + (
            check.warning or "the selected sex column has too many distinct values."
        )

    defaults = default_sex_value_mapping(
        upload.rows,
        sex_column,
        st.session_state.get("sex_value_mapping"),
    )
    if not defaults:
        return {}, "Load disabled: the selected sex column has no values to map."

    mapping = {}
    st.subheader("Map sex values")
    for raw_value, default in defaults.items():
        widget_key = f"sex_widget_{raw_value}"
        default_choice = "unknown" if default in ("", "unknown") else default
        if (
            widget_key not in st.session_state
            or st.session_state[widget_key] not in _SEX_OPTIONS
        ):
            st.session_state[widget_key] = default_choice
        choice = st.selectbox(
            sex_value_label(raw_value),
            options=_SEX_OPTIONS,
            key=widget_key,
        )
        mapping[raw_value] = "" if choice == "unknown" else choice
    st.session_state["sex_value_mapping"] = mapping
    return mapping, None


def _display_load_result(result: UploadLoadResult) -> None:
    view = format_load_result(result)
    st.subheader("Load result")
    st.info(_STRUCTURE_NOTICE)
    st.write(f"Loaded animals: {view.animal_count}")
    if view.error_message:
        st.error(view.error_message)
    st.write(f"Animals with unknown sex: {view.unknown_sex_count}")
    if view.unknown_sex_ids:
        st.write("Unknown-sex animal IDs: " + ", ".join(view.unknown_sex_ids))
    for group, messages in view.grouped_warnings:
        st.markdown(f"**{group}**")
        for message in messages:
            st.warning(message)
    for warning in colony_sex_warnings(result.animals):
        st.warning(warning)
    for gene, values in view.gene_value_summaries:
        st.write(f"{gene}: {', '.join(values)}")


def _choose_workbook_sheet() -> None:
    """Show a sheet picker after discovery and persist its selected sheet."""
    sheet_names = st.session_state.get("available_sheets", ())
    if not sheet_names:
        return

    current_sheet = st.session_state.get("selected_sheet")
    if current_sheet not in sheet_names:
        current_sheet = sheet_names[0]
        st.session_state["sheet_picker"] = current_sheet
    elif "sheet_picker" not in st.session_state:
        st.session_state["sheet_picker"] = current_sheet

    selected_sheet = st.selectbox(
        "Choose a worksheet",
        options=sheet_names,
        key="sheet_picker",
    )
    if selected_sheet != st.session_state.get("selected_sheet"):
        _reset_upload_state()
        st.session_state["selected_sheet"] = selected_sheet


def main() -> None:
    st.set_page_config(page_title="Colony Breeder Upload")
    st.title("Load colony data")
    st.caption(
        "Uploads are processed in memory and are not stored. For sensitive data, "
        "run the app locally."
    )

    uploaded_file = st.file_uploader(
        "Upload a colony CSV or XLSX file",
        type=["csv", "xlsx"],
        key="colony_upload",
    )
    if uploaded_file is None:
        return

    content = uploaded_file.getvalue()
    file_name = uploaded_file.name
    content_hash = hashlib.sha256(content).hexdigest()
    file_identity = (file_name, content_hash)
    if st.session_state.get("file_identity") != file_identity:
        _reset_upload_state()
        for key in (
            "header_row_widget",
            "sheet_picker",
            "available_sheets",
            "selected_sheet",
        ):
            st.session_state.pop(key, None)
        st.session_state["file_identity"] = file_identity
        st.session_state["header_row"] = 1

    selected_sheet = st.session_state.get("selected_sheet")
    known_sheets = st.session_state.get("available_sheets", ())
    if known_sheets:
        _choose_workbook_sheet()
        selected_sheet = st.session_state.get("selected_sheet")

    if "header_row_widget" not in st.session_state:
        st.session_state["header_row_widget"] = int(
            st.session_state.get("header_row", 1)
        )
    header_row = st.number_input(
        "Header row (1-based)",
        min_value=1,
        step=1,
        key="header_row_widget",
    )
    previous_header_row = st.session_state.get("header_row")
    if previous_header_row is not None and header_row != previous_header_row:
        _reset_upload_state()
    st.session_state["header_row"] = int(header_row)

    try:
        upload = read_upload(
            content,
            file_name,
            sheet_name=selected_sheet,
            header_row=int(header_row),
        )
    except SheetChoiceRequired as exc:
        st.session_state["available_sheets"] = exc.sheet_names
        st.session_state["selected_sheet"] = None
        st.session_state.pop("sheet_picker", None)
        _choose_workbook_sheet()
        return
    except ValueError as exc:
        st.error(str(exc))
        return

    if upload.sheet_name is None:
        st.session_state["selected_sheet"] = None
        st.session_state["available_sheets"] = ()
    elif not st.session_state.get("available_sheets"):
        st.session_state["selected_sheet"] = upload.sheet_name

    st.subheader("Preview")
    preview = preview_rows(upload, 10)
    preview = [
        {key: value for key, value in row.items() if key is not None}
        for row in preview
    ]
    st.caption(f"Header row: {header_row}; first data row: {upload.first_data_row}.")
    st.dataframe(preview, width="stretch")

    header_context = (file_name, content_hash, upload.sheet_name, int(header_row))
    if st.session_state.get("header_context") != header_context:
        _reset_upload_state()
        st.session_state["header_context"] = header_context

    mapping = _mapping_selectboxes(
        upload.fieldnames,
        st.session_state.get("mapping_selection"),
    )
    sex_mapping, load_disabled_reason = _sex_mapping_table(
        upload,
        mapping.get("sex"),
    )

    options = gene_column_options(upload.fieldnames, mapping)
    stored_genes = st.session_state.get(
        "gene_selection_widget",
        st.session_state.get("gene_selection", []),
    )
    st.session_state["gene_selection_widget"] = drop_stale_genes(
        stored_genes,
        options,
    )
    gene_selection = st.multiselect(
        "Gene columns to inspect (not parsed in this step)",
        options=options,
        key="gene_selection_widget",
    )
    st.session_state["gene_selection"] = gene_selection
    for gene, values in summarize_distinct_values(
        upload.rows,
        gene_selection,
    ).items():
        st.write(f"{gene}: {', '.join(values)}")

    signature = input_signature(
        file_name=file_name,
        content_hash=content_hash,
        sheet_name=upload.sheet_name,
        header_row=int(header_row),
        column_mapping=mapping,
        sex_mapping=sex_mapping,
        gene_selection=gene_selection,
    )
    st.caption(_STRUCTURE_NOTICE)
    if load_disabled_reason:
        st.caption(load_disabled_reason)
    if st.button("Load colony", disabled=bool(load_disabled_reason)):
        result = load_from_upload(
            upload,
            mapping,
            gene_selection,
            sex_mapping,
            source_name=file_name,
        )
        st.session_state["loaded_result"] = result
        st.session_state["loaded_signature"] = signature

    stored_result = st.session_state.get("loaded_result")
    if stored_result is not None:
        if st.session_state.get("loaded_signature") == signature:
            _display_load_result(stored_result)
        else:
            st.info("Inputs changed; press Load colony again to refresh the result.")


if __name__ == "__main__":
    main()
