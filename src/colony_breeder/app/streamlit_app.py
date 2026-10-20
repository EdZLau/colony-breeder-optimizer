"""Streamlit page for loading and ranking an uploaded colony."""

import hashlib
from typing import Dict, Optional, Sequence, Tuple

import streamlit as st

from colony_breeder.app.streamlit_logic import (
    build_settings_from_widgets,
    colony_sex_warnings,
    default_column_selections,
    default_sex_value_mapping,
    drop_stale_genes,
    drop_stale_run_settings,
    build_ranking_upload_args,
    excluded_animal_rows,
    format_load_result,
    gene_column_options,
    input_signature,
    phase_column_options,
    preview_rows,
    ranked_pair_rows,
    run_result_summary,
    sex_mapping_check,
    sex_value_label,
    skipped_pair_rows,
    split_error_message,
    stale_linked_pair_notices,
)
from colony_breeder.app.run_ranking import (
    RunResult,
    result_is_current,
    run_ranking,
)
from colony_breeder.app.run_settings import (
    RunSettings,
    run_signature,
    validate_run_settings,
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
_SETTINGS_WIDGET_KEYS = (
    "target_genes_widget",
    "lethal_genes_widget",
    "x_linked_genes_widget",
    "linked_pair_count_widget",
    "exclude_related_widget",
    "top_widget",
)
_LINKED_WIDGET_PREFIXES = (
    "linked_gene1_",
    "linked_gene2_",
    "linked_r_",
    "linked_phase_",
)
_SETTINGS_WIDGET_PREFIXES = (
    "target_genotype_",
    "lethal_genotype_",
    *_LINKED_WIDGET_PREFIXES,
)
_STRUCTURE_NOTICE = (
    "Colony structure is loaded. Genotypes are parsed and validated when you "
    "run ranking."
)


def _clear_settings_widgets() -> None:
    for key in _SETTINGS_WIDGET_KEYS:
        st.session_state.pop(key, None)
    for key in list(st.session_state):
        if key.startswith(_SETTINGS_WIDGET_PREFIXES):
            st.session_state.pop(key, None)


def _reset_upload_state() -> None:
    for key in (
        "mapping_selection",
        "sex_value_mapping",
        "gene_selection",
        "loaded_result",
        "loaded_signature",
        "header_context",
        "run_settings_draft",
        "settings_context",
        "settings_signature",
        "run_result",
        "run_result_signature",
        "run_error",
    ):
        st.session_state.pop(key, None)
    for key in _MAPPING_WIDGETS:
        st.session_state.pop(f"mapping_widget_{key}", None)
    _clear_settings_widgets()
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


def _display_run_error(message: str) -> None:
    errors = split_error_message(message)
    if len(errors) > 1:
        st.error("Ranking could not be completed:")
        st.markdown("\n".join(f"- {error}" for error in errors))
    else:
        st.error(message)


def _display_run_result(result: RunResult, settings: RunSettings) -> None:
    st.subheader("Ranking results")
    for line in run_result_summary(result, settings.top):
        st.write(line)

    for note in result.notes:
        st.info(note)
    if result.pedigree_warnings:
        st.markdown("**Pedigree warnings**")
        for severity, message in result.pedigree_warnings:
            st.markdown(f"**{severity}**")
            st.write(message)

    rows = ranked_pair_rows(result, settings)
    if rows:
        st.dataframe(rows, width="stretch")
    else:
        st.write("No nonzero-probability pairs to display.")

    with st.expander(f"Skipped pairs ({len(result.skipped)})"):
        rows = skipped_pair_rows(result.skipped)
        if rows:
            st.dataframe(rows, width="stretch")
        else:
            st.write("No skipped pairs.")

    with st.expander(f"Excluded animals ({len(result.excluded_animals)})"):
        rows = excluded_animal_rows(result.excluded_animals)
        if rows:
            st.dataframe(rows, width="stretch")
        else:
            st.write("No excluded animals.")


def _settings_controls(
    upload: UploadResult,
    column_mapping: Dict[str, Optional[str]],
    gene_pool: Sequence[str],
    settings_context: tuple,
    base_signature: tuple,
    loaded_result: UploadLoadResult,
) -> None:
    st.header("Breeding settings")
    st.caption(
        "Untargeted autosomal lethals are ignored unless the gene is in a linked "
        "pair. Untargeted lethals on X-linked or linked genes affect survival."
    )
    st.caption(
        "In your spreadsheet, X-linked males must be written A/Y and females "
        "A/a or A/A; target entries may use either form. "
        "X-linked genes are treated as unlinked to each other."
    )
    st.caption(
        '"(none)" for mother or father disables parent-based relatedness flags; '
        '"(none)" for cage disables same-cage flags. Linked pairs are assumed '
        "independent of other genes."
    )

    if st.session_state.get("settings_context") != settings_context:
        _clear_settings_widgets()
        st.session_state.pop("run_settings_draft", None)
        st.session_state.pop("settings_signature", None)
        st.session_state["settings_context"] = settings_context

    previous_settings = st.session_state.get("run_settings_draft")
    prior_phases = (
        [pair[3] for pair in previous_settings.linked_pairs if pair[3]]
        if isinstance(previous_settings, RunSettings)
        else []
    )
    available_gene_options = gene_column_options(
        upload.fieldnames,
        column_mapping,
        excluded_columns=prior_phases,
    )
    eligible_gene_pool = drop_stale_genes(gene_pool, available_gene_options)
    if isinstance(previous_settings, RunSettings):
        previous_settings = drop_stale_run_settings(
            previous_settings,
            eligible_gene_pool,
            upload.fieldnames,
            column_mapping,
        )
        for notice in stale_linked_pair_notices(
            st.session_state.get("run_settings_draft"),
            previous_settings,
        ):
            st.info(notice)
        st.session_state["run_settings_draft"] = previous_settings

    phase_options = phase_column_options(
        upload.fieldnames,
        column_mapping,
        eligible_gene_pool,
    )

    prior_targets = (
        dict(previous_settings.targets)
        if isinstance(previous_settings, RunSettings)
        else {}
    )
    prior_lethals = (
        dict(previous_settings.lethals)
        if isinstance(previous_settings, RunSettings)
        else {}
    )
    prior_target_genes = [
        gene for gene in prior_targets if gene in eligible_gene_pool
    ]
    prior_lethal_genes = [
        gene for gene in prior_lethals if gene in eligible_gene_pool
    ]

    target_key = "target_genes_widget"
    if target_key not in st.session_state:
        st.session_state[target_key] = drop_stale_genes(
            prior_target_genes,
            eligible_gene_pool,
        )
    else:
        st.session_state[target_key] = drop_stale_genes(
            st.session_state[target_key],
            eligible_gene_pool,
        )
    target_genes = st.multiselect(
        "Target genes",
        options=eligible_gene_pool,
        key=target_key,
    )
    for key in list(st.session_state):
        if key.startswith("target_genotype_") and key.removeprefix(
            "target_genotype_"
        ) not in target_genes:
            st.session_state.pop(key, None)
    target_texts = {}
    for gene in target_genes:
        widget_key = f"target_genotype_{gene}"
        if widget_key not in st.session_state:
            st.session_state[widget_key] = "\n".join(prior_targets.get(gene, ()))
        examples = summarize_distinct_values(upload.rows, [gene]).get(gene, [])
        target_texts[gene] = st.text_area(
            f"Target genotypes for {gene} (one per line)",
            key=widget_key,
            help=(
                "Enter one genotype per line; multiple lines are alternatives. "
                f"Observed values (examples only): {', '.join(examples) or 'none'}"
            ),
        )

    lethal_key = "lethal_genes_widget"
    if lethal_key not in st.session_state:
        st.session_state[lethal_key] = drop_stale_genes(
            prior_lethal_genes,
            eligible_gene_pool,
        )
    else:
        st.session_state[lethal_key] = drop_stale_genes(
            st.session_state[lethal_key],
            eligible_gene_pool,
        )
    lethal_genes = st.multiselect(
        "Genes with lethal genotypes (optional)",
        options=eligible_gene_pool,
        key=lethal_key,
    )
    for key in list(st.session_state):
        if key.startswith("lethal_genotype_") and key.removeprefix(
            "lethal_genotype_"
        ) not in lethal_genes:
            st.session_state.pop(key, None)
    lethal_texts = {}
    for gene in lethal_genes:
        widget_key = f"lethal_genotype_{gene}"
        if widget_key not in st.session_state:
            st.session_state[widget_key] = "\n".join(prior_lethals.get(gene, ()))
        examples = summarize_distinct_values(upload.rows, [gene]).get(gene, [])
        lethal_texts[gene] = st.text_area(
            f"Lethal genotypes for {gene} (one per line)",
            key=widget_key,
            help=(
                "Optional; enter one genotype per line. "
                f"Observed values (examples only): {', '.join(examples) or 'none'}"
            ),
        )

    prior_x_linked = (
        list(previous_settings.x_linked_genes)
        if isinstance(previous_settings, RunSettings)
        else []
    )
    x_linked_key = "x_linked_genes_widget"
    if x_linked_key not in st.session_state:
        st.session_state[x_linked_key] = drop_stale_genes(
            prior_x_linked,
            eligible_gene_pool,
        )
    else:
        st.session_state[x_linked_key] = drop_stale_genes(
            st.session_state[x_linked_key],
            eligible_gene_pool,
        )
    x_linked_genes = st.multiselect(
        "X-linked genes",
        options=eligible_gene_pool,
        key=x_linked_key,
    )

    prior_pairs = (
        list(previous_settings.linked_pairs)
        if isinstance(previous_settings, RunSettings)
        else []
    )
    pair_count_key = "linked_pair_count_widget"
    if pair_count_key not in st.session_state:
        st.session_state[pair_count_key] = min(len(prior_pairs), 3)
    pair_count = st.number_input(
        "Linked gene pairs (up to 3)",
        min_value=0,
        max_value=3,
        step=1,
        key=pair_count_key,
    )
    for key in list(st.session_state):
        for prefix in _LINKED_WIDGET_PREFIXES:
            if key.startswith(prefix):
                try:
                    row_index = int(key[len(prefix):])
                except ValueError:
                    continue
                if row_index >= int(pair_count):
                    st.session_state.pop(key, None)
    linked_pairs = []
    for index in range(int(pair_count)):
        previous_pair = prior_pairs[index] if index < len(prior_pairs) else (
            "",
            "",
            "",
            "",
        )
        st.markdown(f"**Linked pair {index + 1}**")
        col1, col2 = st.columns(2)
        gene1_key = f"linked_gene1_{index}"
        gene2_key = f"linked_gene2_{index}"
        phase_key = f"linked_phase_{index}"
        r_key = f"linked_r_{index}"

        gene_select_options = [NONE_COLUMN, *eligible_gene_pool]
        phase_select_options = [NONE_COLUMN, *phase_options]
        if r_key not in st.session_state:
            st.session_state[r_key] = previous_pair[2]
        for key, options, previous_value in (
            (gene1_key, gene_select_options, previous_pair[0]),
            (gene2_key, gene_select_options, previous_pair[1]),
            (phase_key, phase_select_options, previous_pair[3]),
        ):
            if st.session_state.get(key) not in options:
                st.session_state[key] = (
                    previous_value if previous_value in options else NONE_COLUMN
                )

        with col1:
            gene1 = st.selectbox(
                "Gene 1",
                options=gene_select_options,
                key=gene1_key,
            )
            gene2 = st.selectbox(
                "Gene 2",
                options=gene_select_options,
                key=gene2_key,
            )
        with col2:
            recombination_text = st.text_input(
                "Recombination fraction r",
                key=r_key,
                help=(
                    "recombination fraction from 0 to 0.5, not map distance; "
                    "one value is used for both sexes."
                ),
            )
            phase_column = st.selectbox(
                "Phase column",
                options=phase_select_options,
                key=phase_key,
            )
        st.caption("Phase text must follow gene1|gene2 order.")
        row_has_setting = (
            gene1 != NONE_COLUMN
            or gene2 != NONE_COLUMN
            or bool(recombination_text.strip())
            or phase_column != NONE_COLUMN
        )
        if row_has_setting:
            linked_pairs.append(
                (
                    "" if gene1 == NONE_COLUMN else gene1,
                    "" if gene2 == NONE_COLUMN else gene2,
                    recombination_text,
                    "" if phase_column == NONE_COLUMN else phase_column,
                )
            )

    if "exclude_related_widget" not in st.session_state:
        st.session_state["exclude_related_widget"] = (
            previous_settings.exclude_related
            if isinstance(previous_settings, RunSettings)
            else False
        )
    exclude_related = st.checkbox(
        "Exclude related pairs",
        key="exclude_related_widget",
    )
    if "top_widget" not in st.session_state:
        st.session_state["top_widget"] = (
            previous_settings.top
            if isinstance(previous_settings, RunSettings)
            and previous_settings.top is not None
            else 10
        )
    top = st.number_input(
        "Top N pairs to show",
        min_value=1,
        step=1,
        key="top_widget",
    )

    settings = build_settings_from_widgets(
        target_texts,
        lethal_texts,
        x_linked_genes,
        linked_pairs,
        exclude_related,
        int(top),
    )
    settings = drop_stale_run_settings(
        settings,
        eligible_gene_pool,
        upload.fieldnames,
        column_mapping,
    )
    st.session_state["run_settings_draft"] = settings
    errors = validate_run_settings(
        settings,
        eligible_gene_pool,
        upload.fieldnames,
        column_mapping,
    )
    if errors:
        st.session_state.pop("settings_signature", None)
        st.error("Settings need attention:")
        st.markdown("\n".join(f"- {error}" for error in errors))
    else:
        st.session_state["settings_signature"] = run_signature(
            base_signature,
            settings,
        )
        st.success("Settings valid.")

    run_clicked = st.button("Run ranking", disabled=bool(errors))
    if run_clicked:
        st.session_state.pop("run_result", None)
        st.session_state.pop("run_result_signature", None)
        st.session_state.pop("run_error", None)
        try:
            upload_args = build_ranking_upload_args(
                upload,
                column_mapping,
                loaded_result.animals,
            )
        except ValueError as exc:
            st.session_state["run_error"] = str(exc)
            st.session_state["run_result_signature"] = st.session_state[
                "settings_signature"
            ]
        else:
            result = run_ranking(
                upload_args,
                settings,
                eligible_gene_pool,
                upload.fieldnames,
                column_mapping,
                base_signature,
            )
            st.session_state["run_result"] = result
            st.session_state["run_result_signature"] = st.session_state[
                "settings_signature"
            ]
            if result.error_message:
                st.session_state["run_error"] = result.error_message

    current_signature = st.session_state.get("settings_signature")
    stored_result_signature = st.session_state.get("run_result_signature")
    if (
        current_signature is not None
        and stored_result_signature == current_signature
        and st.session_state.get("run_error") is not None
    ):
        _display_run_error(st.session_state["run_error"])

    stored_result = st.session_state.get("run_result")
    if (
        stored_result is not None
        and current_signature is not None
        and stored_result_signature == current_signature
        and stored_result.error_message is None
        and result_is_current(stored_result, current_signature)
    ):
        _display_run_result(stored_result, settings)


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

    saved_settings = st.session_state.get("run_settings_draft")
    saved_phase_columns = (
        [pair[3] for pair in saved_settings.linked_pairs if pair[3]]
        if isinstance(saved_settings, RunSettings)
        else []
    )
    options = gene_column_options(
        upload.fieldnames,
        mapping,
        excluded_columns=saved_phase_columns,
    )
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
            if stored_result.error_message is None:
                settings_context = (
                    file_name,
                    content_hash,
                    upload.sheet_name,
                    int(header_row),
                )
                _settings_controls(
                    upload,
                    mapping,
                    gene_selection,
                    settings_context,
                    signature,
                    stored_result,
                )
        else:
            st.info("Inputs changed; press Load colony again to refresh the result.")


if __name__ == "__main__":
    main()
