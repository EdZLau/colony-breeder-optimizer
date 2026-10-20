from dataclasses import FrozenInstanceError
from fractions import Fraction

import pytest

from colony_breeder.app.streamlit_logic import (
    LoadResultView,
    colony_sex_warnings,
    build_ranking_upload_args,
    build_settings_from_widgets,
    default_column_selections,
    default_sex_value_mapping,
    drop_stale_genes,
    drop_stale_run_settings,
    excluded_animal_rows,
    format_flags,
    format_fraction_with_decimal,
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
from colony_breeder.app.ui_helpers import (
    NONE_COLUMN,
    UploadLoadResult,
    guess_column_mapping,
    guess_sex_value_mapping,
    load_from_upload,
)
from colony_breeder.app.run_settings import build_run_settings
from colony_breeder.app.run_ranking import (
    ExcludedAnimalView,
    RankedPairView,
    RunResult,
    SkippedPairView,
)
from colony_breeder.app.upload import UploadResult
from colony_breeder.core.colony import Animal


def _upload(rows, fieldnames=None, warnings=None, first_data_row=2):
    return UploadResult(
        fieldnames=fieldnames or ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
        rows=rows,
        warnings=warnings or [],
        first_data_row=first_data_row,
    )


def _animal(animal_id, sex):
    return Animal(
        animal_id=animal_id,
        sex=sex,
        cage="C1",
        mother_id="",
        father_id="",
    )


def test_default_column_selections_uses_guesses_and_optional_none():
    fieldnames = ["Mouse-ID", "Gender", "Cre"]

    selections = default_column_selections(fieldnames)

    assert selections == {
        "animal_id": "Mouse-ID",
        "sex": "Gender",
        "cage": NONE_COLUMN,
        "mother_id": NONE_COLUMN,
        "father_id": NONE_COLUMN,
    }


def test_default_column_selections_preserves_only_current_previous_choices():
    fieldnames = ["Animal ID", "Gender", "Cage No.", "Dam ID", "Sire-ID"]
    previous = {
        "animal_id": "Animal ID",
        "sex": "old sex",
        "cage": "old cage",
        "mother_id": "Dam ID",
        "father_id": "old father",
    }

    selections = default_column_selections(fieldnames, previous)

    assert selections == {
        "animal_id": "Animal ID",
        "sex": "Gender",
        "cage": "Cage No.",
        "mother_id": "Dam ID",
        "father_id": "Sire-ID",
    }


def test_default_column_selections_preserves_none_for_optional_columns():
    fieldnames = ["Animal ID", "Gender", "Cage", "Dam", "Cre"]
    previous = {
        "animal_id": "Animal ID",
        "sex": "Gender",
        "cage": NONE_COLUMN,
        "mother_id": NONE_COLUMN,
        "father_id": NONE_COLUMN,
    }

    selections = default_column_selections(fieldnames, previous)

    assert selections["cage"] == NONE_COLUMN
    assert selections["mother_id"] == NONE_COLUMN


@pytest.mark.parametrize("key", ["animal_id", "sex"])
def test_default_column_selections_does_not_keep_none_for_required_columns(key):
    fieldnames = ["Animal ID", "Gender", "Cage"]
    guesses = {"animal_id": "Animal ID", "sex": "Gender"}
    previous = {
        "animal_id": "Animal ID",
        "sex": "Gender",
        "cage": "Cage",
        "mother_id": NONE_COLUMN,
        "father_id": NONE_COLUMN,
    }
    previous[key] = NONE_COLUMN

    selections = default_column_selections(fieldnames, previous)

    assert selections[key] == guesses[key]


def test_default_sex_value_mapping_keeps_current_choices_and_guesses_new_values():
    rows = [{"sex": "female"}, {"sex": "male"}, {"sex": ""}]
    previous = {"female": "unknown", "stale": "F"}

    mapping = default_sex_value_mapping(rows, "sex", previous)

    assert mapping == {"female": "unknown", "male": "M", "": ""}
    assert guess_sex_value_mapping(rows, "sex")["female"] == "F"


def test_blank_sex_maps_by_raw_empty_key_and_loads_as_unknown():
    rows = [
        {
            "animal_id": "M01",
            "sex": "",
            "cage": "C1",
            "mother_id": "",
            "father_id": "",
            "Cre": "junk genotype",
        }
    ]
    upload = _upload(rows)
    mapping = guess_column_mapping(upload.fieldnames)
    sex_mapping = default_sex_value_mapping(upload.rows, "sex")

    assert "" in sex_mapping
    assert "(blank)" not in sex_mapping

    result = load_from_upload(
        upload,
        mapping,
        gene_columns=["Cre"],
        sex_value_mapping=sex_mapping,
    )

    assert result.error_message is None
    assert result.unknown_sex_count == 1
    assert result.unknown_sex_ids == ["M01"]
    assert result.animals[0].sex == ""


def test_sex_value_label_is_display_only_for_blank_values():
    assert sex_value_label("") == "(blank)"
    assert sex_value_label("female") == "female"


def test_sex_mapping_check_reports_count_and_displayability():
    rows = [{"sex": f"value-{index}"} for index in range(20)]
    over_limit = [{"sex": f"value-{index}"} for index in range(21)]

    assert sex_mapping_check(rows, "sex") == sex_mapping_check(rows, "sex", limit=20)
    assert sex_mapping_check(rows, "sex").distinct_value_count == 20
    assert sex_mapping_check(rows, "sex").displayable is True
    assert sex_mapping_check(over_limit, "sex").distinct_value_count == 21
    assert sex_mapping_check(over_limit, "sex").displayable is False
    assert "probably isn't a sex column" in sex_mapping_check(
        over_limit, "sex"
    ).warning


def test_gene_column_options_excludes_every_mapped_colony_column():
    fieldnames = ["ID", "Sex", "Cage", "Dam", "Sire", "Cre", "Ob"]
    mapping = {
        "animal_id": "ID",
        "sex": "Sex",
        "cage": "Cage",
        "mother_id": "Dam",
        "father_id": "Sire",
    }

    assert gene_column_options(fieldnames, mapping) == ["Cre", "Ob"]


def test_gene_options_exclude_phase_columns_and_phase_options_exclude_gene_pool():
    fieldnames = ["ID", "Sex", "Cage", "A", "B", "AB_phase", "Cre"]
    mapping = {
        "animal_id": "ID",
        "sex": "Sex",
        "cage": "Cage",
        "mother_id": NONE_COLUMN,
        "father_id": NONE_COLUMN,
    }

    assert gene_column_options(
        fieldnames,
        mapping,
        excluded_columns=["AB_phase"],
    ) == ["A", "B", "Cre"]
    assert phase_column_options(fieldnames, mapping, ["A", "B"]) == [
        "AB_phase",
        "Cre",
    ]


def test_drop_stale_run_settings_when_pool_headers_or_mapping_change():
    settings = build_run_settings(
        targets={"A": ["A/a"], "Cre": ["Cre/+"]},
        lethals={"B": ["B/B"]},
        x_linked_genes=["X", "A"],
        linked_pairs=[
            ("A", "B", "1/10", "AB_phase"),
            ("A", "Cre", "1/10", "Cre_phase"),
        ],
        exclude_related=True,
        top=5,
    )
    fieldnames = ["ID", "Sex", "Cage", "A", "B", "AB_phase", "Cre", "Cre_phase"]
    mapping = {
        "animal_id": "ID",
        "sex": "Sex",
        "cage": "Cage",
        "mother_id": NONE_COLUMN,
        "father_id": NONE_COLUMN,
    }

    current = drop_stale_run_settings(
        settings,
        gene_pool=["A", "B"],
        fieldnames=fieldnames,
        column_mapping=mapping,
    )

    assert dict(current.targets) == {"A": ("A/a",)}
    assert dict(current.lethals) == {"B": ("B/B",)}
    assert current.x_linked_genes == ("A",)
    assert current.linked_pairs == (("A", "B", "1/10", "AB_phase"),)
    assert current.exclude_related is True
    assert current.top == 5


def test_build_settings_from_widget_text_and_linked_rows():
    settings = build_settings_from_widgets(
        target_genotypes={"Cre": "Cre/+\n\n+/+"},
        lethal_genotypes={"Ob": " ob/ob \n"},
        x_linked_genes=["Xgene"],
        linked_pairs=[("Cre", "Ob", " 0.1 ", "Cre_Ob_phase")],
        exclude_related=True,
        top=10,
    )

    assert dict(settings.targets) == {"Cre": ("Cre/+", "+/+")}
    assert dict(settings.lethals) == {"Ob": ("ob/ob",)}
    assert settings.x_linked_genes == ("Xgene",)
    assert settings.linked_pairs == (
        ("Cre", "Ob", "0.1", "Cre_Ob_phase"),
    )
    assert settings.exclude_related is True
    assert settings.top == 10


def test_drop_stale_genes_preserves_order_of_current_options():
    assert drop_stale_genes(["Ob", "old", "Cre"], ["Cre", "Ob"]) == [
        "Ob",
        "Cre",
    ]


def test_input_signature_is_stable_and_changes_for_each_input():
    base = {
        "file_name": "colony.csv",
        "content_hash": "a" * 64,
        "sheet_name": None,
        "header_row": 1,
        "column_mapping": {"animal_id": "ID", "sex": "Sex"},
        "sex_mapping": {"F": "F"},
        "gene_selection": ["Cre"],
    }
    original = input_signature(**base)
    assert input_signature(**base) == original
    changes = [
        {"file_name": "other.csv"},
        {"content_hash": "b" * 64},
        {"sheet_name": "Sheet1"},
        {"header_row": 2},
        {"column_mapping": {"animal_id": "Other ID", "sex": "Sex"}},
        {"sex_mapping": {"F": "unknown"}},
        {"gene_selection": ["Ob"]},
    ]

    for change in changes:
        updated = dict(base)
        updated.update(change)
        assert input_signature(**updated) != original


def test_input_signature_distinguishes_same_name_same_size_content_changes():
    common = {
        "file_name": "colony.csv",
        "sheet_name": None,
        "header_row": 1,
        "column_mapping": {"animal_id": "ID", "sex": "Sex"},
        "sex_mapping": {"F": "F"},
        "gene_selection": ["Cre"],
    }

    first = input_signature(content_hash="a" * 64, **common)
    second = input_signature(content_hash="b" * 64, **common)

    assert first != second


def test_preview_rows_returns_bounded_deep_copies_without_mutating_upload():
    upload = _upload(
        [
            {"animal_id": "M01", "sex": "M", None: ["overflow"]},
            {"animal_id": "F01", "sex": "F"},
        ]
    )

    preview = preview_rows(upload, limit=1)
    preview[0]["animal_id"] = "CHANGED"
    preview[0][None].append("also changed")

    assert len(preview) == 1
    assert upload.rows[0]["animal_id"] == "M01"
    assert upload.rows[0][None] == ["overflow"]


def test_format_load_result_returns_frozen_view_and_preserves_warning_text():
    warning_text = "Relatedness flags unavailable because mother_id is unmapped."
    result = UploadLoadResult(
        animals=[_animal("M01", "")],
        warnings={"Relatedness": [warning_text], "Cre": ["Column warning text."]},
        unknown_sex_count=1,
        unknown_sex_ids=["M01"],
        distinct_values={"Cre": ["junk genotype"]},
        error_message=None,
    )

    view = format_load_result(result)

    assert isinstance(view, LoadResultView)
    assert view.animal_count == 1
    assert view.unknown_sex_count == 1
    assert view.unknown_sex_ids == ("M01",)
    assert view.grouped_warnings == (
        ("Relatedness", (warning_text,)),
        ("Cre", ("Column warning text.",)),
    )
    assert view.gene_value_summaries == (("Cre", ("junk genotype",)),)
    with pytest.raises(FrozenInstanceError):
        view.animal_count = 2


def test_colony_sex_warnings_for_missing_sex_groups():
    assert any("no males" in warning.lower() for warning in colony_sex_warnings([_animal("F01", "F")]))
    assert any("no females" in warning.lower() for warning in colony_sex_warnings([_animal("M01", "M")]))
    assert colony_sex_warnings([_animal("M01", "M"), _animal("F01", "F")]) == []


def test_build_ranking_upload_args_uses_loaded_sex_mapping_without_mutating_upload():
    rows = [
        {
            "ID": "M01",
            "Gender": "male",
            "Cage": "C1",
            "Dam": "",
            "Sire": "",
            "Cre": "Cre/+",
        }
    ]
    fieldnames = ["ID", "Gender", "Cage", "Dam", "Sire", "Cre"]
    upload = _upload(rows, fieldnames=fieldnames)
    mapping = {
        "animal_id": "ID",
        "sex": "Gender",
        "cage": "Cage",
        "mother_id": "Dam",
        "father_id": "Sire",
    }

    args = build_ranking_upload_args(
        upload,
        mapping,
        [_animal("M01", "M")],
    )

    assert args["rows"][0]["Gender"] == "M"
    assert upload.rows[0]["Gender"] == "male"


def test_ranking_formatters_preserve_probabilities_messages_flags_and_all_rows():
    settings = build_run_settings(
        targets={"A": "A/A"},
        linked_pairs=[("A", "B", "1/10", "AB_phase")],
        top=1,
    )
    result = RunResult(
        ranked=(
            RankedPairView(
                "M01",
                "F01",
                Fraction(1, 2),
                (("A", Fraction(1, 2)), ("B", Fraction(1, 4))),
                ((("A", "B"), Fraction(1, 2)),),
                ("full_siblings",),
            ),
            RankedPairView(
                "M02",
                "F02",
                Fraction(1, 4),
                (),
                (),
                (),
            ),
        ),
        skipped=(
            SkippedPairView("M01", "F03", "missing_phase", "Phase is blank.", ()),
            SkippedPairView("M02", "F04", "all_lethal", "All offspring are lethal.", ("same_cage",)),
        ),
        excluded_animals=(ExcludedAnimalView("U01", "sex is blank"),),
        notes=("Core note, unchanged.",),
        pedigree_warnings=(("warning", "Core pedigree warning."),),
        signature=("current",),
        error_message=None,
    )

    assert format_fraction_with_decimal(Fraction(1, 2)) == "1/2 (0.5000)"
    assert split_error_message("First issue; Second issue") == (
        "First issue",
        "Second issue",
    )
    assert format_flags(("full_siblings", "same_cage")) == "full_siblings, same_cage"
    assert format_flags(()) == "-"

    ranked_rows = ranked_pair_rows(result, settings)
    assert len(ranked_rows) == 1
    assert ranked_rows[0]["Male"] == "M01"
    assert ranked_rows[0]["Probability (exact; decimal)"] == "1/2 (0.5000)"
    assert ranked_rows[0]["Flags"] == "full_siblings"
    assert "A=1/2 (descriptive)" in ranked_rows[0][
        "Per-gene values (linked/X values marked descriptive; do not multiply those)"
    ]
    assert ranked_rows[0]["Linked-pair factors"] == "A,B=1/2 (0.5000)"

    skipped_rows = skipped_pair_rows(result.skipped)
    assert len(skipped_rows) == len(result.skipped)
    assert skipped_rows[1]["Flags"] == "same_cage"
    assert excluded_animal_rows(result.excluded_animals) == [
        {"Animal": "U01", "Reason": "sex is blank"}
    ]
    summary = run_result_summary(result, top=1)
    assert summary[:6] == (
        "Ranked pairs: 2",
        "Visible nonzero pairs: 1",
        "Zero-probability pairs hidden: 0",
        "Pairs hidden by top: 1",
        "Skipped pairs: 2",
        "Excluded animals: 1",
    )
    assert "missing_genotype=0" in summary[6]
    assert "missing_phase=1" in summary[6]
    assert "all_lethal=1" in summary[6]


def test_stale_linked_pair_notice_names_the_removed_pair():
    previous = build_run_settings(
        targets={"A": "A/a"},
        linked_pairs=[
            ("A", "B", "1/10", "AB_phase"),
            ("C", "D", "1/10", "CD_phase"),
        ],
    )
    current = build_run_settings(
        targets={"A": "A/a"},
        linked_pairs=[("A", "B", "1/10", "AB_phase")],
    )

    assert stale_linked_pair_notices(previous, current) == (
        "Linked pair 2 was removed because its gene or phase column is no longer available.",
    )
