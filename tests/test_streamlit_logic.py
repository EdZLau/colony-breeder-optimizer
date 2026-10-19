from dataclasses import FrozenInstanceError

import pytest

from colony_breeder.app.streamlit_logic import (
    LoadResultView,
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
    guess_column_mapping,
    guess_sex_value_mapping,
    load_from_upload,
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
