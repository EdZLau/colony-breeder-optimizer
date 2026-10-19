import io

import pytest

from colony_breeder.app.ui_helpers import (
    NONE_COLUMN,
    RESERVED_COLUMNS,
    build_loader_args,
    guess_column_mapping,
    guess_sex_value_mapping,
    load_from_upload,
    summarize_distinct_values,
)
from colony_breeder.app.upload import UploadResult, UploadWarning, read_upload


def _mapping(**overrides):
    mapping = {
        "animal_id": "Animal ID",
        "sex": "Gender",
        "cage": "Cage",
        "mother_id": "Dam",
        "father_id": "Sire",
    }
    mapping.update(overrides)
    return mapping


def test_guess_column_mapping_ignores_case_spaces_and_punctuation():
    assert guess_column_mapping(
        ["Mouse-ID", "Gender", "Cage #", "Mother_ID", "SIRE"]
    ) == {
        "animal_id": "Mouse-ID",
        "sex": "Gender",
        "cage": "Cage #",
        "mother_id": "Mother_ID",
        "father_id": "SIRE",
    }


@pytest.mark.parametrize(
    ("header", "mapping_key"),
    [
        ("Dam ID", "mother_id"),
        ("Sire-ID", "father_id"),
        ("Cage No.", "cage"),
        ("Cage Number", "cage"),
    ],
)
def test_guess_column_mapping_supports_additional_aliases(header, mapping_key):
    assert guess_column_mapping([header])[mapping_key] == header


def test_guess_column_mapping_uses_optional_none_and_required_missing_values():
    assert guess_column_mapping(["Animal", "Sex"]) == {
        "animal_id": "Animal",
        "sex": "Sex",
        "cage": NONE_COLUMN,
        "mother_id": NONE_COLUMN,
        "father_id": NONE_COLUMN,
    }
    assert guess_column_mapping(["Genotype"])["animal_id"] is None


def test_guess_sex_value_mapping_handles_known_unknown_and_blank_values():
    rows = [
        {"Gender": "m"},
        {"Gender": "Male"},
        {"Gender": "f"},
        {"Gender": "female"},
        {"Gender": "unknown"},
        {"Gender": ""},
    ]

    assert guess_sex_value_mapping(rows, "Gender") == {
        "m": "M",
        "Male": "M",
        "f": "F",
        "female": "F",
        "unknown": "",
        "": "",
    }


def test_summarize_distinct_values_caps_at_fifty_and_reports_remainder():
    rows = [{"Cre": f"value-{index:02d}"} for index in range(52)]

    summary = summarize_distinct_values(rows, ["Cre"])

    assert summary["Cre"] == [f"value-{index:02d}" for index in range(50)] + [
        "+2 more"
    ]


def test_build_loader_args_injects_reserved_blank_columns_and_keeps_offset():
    fieldnames = ["Animal ID", "Gender", "Cre"]
    rows = [
        {"Animal ID": "M01", "Gender": "m", "Cre": "not a parsed genotype"}
    ]

    args = build_loader_args(
        fieldnames,
        rows,
        _mapping(cage=NONE_COLUMN, mother_id=NONE_COLUMN, father_id=NONE_COLUMN),
        first_data_row=5,
    )

    assert args["fieldnames"] == fieldnames + [
        RESERVED_COLUMNS["cage"],
        RESERVED_COLUMNS["mother_id"],
        RESERVED_COLUMNS["father_id"],
    ]
    assert args["rows"][0][RESERVED_COLUMNS["cage"]] == ""
    assert args["column_mapping"]["mother_id"] == RESERVED_COLUMNS["mother_id"]
    assert args["gene_columns"] == []
    assert args["first_data_row"] == 5
    assert args["allow_unknown_sex"] is True


def test_build_loader_args_rejects_duplicate_metadata_selections():
    with pytest.raises(ValueError, match="different uploaded column"):
        build_loader_args(
            ["Animal ID", "Gender"],
            [{"Animal ID": "M01", "Gender": "M"}],
            _mapping(sex="Animal ID", cage=NONE_COLUMN, mother_id=NONE_COLUMN, father_id=NONE_COLUMN),
            first_data_row=2,
        )


def test_build_loader_args_rejects_stale_mapping_selection():
    with pytest.raises(ValueError, match="Selected column 'Old ID'.*no longer available"):
        build_loader_args(
            ["Animal ID", "Gender"],
            [{"Animal ID": "M01", "Gender": "M"}],
            _mapping(animal_id="Old ID", cage=NONE_COLUMN, mother_id=NONE_COLUMN, father_id=NONE_COLUMN),
            first_data_row=2,
        )


def test_build_loader_args_rejects_reserved_name_collision():
    reserved = RESERVED_COLUMNS["cage"]
    fieldnames = ["Animal ID", "Gender", reserved]
    rows = [{"Animal ID": "M01", "Gender": "M", reserved: "source value"}]

    with pytest.raises(ValueError, match="Reserved column name"):
        build_loader_args(
            fieldnames,
            rows,
            _mapping(cage=NONE_COLUMN, mother_id=NONE_COLUMN, father_id=NONE_COLUMN),
            first_data_row=2,
        )


def test_load_from_fake_csv_groups_warnings_and_counts_unknown_sex():
    upload = read_upload(
        b"Animal ID,Gender,Cage,Dam,Sire,Cre\n"
        b"M01,male,C1,,F99,not-a-genotype\n"
        b"F01,unreported,C2,,,also-not-a-genotype\n",
        "fake-colony.csv",
    )
    upload.warnings = [
        UploadWarning("Column 'Cre' contains 1 date/time-typed cell(s).", "Cre")
    ]

    result = load_from_upload(
        upload,
        _mapping(),
        gene_columns=["Cre"],
        sex_value_mapping={"male": "M", "unreported": "unknown"},
        source_name="fake-colony.csv",
    )

    assert result.error_message is None
    assert [animal.animal_id for animal in result.animals] == ["M01", "F01"]
    assert result.animals[1].sex == ""
    assert result.unknown_sex_count == 1
    assert result.unknown_sex_ids == ["F01"]
    assert result.warnings["Cre"] == upload.warnings
    assert "1 informational pedigree message(s)." in result.warnings["Pedigree"]
    assert result.distinct_values["Cre"] == ["not-a-genotype", "also-not-a-genotype"]


def test_load_from_upload_only_mentions_parent_columns_for_relatedness():
    upload = UploadResult(
        fieldnames=["Animal ID", "Gender", "Dam", "Sire"],
        rows=[
            {
                "Animal ID": "M01",
                "Gender": "m",
                "Dam": "",
                "Sire": "",
            }
        ],
        warnings=[],
        first_data_row=2,
    )

    result = load_from_upload(
        upload,
        _mapping(cage=NONE_COLUMN),
        gene_columns=[],
        sex_value_mapping={"m": "M"},
    )

    assert result.error_message is None
    assert "Relatedness" not in result.warnings
    assert result.warnings["Cage"] == [
        "The same_cage flag is unavailable because the cage column is unmapped."
    ]


def test_load_from_upload_does_not_modify_upload_rows():
    upload = UploadResult(
        fieldnames=["Animal ID", "Gender", "Cage", "Dam", "Sire", "Cre"],
        rows=[
            {
                "Animal ID": "M01",
                "Gender": "m",
                "Cage": "C1",
                "Dam": "",
                "Sire": "",
                "Cre": "arbitrary text",
            }
        ],
        warnings=[],
        first_data_row=2,
    )
    original_rows = [dict(row) for row in upload.rows]

    result = load_from_upload(
        upload,
        _mapping(),
        gene_columns=["Cre"],
        sex_value_mapping={"m": "M"},
    )

    assert result.error_message is None
    assert upload.rows == original_rows


def test_load_from_upload_counts_info_and_caps_pedigree_warning_details():
    rows = []
    for index in range(55):
        animal_id = f"M{index:02d}"
        rows.append(
            {
                "Animal ID": animal_id,
                "Gender": "m",
                "Cage": "C1",
                "Dam": animal_id,
                "Sire": "missing",
            }
        )
    upload = UploadResult(
        fieldnames=["Animal ID", "Gender", "Cage", "Dam", "Sire"],
        rows=rows,
        warnings=[],
        first_data_row=2,
    )

    result = load_from_upload(
        upload,
        _mapping(),
        gene_columns=[],
        sex_value_mapping={"m": "M"},
    )

    pedigree_messages = result.warnings["Pedigree"]
    assert len([message for message in pedigree_messages if message.startswith("warning:")]) == 50
    assert "+5 more" in pedigree_messages
    assert "55 informational pedigree message(s)." in pedigree_messages


def test_upload_warnings_keep_column_as_structured_metadata():
    warning = UploadWarning("Column 'Cre' contains date cells.", column="Cre")

    result = load_from_upload(
        UploadResult(
            fieldnames=["Animal ID", "Gender", "Cage", "Dam", "Sire", "Cre"],
            rows=[
                {
                    "Animal ID": "M01",
                    "Gender": "m",
                    "Cage": "",
                    "Dam": "",
                    "Sire": "",
                    "Cre": "value",
                }
            ],
            warnings=[warning],
            first_data_row=2,
        ),
        _mapping(),
        gene_columns=[],
        sex_value_mapping={"m": "M"},
    )

    assert warning.column == "Cre"
    assert result.warnings["Cre"] == ["Column 'Cre' contains date cells."]


def test_load_from_fake_xlsx_uses_empty_parent_columns():
    openpyxl = pytest.importorskip("openpyxl")
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.append(["Animal ID", "Gender", "Cre"])
    sheet.append(["M01", "m", "arbitrary symbol"])
    stream = io.BytesIO()
    workbook.save(stream)
    upload = read_upload(stream.getvalue(), "fake-colony.xlsx")

    result = load_from_upload(
        upload,
        _mapping(cage=NONE_COLUMN, mother_id=NONE_COLUMN, father_id=NONE_COLUMN),
        gene_columns=["Cre"],
        sex_value_mapping={"m": "M"},
    )

    assert result.error_message is None
    assert len(result.animals) == 1
    assert result.animals[0].cage == ""
    assert result.animals[0].mother_id == ""
    assert result.animals[0].father_id == ""
    assert "Relatedness" in result.warnings
    assert result.distinct_values["Cre"] == ["arbitrary symbol"]


def test_load_from_upload_returns_stale_selection_as_error_message():
    upload = UploadResult(
        fieldnames=["Animal ID", "Gender"],
        rows=[{"Animal ID": "M01", "Gender": "m"}],
        warnings=[],
        first_data_row=2,
    )

    result = load_from_upload(
        upload,
        _mapping(animal_id="Old ID"),
        gene_columns=[],
        sex_value_mapping={"m": "M"},
    )

    assert result.animals == []
    assert result.error_message is not None
    assert "Old ID" in result.error_message


@pytest.mark.parametrize(
    ("gene_columns", "message"),
    [
        (["Cre", "Cre"], "selected more than once"),
        (["Old gene"], "no longer available"),
        (["Gender"], "already mapped to a colony field"),
    ],
)
def test_load_from_upload_rejects_invalid_gene_selection(gene_columns, message):
    upload = UploadResult(
        fieldnames=["Animal ID", "Gender", "Cre", "Cage", "Dam", "Sire"],
        rows=[
            {
                "Animal ID": "M01",
                "Gender": "m",
                "Cre": "Cre/+",
                "Cage": "C1",
                "Dam": "",
                "Sire": "",
            }
        ],
        warnings=[],
        first_data_row=2,
    )

    result = load_from_upload(
        upload,
        _mapping(),
        gene_columns=gene_columns,
        sex_value_mapping={"m": "M"},
    )

    assert result.animals == []
    assert result.error_message is not None
    assert message in result.error_message
