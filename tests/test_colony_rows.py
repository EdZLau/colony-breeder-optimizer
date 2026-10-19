import csv
import io
from pathlib import Path

import pytest

from colony_breeder.core.colony import load_colony, load_colony_rows

DATA = Path(__file__).resolve().parent.parent / "data"


def test_load_colony_rows_matches_load_colony(tmp_path):
    csv_path = tmp_path / "colony.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        "F01,F,C2,,,+/+\n",
        encoding="utf-8",
    )
    loaded_from_file = load_colony(csv_path, gene_columns=["Cre"])

    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [name.strip() if name is not None else "" for name in reader.fieldnames]
        loaded_from_rows = load_colony_rows(
            fieldnames,
            reader,
            gene_columns=["Cre"],
            source_name=str(csv_path),
        )

    assert loaded_from_rows == loaded_from_file


@pytest.mark.parametrize(
    ("filename", "settings"),
    [
        (
            "sample_mice.csv",
            {
                "gene_columns": ["Cre", "Ob", "Agouti", "Xgene"],
                "x_linked_genes": ["Xgene"],
            },
        ),
        (
            "sample_linked_FAKE.csv",
            {
                "gene_columns": ["A", "B", "C", "D", "E"],
                "linked_pairs": [("A", "B", "1/10", "AB_phase")],
            },
        ),
    ],
)
def test_load_colony_rows_matches_file_loader_for_sample_configs(filename, settings):
    csv_path = DATA / filename
    loaded_from_file = load_colony(csv_path, **settings)

    with csv_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = [
            fieldname.strip() if fieldname is not None else ""
            for fieldname in reader.fieldnames
        ]
        loaded_from_rows = load_colony_rows(
            fieldnames,
            reader,
            source_name=str(csv_path),
            **settings,
        )

    assert loaded_from_rows == loaded_from_file
    if filename == "sample_linked_FAKE.csv":
        assert loaded_from_rows[0].phases[("A", "B")] == "A|B/a|b"


def test_load_colony_rows_supports_custom_column_mapping():
    fieldnames = ["Mouse ID", "Gender", "Cage No", "Dam", "Sire", "Cre"]
    rows = [
        {
            "Mouse ID": "M01",
            "Gender": "M",
            "Cage No": "C1",
            "Dam": "F02",
            "Sire": "M02",
            "Cre": "Cre/+",
        }
    ]
    mapping = {
        "animal_id": "Mouse ID",
        "sex": "Gender",
        "cage": "Cage No",
        "mother_id": "Dam",
        "father_id": "Sire",
    }

    colony = load_colony_rows(
        fieldnames,
        rows,
        column_mapping=mapping,
        gene_columns=["Cre"],
    )

    assert colony[0].animal_id == "M01"
    assert colony[0].sex == "M"
    assert colony[0].cage == "C1"
    assert colony[0].mother_id == "F02"
    assert colony[0].father_id == "M02"
    assert colony[0].genotypes["Cre"] == "Cre/+"


def test_load_colony_rows_reports_missing_mapped_column():
    mapping = {"animal_id": "Mouse ID"}

    with pytest.raises(ValueError, match="Required columns missing.*Mouse ID"):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            [],
            column_mapping=mapping,
            gene_columns=["Cre"],
        )


def test_load_colony_rows_duplicate_id_error_names_source():
    rows = [
        {
            "animal_id": "M01",
            "sex": "M",
            "cage": "C1",
            "mother_id": "",
            "father_id": "",
            "Cre": "Cre/+",
        },
        {
            "animal_id": "M01",
            "sex": "F",
            "cage": "C2",
            "mother_id": "",
            "father_id": "",
            "Cre": "Cre/+",
        },
    ]

    with pytest.raises(ValueError, match="Duplicate animal_id.*uploaded batch.csv"):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            rows,
            gene_columns=["Cre"],
            source_name="uploaded batch.csv",
        )


def test_load_colony_rows_rejects_non_text_cell_with_row_and_column():
    rows = [
        {
            "animal_id": "M01",
            "sex": 1,
            "cage": "C1",
            "mother_id": "",
            "father_id": "",
            "Cre": "Cre/+",
        }
    ]

    with pytest.raises(ValueError, match=r"Row 2.*column 'sex'"):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            rows,
            gene_columns=["Cre"],
        )


def test_load_colony_rows_skips_blank_lines_and_starts_data_at_row_two():
    content = (
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "\n"
        ",M,C1,,,Cre/+\n"
    )
    reader = csv.DictReader(io.StringIO(content))

    with pytest.raises(ValueError, match=r"Row 2 is missing an animal_id"):
        load_colony_rows(
            reader.fieldnames,
            reader,
            gene_columns=["Cre"],
        )


def test_load_colony_rows_uses_custom_first_data_row_for_validation():
    rows = [
        {
            "animal_id": "",
            "sex": "F",
            "cage": "C1",
            "mother_id": "",
            "father_id": "",
            "Cre": "Cre/+",
        }
    ]

    with pytest.raises(ValueError, match=r"Row 9 is missing an animal_id"):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            rows,
            gene_columns=["Cre"],
            first_data_row=9,
        )


def _row_for_genotype_error(animal_id, genotype):
    return {
        "animal_id": animal_id,
        "sex": "M",
        "cage": "C1",
        "mother_id": "",
        "father_id": "",
        "Cre": genotype,
    }


def _row_for_phase_error(animal_id, phase):
    return {
        "animal_id": animal_id,
        "sex": "M",
        "cage": "C1",
        "mother_id": "",
        "father_id": "",
        "A": "A/a",
        "B": "B/b",
        "AB_phase": phase,
    }


def test_malformed_genotype_error_uses_non_first_header_source_row():
    with pytest.raises(
        ValueError,
        match=r"^Row 4: Invalid genotype for animal 'M01' in column 'Cre': 'Cre\+'\.",
    ):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            [_row_for_genotype_error("M01", "Cre+")],
            gene_columns=["Cre"],
            first_data_row=4,
        )


def test_malformed_genotype_error_counts_an_interior_blank_row():
    rows = [
        _row_for_genotype_error("M00", "Cre/+"),
        {
            "animal_id": "",
            "sex": "",
            "cage": "",
            "mother_id": "",
            "father_id": "",
            "Cre": "",
        },
        _row_for_genotype_error("M01", "Cre+"),
    ]

    with pytest.raises(
        ValueError,
        match=r"^Row 4: Invalid genotype for animal 'M01' in column 'Cre': 'Cre\+'\.",
    ):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            rows,
            gene_columns=["Cre"],
            first_data_row=2,
        )


def test_malformed_genotype_error_uses_plain_csv_source_row(tmp_path):
    csv_path = tmp_path / "malformed_genotype.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre+\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match=r"^Row 2: Invalid genotype for animal 'M01' in column 'Cre': 'Cre\+'\.",
    ):
        load_colony(csv_path, gene_columns=["Cre"])


def test_malformed_phase_error_uses_non_first_header_source_row():
    with pytest.raises(
        ValueError,
        match=(
            r"^Row 4: Invalid phase for animal 'M01' in gene pair \('A', 'B'\): "
            r"'bad phase'\."
        ),
    ):
        load_colony_rows(
            [
                "animal_id",
                "sex",
                "cage",
                "mother_id",
                "father_id",
                "A",
                "B",
                "AB_phase",
            ],
            [_row_for_phase_error("M01", "bad phase")],
            gene_columns=["A", "B"],
            linked_pairs=[("A", "B", "1/10", "AB_phase")],
            first_data_row=4,
        )


def test_malformed_phase_error_counts_an_interior_blank_row():
    rows = [
        _row_for_phase_error("M00", "A|B/a|b"),
        {
            "animal_id": "",
            "sex": "",
            "cage": "",
            "mother_id": "",
            "father_id": "",
            "A": "",
            "B": "",
            "AB_phase": "",
        },
        _row_for_phase_error("M01", "bad phase"),
    ]

    with pytest.raises(
        ValueError,
        match=(
            r"^Row 4: Invalid phase for animal 'M01' in gene pair \('A', 'B'\): "
            r"'bad phase'\."
        ),
    ):
        load_colony_rows(
            [
                "animal_id",
                "sex",
                "cage",
                "mother_id",
                "father_id",
                "A",
                "B",
                "AB_phase",
            ],
            rows,
            gene_columns=["A", "B"],
            linked_pairs=[("A", "B", "1/10", "AB_phase")],
            first_data_row=2,
        )


def test_malformed_phase_error_uses_plain_csv_source_row(tmp_path):
    csv_path = tmp_path / "malformed_phase.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,A,B,AB_phase\n"
        "M01,M,C1,,,A/a,B/b,bad phase\n",
        encoding="utf-8",
    )

    with pytest.raises(
        ValueError,
        match=(
            r"^Row 2: Invalid phase for animal 'M01' in gene pair \('A', 'B'\): "
            r"'bad phase'\."
        ),
    ):
        load_colony(
            csv_path,
            gene_columns=["A", "B"],
            linked_pairs=[("A", "B", "1/10", "AB_phase")],
        )


def test_load_colony_rows_allows_unknown_sex_when_requested():
    colony = load_colony_rows(
        ["animal_id", "sex", "cage", "mother_id", "father_id"],
        [
            {
                "animal_id": "M01",
                "sex": "Banana",
                "cage": "C1",
                "mother_id": "",
                "father_id": "",
            }
        ],
        gene_columns=[],
        allow_unknown_sex=True,
    )

    assert colony[0].sex == ""


def test_load_colony_rows_rejects_unknown_sex_by_default():
    with pytest.raises(ValueError, match="Invalid sex value 'Banana'"):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id"],
            [
                {
                    "animal_id": "M01",
                    "sex": "Banana",
                    "cage": "C1",
                    "mother_id": "",
                    "father_id": "",
                }
            ],
            gene_columns=[],
        )


def test_load_colony_rows_rejects_unknown_column_mapping_key():
    with pytest.raises(
        ValueError,
        match=(
            r"Unknown column_mapping key\(s\): animal\. "
            r"Valid keys are: animal_id, sex, cage, mother_id, father_id"
        ),
    ):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            [],
            column_mapping={"animal": "animal_id"},
            gene_columns=["Cre"],
        )


def test_load_colony_rows_rejects_duplicate_header_name():
    with pytest.raises(ValueError, match=r"Duplicate header 'animal_id'"):
        load_colony_rows(
            ["animal_id", "animal_id", "cage", "mother_id", "father_id", "Cre"],
            [],
            gene_columns=["Cre"],
        )


def test_load_colony_rows_ignores_blank_header_name():
    rows = [
        {
            "animal_id": "M01",
            "sex": "M",
            "cage": "C1",
            "mother_id": "",
            "father_id": "",
            "Cre": "Cre/+",
            "": "ignored",
        }
    ]

    colony = load_colony_rows(
        ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre", ""],
        rows,
        gene_columns=["Cre"],
    )

    assert len(colony) == 1
    assert colony[0].animal_id == "M01"
    assert colony[0].genotypes["Cre"] == "Cre/+"


def test_load_colony_accepts_trailing_commas_in_header_and_data_rows(tmp_path):
    csv_path = tmp_path / "trailing_commas.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre,\n"
        "M01,M,C1,,,Cre/+, \n"
        "F01,F,C2,,,+/+, \n",
        encoding="utf-8",
    )

    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert [animal.animal_id for animal in colony] == ["M01", "F01"]
    assert colony[0].genotypes["Cre"] == "Cre/+"
    assert colony[1].genotypes["Cre"] == "+/+"


def test_load_colony_rows_requires_row_keys_to_match_fieldnames():
    rows = [
        {
            "animal_id": "M01",
            "sex": "M",
            "cage": "C1",
            "mother_id": "",
            "Cre": "Cre/+",
        }
    ]

    with pytest.raises(ValueError, match=r"Row 2 keys must match the header names.*father_id"):
        load_colony_rows(
            ["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"],
            rows,
            gene_columns=["Cre"],
        )
