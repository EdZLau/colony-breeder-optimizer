import codecs
from pathlib import Path

import pytest

from colony_breeder.core.colony import load_colony, pedigree_warnings

SAMPLE = Path(__file__).resolve().parent.parent / "data" / "sample_mice.csv"


def test_load_colony_uses_default_column_mapping_and_reads_sample_data():
    colony = load_colony(str(SAMPLE), gene_columns=["Cre", "Ob", "Agouti"])

    assert len(colony) == 12

    mouse = next(item for item in colony if item.animal_id == "M03")
    assert mouse.sex == "M"
    assert mouse.cage == "C4"
    assert mouse.mother_id == "F01"
    assert mouse.father_id == "M01"
    assert mouse.genotypes["Cre"] == "Cre/+"
    assert mouse.genotypes["Ob"] == "ob/ob"
    assert mouse.genotypes["Agouti"] == "Ay/a"


def test_load_colony_raises_for_duplicate_ids(tmp_path):
    csv_path = tmp_path / "duplicate.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        "M01,F,C2,,,Cre/+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Duplicate animal_id.*M01"):
        load_colony(str(csv_path), gene_columns=["Cre"])


def test_load_colony_rejects_missing_required_columns(tmp_path):
    csv_path = tmp_path / "missing_columns.csv"
    csv_path.write_text(
        "animal_id,sex,cage,Cre\n"
        "M01,M,C1,Cre/+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Required columns missing"):
        load_colony(str(csv_path), gene_columns=["Cre"])


def test_load_colony_supports_non_default_column_mapping(tmp_path):
    csv_path = tmp_path / "mapped_headers.csv"
    csv_path.write_text(
        "Mouse ID,Gender,Cage No,Dam,Sire,Cre\n"
        "M01,M,C1,F02,M02,Cre/+\n",
        encoding="utf-8",
    )
    mapping = {
        "animal_id": "Mouse ID",
        "sex": "Gender",
        "cage": "Cage No",
        "mother_id": "Dam",
        "father_id": "Sire",
    }

    colony = load_colony(csv_path, column_mapping=mapping, gene_columns=["Cre"])

    assert len(colony) == 1
    assert colony[0].animal_id == "M01"
    assert colony[0].sex == "M"
    assert colony[0].cage == "C1"
    assert colony[0].mother_id == "F02"
    assert colony[0].father_id == "M02"
    assert colony[0].genotypes["Cre"] == "Cre/+"


def test_load_colony_names_missing_mapped_header(tmp_path):
    csv_path = tmp_path / "missing_mapped_header.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,F02,M02,Cre/+\n",
        encoding="utf-8",
    )
    mapping = {"animal_id": "Mouse ID"}

    with pytest.raises(ValueError, match="Required columns missing.*Mouse ID"):
        load_colony(csv_path, column_mapping=mapping, gene_columns=["Cre"])


def test_load_colony_reports_missing_gene_column_and_available_columns(tmp_path):
    csv_path = tmp_path / "missing_gene.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Genotype\n"
        "M01,M,C1,,,Cre/+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"Cre.*Available columns:.*Genotype"):
        load_colony(csv_path, gene_columns=["Cre"])


def test_load_colony_rejects_invalid_or_blank_sex(tmp_path):
    csv_path = tmp_path / "bad_sex.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,?,C1,,,Cre/+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Invalid sex.*M01"):
        load_colony(str(csv_path), gene_columns=["Cre"])

    blank_path = tmp_path / "blank_sex.csv"
    blank_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M02,,C1,,,Cre/+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="Invalid sex.*M02"):
        load_colony(str(blank_path), gene_columns=["Cre"])


def test_load_colony_normalises_sex_and_trims_whitespace_in_ids_and_genotypes(tmp_path):
    csv_path = tmp_path / "normalise.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "  M04  ,  m  ,  C9  ,  F02  ,  M02  ,  Cre/+  \n",
        encoding="utf-8",
    )

    colony = load_colony(str(csv_path), gene_columns=["Cre"])
    assert len(colony) == 1
    animal = colony[0]
    assert animal.animal_id == "M04"
    assert animal.sex == "M"
    assert animal.cage == "C9"
    assert animal.mother_id == "F02"
    assert animal.father_id == "M02"
    assert animal.genotypes["Cre"] == "Cre/+"


def test_load_colony_validates_all_non_blank_genotypes_and_names_the_animal_and_column(tmp_path):
    csv_path = tmp_path / "bad_genotype.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,badvalue\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="M01.*Cre"):
        load_colony(str(csv_path), gene_columns=["Cre"])


def test_load_colony_includes_genotype_parser_error_in_message(tmp_path):
    csv_path = tmp_path / "bad_genotype_detail.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match=r"Invalid genotype.*Cannot read genotype 'Cre\+'"):
        load_colony(csv_path, gene_columns=["Cre"])


def test_load_colony_keeps_blank_genotypes_as_unknown_values(tmp_path):
    csv_path = tmp_path / "blank_genotype.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre,Ob\n"
        "M01,M,C1,,, ,\n",
        encoding="utf-8",
    )

    colony = load_colony(str(csv_path), gene_columns=["Cre", "Ob"])
    mouse = colony[0]
    assert mouse.genotypes["Cre"] == ""
    assert mouse.genotypes["Ob"] == ""


def test_load_colony_skips_rows_where_every_field_is_blank(tmp_path):
    csv_path = tmp_path / "trailing_blank_row.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        ",,,,,\n",
        encoding="utf-8",
    )

    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert [animal.animal_id for animal in colony] == ["M01"]


def test_load_colony_skips_rows_with_only_whitespace_and_commas(tmp_path):
    csv_path = tmp_path / "whitespace_blank_row.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        " , , , , ,\n",
        encoding="utf-8",
    )

    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert [animal.animal_id for animal in colony] == ["M01"]


def test_load_colony_rejects_duplicate_gene_column_names(tmp_path):
    csv_path = tmp_path / "duplicate_gene_names.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="duplicate gene name.*Cre"):
        load_colony(csv_path, gene_columns=["Cre", " Cre "])


def test_load_colony_accepts_bom_prefixed_files(tmp_path):
    csv_path = tmp_path / "bom.csv"
    content = (
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
    )
    csv_path.write_bytes(codecs.BOM_UTF8 + content.encode("utf-8"))

    colony = load_colony(str(csv_path), gene_columns=["Cre"])
    assert len(colony) == 1
    assert colony[0].animal_id == "M01"


def test_pedigree_warnings_detects_animal_listed_as_own_parent(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "self_parent.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,M01,,Cre/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert pedigree_warnings(colony) == [
        ("warning", "Animal 'M01' lists itself in mother_id."),
        (
            "warning",
            "Animal 'M01' has mother_id 'M01'; the recorded parent sex is 'M', expected 'F'.",
        ),
    ]


def test_pedigree_warnings_detects_same_id_as_both_parents(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "same_parents.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "P01,M,C1,,,Cre/+\n"
        "M01,M,C2,P01,P01,Cre/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert pedigree_warnings(colony) == [
        (
            "warning",
            "Animal 'M01' has the same ID in mother_id and father_id.",
        ),
        (
            "warning",
            "Animal 'M01' has mother_id 'P01'; the recorded parent sex is 'M', expected 'F'.",
        ),
    ]


def test_pedigree_warnings_detects_mother_with_non_female_sex(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "bad_mother_sex.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "P01,M,C1,,,Cre/+\n"
        "F01,F,C2,P01,,Cre/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert pedigree_warnings(colony) == [
        (
            "warning",
            "Animal 'F01' has mother_id 'P01'; the recorded parent sex is 'M', expected 'F'.",
        )
    ]


def test_pedigree_warnings_detects_father_with_non_male_sex(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "bad_father_sex.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "P01,F,C1,,,Cre/+\n"
        "F01,F,C2,,P01,Cre/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert pedigree_warnings(colony) == [
        (
            "warning",
            "Animal 'F01' has father_id 'P01'; the recorded parent sex is 'F', expected 'M'.",
        )
    ]


def test_pedigree_warnings_detects_two_animal_parent_loop(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "parent_loop.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,F01,,Cre/+\n"
        "F01,F,C2,,M01,Cre/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert pedigree_warnings(colony) == [
        (
            "warning",
            "Animals 'F01' and 'M01' form a two-animal parent loop; "
            "check 'M01'.mother_id and 'F01'.father_id.",
        )
    ]


def test_pedigree_warnings_reports_parent_not_in_colony_as_info(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "missing_parent.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,F00,M00,Cre/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert pedigree_warnings(colony) == [
        (
            "info",
            "Animal 'M01' has mother_id 'F00', which is not in the colony.",
        ),
        (
            "info",
            "Animal 'M01' has father_id 'M00', which is not in the colony.",
        ),
    ]


def test_pedigree_warnings_reports_no_issues_for_sample_colony():
    colony = load_colony(SAMPLE, gene_columns=["Cre", "Ob", "Agouti"])

    assert pedigree_warnings(colony) == []


def test_spaces_only_parent_and_cage_cells_are_blank_for_pedigree_warnings(tmp_path):
    csv_path = _write_pedigree_csv(
        tmp_path,
        "spaces_are_blank.csv",
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,   ,   ,   ,Cre/+\n"
        "F01,F,,,,+/+\n",
    )
    colony = load_colony(csv_path, gene_columns=["Cre"])

    assert colony[0].cage == ""
    assert colony[0].mother_id == ""
    assert colony[0].father_id == ""
    assert pedigree_warnings(colony) == []


def _write_pedigree_csv(tmp_path, filename, contents):
    csv_path = tmp_path / filename
    csv_path.write_text(contents, encoding="utf-8")
    return csv_path
