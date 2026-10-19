import csv
import io
import time as time_module
import zipfile
from datetime import date, timedelta

import pytest

openpyxl = pytest.importorskip("openpyxl")

from colony_breeder.app import upload
from colony_breeder.app.upload import SheetChoiceRequired, read_upload
from colony_breeder.core.colony import load_colony_rows


def _workbook_bytes(build):
    workbook = openpyxl.Workbook()
    build(workbook)
    stream = io.BytesIO()
    workbook.save(stream)
    return stream.getvalue()


def test_xlsx_values_are_text_and_trailing_blank_rows_and_columns_are_trimmed():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "number", "leading_zero", "genotype", None, None])
        sheet.append(["M01", 3, "001", "+/+", None, None])
        sheet.append([None, None, None, None, None, None])
        sheet.append(["F01", 3.0, "002", "+/+", None, None])
        sheet.append([None, None, None, None, None, None])

    result = read_upload(_workbook_bytes(build), "records.csv")

    assert result.fieldnames == ["id", "number", "leading_zero", "genotype"]
    assert result.rows == [
        {
            "id": "M01",
            "number": "3",
            "leading_zero": "001",
            "genotype": "+/+",
        },
        {
            "id": "",
            "number": "",
            "leading_zero": "",
            "genotype": "",
        },
        {
            "id": "F01",
            "number": "3",
            "leading_zero": "002",
            "genotype": "+/+",
        },
    ]
    assert result.sheet_name == "Sheet"
    assert result.warnings == []


def test_xlsx_date_warning_is_grouped_once_per_column():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "date"])
        sheet.append(["M01", date(2020, 1, 2)])
        sheet.append(["F01", date(2020, 1, 3)])

    result = read_upload(_workbook_bytes(build), "dates.xlsx")

    assert len(result.warnings) == 1
    assert "Column 'date'" in result.warnings[0]
    assert "2 date/time-typed cell(s)" in result.warnings[0]
    assert "first on row 2" in result.warnings[0]
    assert "not repaired" in result.warnings[0]
    assert result.warnings[0].column == "date"
    assert result.rows[0]["date"].startswith("2020-01-02")


def test_xlsx_timedelta_cells_are_included_in_date_warning():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "duration"])
        sheet.append(["M01", timedelta(hours=3)])

    result = read_upload(_workbook_bytes(build), "durations.xlsx")

    assert len(result.warnings) == 1
    assert "Column 'duration'" in result.warnings[0]
    assert "1 date/time-typed cell(s)" in result.warnings[0]


def test_xlsx_boolean_cells_are_returned_as_text_with_warning():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "flag"])
        sheet.append(["M01", True])

    result = read_upload(_workbook_bytes(build), "boolean.xlsx")

    assert result.rows[0]["flag"] == "True"
    assert any("Column 'flag' contains 1 boolean cell(s)" in warning for warning in result.warnings)


def test_xlsx_multiple_sheets_require_a_choice():
    def build(workbook):
        workbook.active.title = "Mice"
        workbook.active.append(["id"])
        workbook.create_sheet("Rats").append(["id"])

    content = _workbook_bytes(build)

    with pytest.raises(SheetChoiceRequired) as exc_info:
        read_upload(content, "colony.xlsx")

    assert exc_info.value.sheet_names == ("Mice", "Rats")
    result = read_upload(content, "colony.xlsx", sheet_name="Rats")
    assert result.sheet_name == "Rats"
    assert result.fieldnames == ["id"]


def test_xlsx_header_row_is_one_based_and_preserves_interior_blank_rows():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["instructions"])
        sheet.append(["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"])
        sheet.append(["M01", "M", "C1", "", "", "Cre/+"])
        sheet.append([None, None, None, None, None, None])
        sheet.append([None, "F", "C2", "", "", "+/+"])

    result = read_upload(
        _workbook_bytes(build),
        "header.xlsx",
        header_row=2,
    )

    assert result.first_data_row == 3
    assert len(result.rows) == 3
    with pytest.raises(ValueError, match=r"Row 5 is missing an animal_id"):
        load_colony_rows(
            result.fieldnames,
            result.rows,
            gene_columns=["Cre"],
            first_data_row=result.first_data_row,
        )


def test_xlsx_upload_to_colony_rows_preserves_sheet_row_numbers():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["Instructions"])
        sheet.append(["For fake test data only"])
        sheet.append(["animal_id", "sex", "cage", "mother_id", "father_id", "Cre"])
        sheet.append(["M01", "M", "C1", "", "", "Cre/+"])
        sheet.append(["", "F", "C2", "", "", "+/+"])

    result = read_upload(
        _workbook_bytes(build),
        "colony.xlsx",
        header_row=3,
    )

    with pytest.raises(ValueError, match=r"Row 5 is missing an animal_id"):
        load_colony_rows(
            result.fieldnames,
            result.rows,
            gene_columns=["Cre"],
            first_data_row=result.first_data_row,
        )


def test_xlsx_streaming_trims_styled_stray_cell_at_far_bottom_right():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "genotype"])
        sheet.append(["M01", "A/a"])
        sheet["Z100"].number_format = "0.00"

    result = read_upload(_workbook_bytes(build), "styled-cell.xlsx")

    assert result.fieldnames == ["id", "genotype"]
    assert result.rows == [{"id": "M01", "genotype": "A/a"}]


def test_xlsx_styled_cell_at_sheet_limit_is_bounded_and_warns():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "genotype"])
        sheet.append(["M01", "A/a"])
        sheet.cell(row=1_048_576, column=16_384).number_format = "0.00"

    content = _workbook_bytes(build)
    started_at = time_module.monotonic()
    result = read_upload(content, "last-cell.xlsx")
    elapsed = time_module.monotonic() - started_at

    assert elapsed < 10
    assert result.fieldnames == ["id", "genotype"]
    assert result.rows == [{"id": "M01", "genotype": "A/a"}]
    assert any("Data below worksheet row" in warning for warning in result.warnings)


def test_xlsx_column_limit_is_checked_after_trailing_blank_columns_are_trimmed(
    monkeypatch,
):
    monkeypatch.setattr(upload, "COLUMN_LIMIT", 2)

    valid = _workbook_bytes(
        lambda workbook: workbook.active.append(["id", "value", None, None])
    )
    assert read_upload(valid, "trailing-columns.xlsx").fieldnames == ["id", "value"]

    invalid = _workbook_bytes(
        lambda workbook: workbook.active.append(["id", "value", "extra"])
    )
    with pytest.raises(ValueError, match="2-column upload limit"):
        read_upload(invalid, "too-wide.xlsx")


def test_upload_byte_size_limit_rejects_large_input(monkeypatch):
    monkeypatch.setattr(upload, "UPLOAD_BYTE_LIMIT", 1024 * 1024)

    with pytest.raises(ValueError, match="1 MB file-size limit"):
        read_upload(b"x" * (1024 * 1024 + 1), "large.csv")


def test_xlsx_expanded_size_limit_rejects_large_archive(monkeypatch):
    monkeypatch.setattr(upload, "XLSX_EXPANDED_BYTE_LIMIT", 1)
    content = _workbook_bytes(lambda workbook: None)

    with pytest.raises(ValueError, match="expanded-size limit"):
        read_upload(content, "expanded.xlsx")


def test_xlsx_extra_cells_use_none_key_and_short_rows_are_padded():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id", "value"])
        sheet.append(["M01"])
        sheet.append(["F01", "text", "extra", ""])

    result = read_upload(_workbook_bytes(build), "cells.xlsx")

    assert result.rows[0] == {"id": "M01", "value": ""}
    assert result.rows[1] == {"id": "F01", "value": "text", None: ["extra"]}


@pytest.mark.parametrize(
    ("content", "expected_delimiter"),
    [
        (b"id,genotype\nM01,+/+\n", ","),
        (b"id;genotype\nM01;+/+\n", ";"),
        (b"id\tgenotype\nM01\t+/+\n", "\t"),
    ],
)
def test_csv_delimiter_detection_and_genotype_preservation(content, expected_delimiter):
    result = read_upload(content, "data.xlsx")

    assert result.delimiter == expected_delimiter
    assert result.fieldnames == ["id", "genotype"]
    assert result.rows[0]["genotype"] == "+/+"
    assert result.encoding == "utf-8-sig"


def test_csv_utf8_bom_and_stripped_headers():
    result = read_upload(
        b"\xef\xbb\xbf id , genotype \nM01,+/+\n",
        "utf8.csv",
    )

    assert result.fieldnames == ["id", "genotype"]
    assert result.rows[0]["genotype"] == "+/+"
    assert result.encoding == "utf-8-sig"


def test_csv_utf16_bom_is_detected():
    content = "id,genotype\nM01,+/+\n".encode("utf-16")

    result = read_upload(content, "utf16.csv")

    assert result.encoding == "utf-16"
    assert result.fieldnames == ["id", "genotype"]
    assert result.rows[0]["genotype"] == "+/+"


def test_csv_cp1252_fallback_reports_non_utf8_warning():
    content = "id,note\nM01,café\n".encode("cp1252")

    result = read_upload(content, "legacy.csv")

    assert result.encoding == "cp1252"
    assert result.rows[0]["note"] == "café"
    assert any("not UTF-8" in warning for warning in result.warnings)


def test_csv_decode_failure_has_clear_value_error():
    with pytest.raises(ValueError, match="neither valid UTF-8 nor cp1252"):
        read_upload(b"id,note\nM01,\x81\n", "bad.csv")


def test_csv_keeps_interior_blank_rows_and_trims_trailing_blank_rows_and_columns():
    result = read_upload(
        b"id,genotype,,\nM01,+/+\n\nF01,A/a,,\n\n",
        "rows.csv",
    )

    assert result.fieldnames == ["id", "genotype"]
    assert len(result.rows) == 3
    assert result.rows[1] == {"id": "", "genotype": ""}
    assert result.rows[2]["id"] == "F01"


def test_csv_pads_short_rows_and_places_extra_cells_under_none():
    result = read_upload(
        b"id,value\nM01\nF01,text,extra,\n",
        "width.csv",
    )

    assert result.rows == [
        {"id": "M01", "value": ""},
        {"id": "F01", "value": "text", None: ["extra"]},
    ]


def test_csv_header_row_and_first_data_row_metadata():
    result = read_upload(
        b"Instructions\nanimal_id,sex,cage,mother_id,father_id,Cre\n"
        b"M01,M,C1,,,Cre/+\n",
        "header.csv",
        header_row=2,
    )

    assert result.first_data_row == 3
    assert result.rows[0]["Cre"] == "Cre/+"


def test_csv_delimiter_detection_ignores_title_line_and_handles_decimal_commas():
    result = read_upload(
        b"Colony upload, not the header;ignore this line\n"
        b"animal_id;score;genotype\n"
        b"M01;1,25;+/+\n"
        b"F01;2,50;A/a\n",
        "semicolon.csv",
        header_row=2,
    )

    assert result.delimiter == ";"
    assert result.fieldnames == ["animal_id", "score", "genotype"]
    assert result.rows[0]["score"] == "1,25"
    assert result.rows[0]["genotype"] == "+/+"


def test_csv_column_limit_rejects_wide_header(monkeypatch):
    monkeypatch.setattr(upload, "COLUMN_LIMIT", 2)

    with pytest.raises(ValueError, match="CSV file exceeds the 2-column upload limit"):
        read_upload(b"id,value,extra\nM01,text,more\n", "too-wide.csv")


def test_empty_upload_has_clear_value_error():
    with pytest.raises(ValueError, match="uploaded file is empty"):
        read_upload(b"", "empty.csv")


def test_corrupt_xlsx_has_clear_value_error():
    with pytest.raises(ValueError, match="valid .xlsx workbook"):
        read_upload(b"PK\x03\x04not a workbook", "broken.xlsx")


def test_corrupt_zip_workbook_has_clear_value_error():
    stream = io.BytesIO()
    with zipfile.ZipFile(stream, "w") as archive:
        archive.writestr("not-a-workbook.txt", "invalid workbook contents")

    with pytest.raises(ValueError, match="valid .xlsx workbook"):
        read_upload(stream.getvalue(), "broken.xlsx")


def test_corrupt_sheet_xml_has_clear_value_error():
    content = _workbook_bytes(
        lambda workbook: workbook.active.append(["id", "value"])
    )
    source = io.BytesIO(content)
    corrupted = io.BytesIO()
    with zipfile.ZipFile(source, "r") as original:
        with zipfile.ZipFile(corrupted, "w") as replacement:
            for member in original.infolist():
                data = original.read(member.filename)
                if member.filename == "xl/worksheets/sheet1.xml":
                    data = b"<worksheet><sheetData><row>"
                replacement.writestr(member, data)

    with pytest.raises(ValueError, match="valid .xlsx workbook"):
        read_upload(corrupted.getvalue(), "bad-sheet.xml.xlsx")


def test_empty_worksheet_has_clear_value_error():
    content = _workbook_bytes(lambda workbook: None)

    with pytest.raises(ValueError, match="empty header row"):
        read_upload(content, "empty_sheet.xlsx")


def test_malformed_csv_has_clear_value_error():
    with pytest.raises(ValueError, match="Could not parse CSV"):
        read_upload(b'id,note\nM01,"unterminated\n', "malformed.csv")


@pytest.mark.parametrize("filename", ["old.xls", "protected.xlsx"])
def test_ole_excel_files_request_plain_xlsx(filename):
    ole_content = bytes.fromhex("D0 CF 11 E0 A1 B1 1A E1") + b"encrypted content"

    with pytest.raises(ValueError, match="old .xls or password-protected"):
        read_upload(ole_content, filename)


def test_row_limit_rejects_more_than_twenty_thousand_rows():
    content = b"id\n" + b"M01\n" * 20_001

    with pytest.raises(ValueError, match="20,000-row upload limit"):
        read_upload(content, "too_many.csv")


def test_xlsx_row_limit_rejects_more_than_twenty_thousand_rows():
    def build(workbook):
        sheet = workbook.active
        sheet.append(["id"])
        for _ in range(20_001):
            sheet.append(["M01"])

    with pytest.raises(ValueError, match="20,000-row upload limit"):
        read_upload(_workbook_bytes(build), "too_many.xlsx")
