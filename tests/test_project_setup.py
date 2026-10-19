"""Smoke tests: the package imports and the sample data is well-formed."""
import csv
from pathlib import Path

import colony_breeder

SAMPLE = Path(__file__).resolve().parent.parent / "data" / "sample_mice.csv"


def test_package_has_version():
    assert colony_breeder.__version__ == "0.1.0"


def test_sample_csv_parents_exist_and_have_right_sex():
    with open(SAMPLE, newline="") as f:
        rows = list(csv.DictReader(f))
    by_id = {r["animal_id"]: r for r in rows}
    assert len(by_id) == len(rows), "animal_id values must be unique"
    for r in rows:
        if r["mother_id"]:
            assert by_id[r["mother_id"]]["sex"] == "F"
        if r["father_id"]:
            assert by_id[r["father_id"]]["sex"] == "M"
