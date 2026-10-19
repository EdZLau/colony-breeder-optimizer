"""Colony loading and validation for CSV-backed breeding data."""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Union

from colony_breeder.core.genetics import parse_genotype

DEFAULT_COLUMN_MAPPING = {
    "animal_id": "animal_id",
    "sex": "sex",
    "cage": "cage",
    "mother_id": "mother_id",
    "father_id": "father_id",
}


@dataclass(frozen=True)
class Animal:
    """A colony animal with identifiers and gene genotypes."""

    animal_id: str
    sex: str
    cage: str
    mother_id: str
    father_id: str
    genotypes: Dict[str, str] = field(default_factory=dict)


def _normalise_mapping(column_mapping: Optional[Mapping[str, str]] = None) -> Dict[str, str]:
    """Merge user overrides into the default column mapping."""
    resolved = DEFAULT_COLUMN_MAPPING.copy()
    if column_mapping:
        resolved.update(dict(column_mapping))
    return resolved


def _validate_gene_columns(gene_columns: Optional[Sequence[str]]) -> List[str]:
    if gene_columns is None:
        raise ValueError("gene_columns is required; provide the columns that hold genotype data.")

    columns = []
    for gene in gene_columns:
        text = str(gene).strip()
        if not text:
            raise ValueError("gene_columns contains an empty gene name.")
        columns.append(text)
    return columns


def load_colony(
    csv_path: Union[str, Path],
    column_mapping: Optional[Mapping[str, str]] = None,
    gene_columns: Optional[Sequence[str]] = None,
) -> List[Animal]:
    """Load a colony CSV and validate animal IDs, sex values, and genotype strings.

    The default column mapping matches the sample data in this project:
    animal_id, sex, cage, mother_id, father_id.
    """
    resolved_mapping = _normalise_mapping(column_mapping)
    gene_names = _validate_gene_columns(gene_columns)
    csv_path = Path(csv_path)

    with open(csv_path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError(f"CSV file {csv_path} is empty or missing a header row.")

        cleaned_fieldnames = [field_name.strip() if field_name is not None else "" for field_name in reader.fieldnames]
        reader.fieldnames = cleaned_fieldnames

        missing_columns = []
        for canonical_name in DEFAULT_COLUMN_MAPPING:
            actual_name = resolved_mapping.get(canonical_name, DEFAULT_COLUMN_MAPPING[canonical_name])
            if actual_name not in reader.fieldnames:
                missing_columns.append(actual_name)
        if missing_columns:
            raise ValueError(
                "Required columns missing from the CSV: " + ", ".join(sorted(set(missing_columns)))
            )

        animals: List[Animal] = []
        seen_ids = set()

        for row_index, row in enumerate(reader, start=2):
            raw_id = (row.get(resolved_mapping["animal_id"], "") or "").strip()
            if not raw_id:
                raise ValueError(f"Row {row_index} is missing an animal_id.")
            if raw_id in seen_ids:
                raise ValueError(f"Duplicate animal_id {raw_id!r} found in {csv_path}.")
            seen_ids.add(raw_id)

            sex_value = (row.get(resolved_mapping["sex"], "") or "").strip()
            sex = sex_value.upper()
            if sex not in {"M", "F"}:
                raise ValueError(f"Invalid sex value {sex_value!r} for animal {raw_id!r}. Expected 'M' or 'F'.")

            cage = (row.get(resolved_mapping["cage"], "") or "").strip()
            mother_id = (row.get(resolved_mapping["mother_id"], "") or "").strip()
            father_id = (row.get(resolved_mapping["father_id"], "") or "").strip()

            genotypes: Dict[str, str] = {}
            for gene_name in gene_names:
                raw_value = row.get(gene_name, "")
                if raw_value is None:
                    raw_value = ""
                cleaned = raw_value.strip()
                if cleaned:
                    try:
                        parse_genotype(cleaned)
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid genotype for animal {raw_id!r} in column {gene_name!r}: {cleaned!r}"
                        ) from exc
                genotypes[gene_name] = cleaned

            animals.append(
                Animal(
                    animal_id=raw_id,
                    sex=sex,
                    cage=cage,
                    mother_id=mother_id,
                    father_id=father_id,
                    genotypes=genotypes,
                )
            )

    return animals
