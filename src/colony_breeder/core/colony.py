"""Colony loading and validation for CSV-backed breeding data."""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from colony_breeder.core.genetics import (
    format_x_linked_genotype,
    parse_genotype,
    parse_x_linked_genotype,
)

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
        if text in columns:
            raise ValueError(f"gene_columns contains duplicate gene name {text!r}.")
        columns.append(text)
    return columns


def load_colony(
    csv_path: Union[str, Path],
    column_mapping: Optional[Mapping[str, str]] = None,
    gene_columns: Optional[Sequence[str]] = None,
    x_linked_genes: Optional[Sequence[str]] = None,
) -> List[Animal]:
    """Load a colony CSV and validate animal IDs, sex values, and genotype strings.

    The default column mapping matches the sample data in this project:
    animal_id, sex, cage, mother_id, father_id.
    Sex is normalized by trimming whitespace and uppercasing; only "M" and
    "F" are usable for validating an X-linked genotype.
    """
    resolved_mapping = _normalise_mapping(column_mapping)
    gene_names = _validate_gene_columns(gene_columns)
    x_linked_gene_names = set(_validate_gene_columns(x_linked_genes or []))
    unknown_x_linked_genes = x_linked_gene_names.difference(gene_names)
    if unknown_x_linked_genes:
        raise ValueError(
            "X-linked gene name(s) must also be included in gene_columns: "
            + ", ".join(sorted(unknown_x_linked_genes))
        )
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

        missing_gene_columns = [gene for gene in gene_names if gene not in reader.fieldnames]
        if missing_gene_columns:
            raise ValueError(
                "Requested gene column(s) missing from the CSV: "
                + ", ".join(missing_gene_columns)
                + ". Available columns: "
                + ", ".join(reader.fieldnames)
            )

        animals: List[Animal] = []
        seen_ids = set()

        for row_index, row in enumerate(reader, start=2):
            extra_values = row.get(None, [])
            if any(value.strip() for value in extra_values if value is not None):
                raise ValueError(f"Row {row_index} has more values than headers.")

            if not any(
                (value or "").strip()
                for key, value in row.items()
                if key is not None and value is not None
            ):
                continue

            raw_id = (row.get(resolved_mapping["animal_id"], "") or "").strip()
            if not raw_id:
                raise ValueError(f"Row {row_index} is missing an animal_id.")
            if raw_id in seen_ids:
                raise ValueError(f"Duplicate animal_id {raw_id!r} found in {csv_path}.")
            seen_ids.add(raw_id)

            sex_value = (row.get(resolved_mapping["sex"], "") or "").strip()
            sex = sex_value.upper()
            if sex not in {"M", "F"} and not x_linked_gene_names:
                raise ValueError(
                    f"Invalid sex value {sex_value!r} for animal {raw_id!r}. Expected 'M' or 'F'."
                )

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
                        if gene_name in x_linked_gene_names:
                            if sex in {"M", "F"}:
                                parsed = parse_x_linked_genotype(cleaned, sex)
                                cleaned = format_x_linked_genotype(parsed)
                            else:
                                parts = [part.strip() for part in cleaned.split("/")]
                                if len(parts) != 2 or not all(parts):
                                    raise ValueError(
                                        f"Cannot read X-linked genotype {cleaned!r}: "
                                        "expected slash style with two non-empty parts, "
                                        "like 'A/a' or 'A/Y'."
                                    )
                        else:
                            parse_genotype(cleaned)
                    except ValueError as exc:
                        raise ValueError(
                            f"Invalid genotype for animal {raw_id!r} in column {gene_name!r}: "
                            f"{cleaned!r}. {exc}"
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


def pedigree_warnings(colony: Sequence[Animal]) -> List[Tuple[str, str]]:
    """Return (severity, message) entries for detectable pedigree issues."""
    animals_by_id = {animal.animal_id.strip(): animal for animal in colony}
    warnings = []
    reported_loops = set()

    for animal in colony:
        animal_id = animal.animal_id.strip()
        mother_id = animal.mother_id.strip()
        father_id = animal.father_id.strip()

        if mother_id and mother_id == animal_id:
            warnings.append(
                ("warning", f"Animal {animal_id!r} lists itself in mother_id.")
            )
        if father_id and father_id == animal_id:
            warnings.append(
                ("warning", f"Animal {animal_id!r} lists itself in father_id.")
            )
        if mother_id and mother_id == father_id:
            warnings.append(
                (
                    "warning",
                    f"Animal {animal_id!r} has the same ID in mother_id and father_id.",
                )
            )

        for field_name, parent_id, expected_sex in (
            ("mother_id", mother_id, "F"),
            ("father_id", father_id, "M"),
        ):
            if not parent_id:
                continue
            parent = animals_by_id.get(parent_id)
            if parent is None:
                warnings.append(
                    (
                        "info",
                        f"Animal {animal_id!r} has {field_name} {parent_id!r}, "
                        "which is not in the colony.",
                    )
                )
                continue
            if parent_id == animal_id:
                continue
            if parent.sex != expected_sex:
                warnings.append(
                    (
                        "warning",
                        f"Animal {animal_id!r} has {field_name} {parent_id!r}; "
                        f"the recorded parent sex is {parent.sex!r}, expected {expected_sex!r}; "
                        f"check {parent_id!r}.sex or {animal_id!r}.{field_name}.",
                    )
                )

            parent_mother_id = parent.mother_id.strip()
            parent_father_id = parent.father_id.strip()
            reciprocal_fields = []
            if parent_mother_id == animal_id:
                reciprocal_fields.append("mother_id")
            if parent_father_id == animal_id:
                reciprocal_fields.append("father_id")
            if reciprocal_fields:
                loop = tuple(sorted((animal_id, parent_id)))
                if loop not in reported_loops:
                    reported_loops.add(loop)
                    reciprocal_field = reciprocal_fields[0]
                    warnings.append(
                        (
                            "warning",
                            f"Animals {loop[0]!r} and {loop[1]!r} form a two-animal "
                            f"parent loop; check {animal_id!r}.{field_name} and "
                            f"{parent_id!r}.{reciprocal_field}.",
                        )
                    )

    return warnings
