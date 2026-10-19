"""Colony loading and validation for CSV-backed breeding data."""

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from colony_breeder.core.genetics import (
    format_x_linked_genotype,
    parse_genotype,
    parse_x_linked_genotype,
)
from colony_breeder.core.linkage import (
    RecombinationInput,
    linked_offspring_distribution,
    parse_phased_genotype,
    validate_phase_matches_genotypes,
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
    phases: Dict[Tuple[str, str], str] = field(default_factory=dict)


def _normalise_mapping(column_mapping: Optional[Mapping[str, str]] = None) -> Dict[str, str]:
    """Merge user overrides into the default column mapping."""
    resolved = DEFAULT_COLUMN_MAPPING.copy()
    if column_mapping:
        unknown_keys = set(column_mapping).difference(DEFAULT_COLUMN_MAPPING)
        if unknown_keys:
            valid_keys = ", ".join(DEFAULT_COLUMN_MAPPING)
            raise ValueError(
                "Unknown column_mapping key(s): "
                + ", ".join(sorted(unknown_keys))
                + f". Valid keys are: {valid_keys}."
            )
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
    linked_pairs: Optional[
        Sequence[Tuple[str, str, RecombinationInput, str]]
    ] = None,
) -> List[Animal]:
    """Load a colony CSV and validate animal IDs, sex values, and genotype strings.

    The default column mapping matches the sample data in this project:
    animal_id, sex, cage, mother_id, father_id.
    Sex is normalized by trimming whitespace and uppercasing; only "M" and
    "F" are usable for validating an X-linked genotype.
    Each linked pair is configured as (gene1, gene2, recombination frequency,
    phase column). Phase text follows that gene order, for example a pair
    configured as ("B", "A", ...) uses "B|A/b|a".
    """
    csv_path = Path(csv_path)

    with open(csv_path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is not None:
            reader.fieldnames = [
                field_name.strip() if field_name is not None else ""
                for field_name in reader.fieldnames
            ]
        return load_colony_rows(
            reader.fieldnames,
            reader,
            column_mapping=column_mapping,
            gene_columns=gene_columns,
            x_linked_genes=x_linked_genes,
            linked_pairs=linked_pairs,
            source_name=str(csv_path),
        )


def load_colony_rows(
    fieldnames: Optional[Sequence[Optional[str]]],
    rows: Iterable[Mapping[Optional[str], object]],
    column_mapping: Optional[Mapping[str, str]] = None,
    gene_columns: Optional[Sequence[str]] = None,
    x_linked_genes: Optional[Sequence[str]] = None,
    linked_pairs: Optional[
        Sequence[Tuple[str, str, RecombinationInput, str]]
    ] = None,
    source_name: str = "uploaded data",
    first_data_row: int = 2,
) -> List[Animal]:
    """Validate and load colony rows using their CSV-style field names."""
    if fieldnames is None:
        raise ValueError(
            f"CSV file {source_name} is empty or missing a header row."
        )

    seen_fieldnames = set()
    for column_index, fieldname in enumerate(fieldnames, start=1):
        if not isinstance(fieldname, str):
            raise ValueError(
                f"Blank header at column {column_index} in {source_name}."
            )
        if not fieldname.strip():
            continue
        if fieldname in seen_fieldnames:
            raise ValueError(
                f"Duplicate header {fieldname!r} in {source_name}."
            )
        seen_fieldnames.add(fieldname)

    resolved_mapping = _normalise_mapping(column_mapping)
    gene_names = _validate_gene_columns(gene_columns)
    x_linked_gene_names = set(_validate_gene_columns(x_linked_genes or []))
    unknown_x_linked_genes = x_linked_gene_names.difference(gene_names)
    if unknown_x_linked_genes:
        raise ValueError(
            "X-linked gene name(s) must also be included in gene_columns: "
            + ", ".join(sorted(unknown_x_linked_genes))
        )
    normalized_linked_pairs = _validate_linked_pairs(
        linked_pairs or (),
        gene_names,
        x_linked_gene_names,
    )

    missing_columns = []
    for canonical_name in DEFAULT_COLUMN_MAPPING:
        actual_name = resolved_mapping.get(canonical_name, DEFAULT_COLUMN_MAPPING[canonical_name])
        if actual_name not in fieldnames:
            missing_columns.append(actual_name)
    if missing_columns:
        raise ValueError(
            "Required columns missing from the CSV: " + ", ".join(sorted(set(missing_columns)))
        )

    missing_gene_columns = [gene for gene in gene_names if gene not in fieldnames]
    if missing_gene_columns:
        raise ValueError(
            "Requested gene column(s) missing from the CSV: "
            + ", ".join(missing_gene_columns)
            + ". Available columns: "
            + ", ".join(fieldname or "" for fieldname in fieldnames)
        )

    missing_phase_columns = [
        phase_column
        for _, _, _, phase_column in normalized_linked_pairs
        if phase_column not in fieldnames
    ]
    if missing_phase_columns:
        raise ValueError(
            "Linked-pair phase column(s) missing from the CSV: "
            + ", ".join(dict.fromkeys(missing_phase_columns))
            + ". Available columns: "
            + ", ".join(fieldname or "" for fieldname in fieldnames)
        )

    animals: List[Animal] = []
    seen_ids = set()

    for row_index, row in enumerate(rows, start=first_data_row):
        for column_name, value in row.items():
            if column_name is None:
                extra_values = value if isinstance(value, list) else []
                for extra_value in extra_values:
                    if extra_value is not None and not isinstance(extra_value, str):
                        raise ValueError(
                            f"Row {row_index} has a non-text value in column "
                            "'<extra>'. Expected text or None."
                        )
                continue
            if value is not None and not isinstance(value, str):
                raise ValueError(
                    f"Row {row_index} has a non-text value in column "
                    f"{column_name!r}. Expected text or None."
                )

        extra_values = row.get(None, [])
        if isinstance(extra_values, list) and any(
            value.strip() for value in extra_values if value is not None
        ):
            raise ValueError(f"Row {row_index} has more values than headers.")

        row_keys = {
            key for key in row
            if isinstance(key, str) and key.strip()
        }
        expected_keys = {
            fieldname for fieldname in fieldnames
            if isinstance(fieldname, str) and fieldname.strip()
        }
        if row_keys != expected_keys:
            missing_keys = sorted(expected_keys - row_keys)
            extra_keys = sorted(
                (key for key in row_keys - expected_keys),
                key=lambda key: "" if key is None else str(key),
            )
            raise ValueError(
                f"Row {row_index} keys must match the header names. "
                f"Missing: {missing_keys!r}; unexpected: {extra_keys!r}."
            )

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
            raise ValueError(f"Duplicate animal_id {raw_id!r} found in {source_name}.")
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

        phases: Dict[Tuple[str, str], str] = {}
        for gene1, gene2, _, phase_column in normalized_linked_pairs:
            raw_phase = row.get(phase_column, "")
            phase = (raw_phase or "").strip()
            if phase and genotypes[gene1] and genotypes[gene2]:
                try:
                    parse_phased_genotype(phase)
                    validate_phase_matches_genotypes(
                        phase,
                        genotypes[gene1],
                        genotypes[gene2],
                    )
                except ValueError as exc:
                    raise ValueError(
                        f"Invalid phase for animal {raw_id!r} in gene pair "
                        f"({gene1!r}, {gene2!r}): {phase!r}. {exc}"
                    ) from exc
            phases[(gene1, gene2)] = phase

        animals.append(
            Animal(
                animal_id=raw_id,
                sex=sex,
                cage=cage,
                mother_id=mother_id,
                father_id=father_id,
                genotypes=genotypes,
                phases=phases,
            )
        )

    return animals


def _validate_linked_pairs(
    linked_pairs: Sequence[Tuple[str, str, RecombinationInput, str]],
    gene_names: Sequence[str],
    x_linked_gene_names: set,
) -> List[Tuple[str, str, RecombinationInput, str]]:
    """Validate linked-pair settings and normalize their names."""
    normalized = []
    linked_genes = set()
    known_genes = set(gene_names)
    for pair in linked_pairs:
        try:
            gene1, gene2, recombination_frequency, phase_column = pair
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Each linked pair must contain gene1, gene2, recombination "
                "frequency, and phase column."
            ) from exc

        gene1 = str(gene1).strip()
        gene2 = str(gene2).strip()
        phase_column = str(phase_column).strip()
        pair_label = f"({gene1!r}, {gene2!r})"
        if not gene1 or not gene2:
            raise ValueError(f"Linked pair {pair_label} must name two genes.")
        if gene1 == gene2:
            raise ValueError(
                f"Linked pair {pair_label} cannot contain the same gene twice."
            )
        if gene1 not in known_genes or gene2 not in known_genes:
            unknown = sorted({gene1, gene2} - known_genes)
            raise ValueError(
                f"Linked pair {pair_label} names gene(s) not in gene_columns: "
                + ", ".join(unknown)
            )
        if gene1 in x_linked_gene_names or gene2 in x_linked_gene_names:
            raise ValueError(
                f"Linked pair {pair_label} cannot include an X-linked gene."
            )
        duplicated = {gene1, gene2} & linked_genes
        if duplicated:
            raise ValueError(
                "A gene cannot appear in more than one linked pair: "
                + ", ".join(sorted(duplicated))
            )
        if not phase_column:
            raise ValueError(f"Linked pair {pair_label} needs a phase column.")

        try:
            linked_offspring_distribution(
                "A|B/A|B",
                "A|B/A|B",
                recombination_frequency,
            )
        except ValueError as exc:
            raise ValueError(
                f"Invalid recombination frequency for linked pair {pair_label}: {exc}"
            ) from exc

        normalized.append(
            (gene1, gene2, recombination_frequency, phase_column)
        )
        linked_genes.update((gene1, gene2))
    return normalized


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
