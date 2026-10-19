"""Pure helpers for mapping uploaded colony data before the Streamlit UI."""

from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from colony_breeder.core.colony import Animal, load_colony_rows, pedigree_warnings
from colony_breeder.app.upload import UploadResult

NONE_COLUMN = "(none)"
RESERVED_COLUMNS = {
    "cage": "__colony_breeder_unmapped_cage__",
    "mother_id": "__colony_breeder_unmapped_mother__",
    "father_id": "__colony_breeder_unmapped_father__",
}
MAPPING_KEYS = ("animal_id", "sex", "cage", "mother_id", "father_id")
OPTIONAL_MAPPING_KEYS = ("cage", "mother_id", "father_id")
DISTINCT_VALUE_LIMIT = 50

_HEADER_ALIASES = {
    "animal_id": ("animalid", "mouseid", "id", "animal", "mouse"),
    "sex": ("sex", "gender"),
    "cage": ("cage", "cageno", "cagenumber"),
    "mother_id": ("motherid", "mother", "dam", "damid"),
    "father_id": ("fatherid", "father", "sire", "sireid"),
}


@dataclass
class UploadLoadResult:
    """Loaded animals and summaries for a parsed upload."""

    animals: List[Animal]
    warnings: Dict[str, List[str]]
    unknown_sex_count: int
    unknown_sex_ids: List[str]
    distinct_values: Dict[str, List[str]]
    error_message: Optional[str] = None


def _normalize_header(value: str) -> str:
    return "".join(character.lower() for character in value if character.isalnum())


def guess_column_mapping(
    fieldnames: Sequence[str],
) -> Dict[str, Optional[str]]:
    """Guess standard colony columns from names, ignoring punctuation and case."""
    normalized = {
        fieldname: _normalize_header(fieldname)
        for fieldname in fieldnames
    }
    guesses: Dict[str, Optional[str]] = {}
    for key, aliases in _HEADER_ALIASES.items():
        guesses[key] = next(
            (
                fieldname
                for alias in aliases
                for fieldname in fieldnames
                if normalized[fieldname] == alias
            ),
            NONE_COLUMN if key in OPTIONAL_MAPPING_KEYS else None,
        )
    return guesses


def guess_sex_value_mapping(
    rows: Sequence[Mapping[Optional[str], object]],
    sex_column: str,
) -> Dict[str, str]:
    """Guess M/F for observed sex values; map blank and unknown values to unknown."""
    values = _distinct_raw_values(rows, sex_column)
    result: Dict[str, str] = {}
    for value in values:
        normalized = value.strip().lower()
        if normalized in {"m", "male"}:
            result[value] = "M"
        elif normalized in {"f", "female"}:
            result[value] = "F"
        else:
            result[value] = ""
    return result


def summarize_distinct_values(
    rows: Sequence[Mapping[Optional[str], object]],
    columns: Sequence[str],
    limit: int = DISTINCT_VALUE_LIMIT,
) -> Dict[str, List[str]]:
    """Summarize distinct values in first-seen order, adding a hidden-count marker."""
    summaries: Dict[str, List[str]] = {}
    for column in columns:
        values = _distinct_raw_values(rows, column)
        visible = [
            "(blank)" if value == "" else value
            for value in values[:limit]
        ]
        remaining = len(values) - limit
        if remaining > 0:
            visible.append(f"+{remaining} more")
        summaries[column] = visible
    return summaries


def _distinct_raw_values(
    rows: Sequence[Mapping[Optional[str], object]],
    column: str,
) -> List[str]:
    seen = set()
    values = []
    for row in rows:
        raw_value = row.get(column, "")
        value = "" if raw_value is None else str(raw_value)
        if value not in seen:
            seen.add(value)
            values.append(value)
    return values


def build_loader_args(
    fieldnames: Sequence[str],
    rows: Sequence[Mapping[Optional[str], object]],
    column_mapping: Mapping[str, Optional[str]],
    first_data_row: int,
) -> Dict[str, object]:
    """Validate current mapping selections and build loader input arguments."""
    available = set(fieldnames)
    normalized_mapping: Dict[str, str] = {}
    for key in MAPPING_KEYS:
        selected = column_mapping.get(key)
        if key in OPTIONAL_MAPPING_KEYS and selected in (None, NONE_COLUMN):
            selected = RESERVED_COLUMNS[key]
            if selected in available:
                raise ValueError(
                    f"Reserved column name {selected!r} already exists in the uploaded headers."
                )
        elif selected is None or selected == NONE_COLUMN:
            raise ValueError(f"A column must be selected for {key}.")
        elif selected not in available:
            raise ValueError(
                f"Selected column {selected!r} for {key} is no longer available. "
                "Choose a column from the current upload."
            )
        normalized_mapping[key] = selected

    actual_selections = [
        normalized_mapping[key]
        for key in MAPPING_KEYS
        if normalized_mapping[key] not in RESERVED_COLUMNS.values()
    ]
    if len(actual_selections) != len(set(actual_selections)):
        raise ValueError("Each mapped field must use a different uploaded column.")

    loader_fieldnames = list(fieldnames)
    loader_rows = [dict(row) for row in rows]
    for key in OPTIONAL_MAPPING_KEYS:
        reserved_name = RESERVED_COLUMNS[key]
        if normalized_mapping[key] == reserved_name:
            loader_fieldnames.append(reserved_name)
            for row in loader_rows:
                row[reserved_name] = ""

    return {
        "fieldnames": loader_fieldnames,
        "rows": loader_rows,
        "column_mapping": normalized_mapping,
        "gene_columns": [],
        "first_data_row": first_data_row,
        "allow_unknown_sex": True,
    }


def load_from_upload(
    upload: UploadResult,
    column_mapping: Mapping[str, Optional[str]],
    gene_columns: Sequence[str],
    sex_value_mapping: Mapping[str, str],
    source_name: str = "uploaded data",
) -> UploadLoadResult:
    """Validate and load a colony from parsed upload rows without UI dependencies."""
    grouped_warnings = _group_upload_warnings(upload.warnings)
    distinct_values = summarize_distinct_values(upload.rows, gene_columns)
    try:
        args = build_loader_args(
            upload.fieldnames,
            upload.rows,
            column_mapping,
            upload.first_data_row,
        )
        _validate_gene_selections(
            upload.fieldnames,
            gene_columns,
            args["column_mapping"],
        )
        transformed_rows = _apply_sex_values(
            args["rows"],
            args["column_mapping"]["sex"],
            sex_value_mapping,
        )
        animals = load_colony_rows(
            args["fieldnames"],
            transformed_rows,
            column_mapping=args["column_mapping"],
            gene_columns=args["gene_columns"],
            source_name=source_name,
            first_data_row=args["first_data_row"],
            allow_unknown_sex=args["allow_unknown_sex"],
        )
    except ValueError as exc:
        return UploadLoadResult(
            animals=[],
            warnings=grouped_warnings,
            unknown_sex_count=0,
            unknown_sex_ids=[],
            distinct_values=distinct_values,
            error_message=str(exc),
        )

    pedigree = pedigree_warnings(animals)
    pedigree_warnings_only = [
        f"warning: {message}"
        for severity, message in pedigree
        if severity == "warning"
    ]
    info_count = sum(severity == "info" for severity, _ in pedigree)
    displayed_pedigree = pedigree_warnings_only[:DISTINCT_VALUE_LIMIT]
    undisplayed_count = len(pedigree_warnings_only) - len(displayed_pedigree)
    if undisplayed_count:
        displayed_pedigree.append(f"+{undisplayed_count} more")
    if info_count:
        displayed_pedigree.append(
            f"{info_count} informational pedigree message(s)."
        )
    if displayed_pedigree:
        grouped_warnings["Pedigree"] = displayed_pedigree

    missing_parent_columns = [
        key
        for key in ("mother_id", "father_id")
        if column_mapping.get(key) in (None, NONE_COLUMN)
    ]
    if missing_parent_columns:
        grouped_warnings.setdefault("Relatedness", []).append(
            "Relatedness flags are unavailable because these parent columns are "
            "unmapped: "
            + ", ".join(missing_parent_columns)
            + "."
        )
    if column_mapping.get("cage") in (None, NONE_COLUMN):
        grouped_warnings.setdefault("Cage", []).append(
            "The same_cage flag is unavailable because the cage column is unmapped."
        )
    unknown_sex_ids = [animal.animal_id for animal in animals if not animal.sex]
    return UploadLoadResult(
        animals=animals,
        warnings=grouped_warnings,
        unknown_sex_count=len(unknown_sex_ids),
        unknown_sex_ids=unknown_sex_ids,
        distinct_values=distinct_values,
    )


def _validate_gene_selections(
    fieldnames: Sequence[str],
    gene_columns: Sequence[str],
    column_mapping: Mapping[str, str],
) -> None:
    available = set(fieldnames)
    seen = set()
    mapped = set(column_mapping.values())
    for gene in gene_columns:
        if gene not in available:
            raise ValueError(
                f"Selected gene column {gene!r} is no longer available. "
                "Choose a column from the current upload."
            )
        if gene in seen:
            raise ValueError(f"Gene column {gene!r} was selected more than once.")
        if gene in mapped:
            raise ValueError(
                f"Column {gene!r} is already mapped to a colony field and cannot also be a gene."
            )
        seen.add(gene)


def _apply_sex_values(
    rows: Sequence[Mapping[Optional[str], object]],
    sex_column: str,
    sex_value_mapping: Mapping[str, str],
) -> List[Dict[Optional[str], object]]:
    values = _distinct_raw_values(rows, sex_column)
    stale_values = set(sex_value_mapping).difference(values)
    if stale_values:
        raise ValueError(
            "Sex mapping contains values no longer present in the upload: "
            + ", ".join(repr(value) for value in sorted(stale_values))
        )
    unmapped_values = set(values).difference(sex_value_mapping)
    if unmapped_values:
        raise ValueError(
            "Choose M, F, or unknown for every sex-column value. Unmapped: "
            + ", ".join(repr(value) for value in sorted(unmapped_values))
        )

    transformed = []
    for row in rows:
        copy = dict(row)
        raw_value = copy.get(sex_column, "")
        value = "" if raw_value is None else str(raw_value)
        mapped = sex_value_mapping[value]
        normalized = mapped.strip().upper()
        if normalized not in {"M", "F", "", "UNKNOWN"}:
            raise ValueError(
                f"Invalid sex mapping {mapped!r} for value {value!r}; "
                "choose M, F, or unknown."
            )
        copy[sex_column] = "" if normalized in {"", "UNKNOWN"} else normalized
        transformed.append(copy)
    return transformed


def _group_upload_warnings(warnings: Sequence[str]) -> Dict[str, List[str]]:
    grouped: Dict[str, List[str]] = {}
    for warning in warnings:
        column = getattr(warning, "column", None)
        key = column if column is not None else "Upload"
        grouped.setdefault(key, []).append(str(warning))
    return grouped
