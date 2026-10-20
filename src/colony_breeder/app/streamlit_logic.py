"""Pure presentation logic for the Streamlit upload workflow."""

from copy import deepcopy
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

from colony_breeder.app.ui_helpers import (
    NONE_COLUMN,
    UploadLoadResult,
    guess_column_mapping,
    guess_sex_value_mapping,
)
from colony_breeder.app.upload import UploadResult
from colony_breeder.app.run_settings import RunSettings, build_run_settings
from colony_breeder.core.colony import Animal

OPTIONAL_COLUMN_KEYS = ("cage", "mother_id", "father_id")
SEX_MAPPING_DISPLAY_LIMIT = 20


@dataclass(frozen=True)
class SexMappingCheck:
    """Summary of whether distinct sex values can reasonably be mapped."""

    distinct_value_count: int
    displayable: bool
    warning: Optional[str]


@dataclass(frozen=True)
class LoadResultView:
    """Immutable, display-ready summary of an upload load result."""

    animal_count: int
    unknown_sex_count: int
    unknown_sex_ids: Tuple[str, ...]
    grouped_warnings: Tuple[Tuple[str, Tuple[str, ...]], ...]
    error_message: Optional[str]
    gene_value_summaries: Tuple[Tuple[str, Tuple[str, ...]], ...]


def default_column_selections(
    fieldnames: Sequence[str],
    previous: Optional[Mapping[str, Optional[str]]] = None,
) -> Dict[str, Optional[str]]:
    """Keep previous column choices that still exist; otherwise use name guesses."""
    guesses = guess_column_mapping(fieldnames)
    available = set(fieldnames)
    if previous is None:
        return guesses

    selections = {}
    for key, guess in guesses.items():
        prior = previous.get(key)
        if key in OPTIONAL_COLUMN_KEYS and prior == NONE_COLUMN:
            selections[key] = NONE_COLUMN
        elif prior in available:
            selections[key] = prior
        else:
            selections[key] = guess
    return selections


def default_sex_value_mapping(
    rows: Sequence[Mapping[Optional[str], object]],
    sex_column: str,
    previous: Optional[Mapping[str, str]] = None,
) -> Dict[str, str]:
    """Preserve mappings for current raw values and guess mappings for new values."""
    guesses = guess_sex_value_mapping(rows, sex_column)
    if previous is None:
        return guesses
    return {
        value: previous.get(value, guess)
        for value, guess in guesses.items()
    }


def sex_value_label(raw_value: str) -> str:
    """Return a display-only label for a raw sex value."""
    return "(blank)" if raw_value == "" else raw_value


def sex_mapping_check(
    rows: Sequence[Mapping[Optional[str], object]],
    sex_column: str,
    limit: int = SEX_MAPPING_DISPLAY_LIMIT,
) -> SexMappingCheck:
    """Count raw distinct sex values and say whether to display the mapping table."""
    distinct_value_count = len(guess_sex_value_mapping(rows, sex_column))
    displayable = distinct_value_count <= limit
    warning = None
    if not displayable:
        warning = (
            f"This column has {distinct_value_count} distinct values, more than "
            f"the limit of {limit}; it probably isn't a sex column."
        )
    return SexMappingCheck(
        distinct_value_count=distinct_value_count,
        displayable=displayable,
        warning=warning,
    )


def gene_column_options(
    fieldnames: Sequence[str],
    column_mapping: Mapping[str, Optional[str]],
    excluded_columns: Sequence[str] = (),
) -> List[str]:
    """Return unmapped header columns, excluding explicitly reserved columns."""
    mapped_columns = {
        column
        for column in column_mapping.values()
        if column is not None and column != NONE_COLUMN
    }
    excluded = mapped_columns.union(excluded_columns)
    return [column for column in fieldnames if column not in excluded]


def phase_column_options(
    fieldnames: Sequence[str],
    column_mapping: Mapping[str, Optional[str]],
    gene_pool: Sequence[str],
) -> List[str]:
    """Return unmapped header columns that are not selected as genes."""
    return gene_column_options(
        fieldnames,
        column_mapping,
        excluded_columns=gene_pool,
    )


def drop_stale_run_settings(
    settings: RunSettings,
    gene_pool: Sequence[str],
    fieldnames: Sequence[str],
    column_mapping: Mapping[str, Optional[str]],
) -> RunSettings:
    """Remove settings whose genes or phase columns are no longer available."""
    mapped_columns = {
        column
        for column in column_mapping.values()
        if column is not None and column != NONE_COLUMN
    }
    available_genes = set(gene_pool).intersection(fieldnames) - mapped_columns
    available_phases = set(phase_column_options(fieldnames, column_mapping, gene_pool))
    targets = {
        gene: genotypes
        for gene, genotypes in settings.targets.items()
        if gene in available_genes
    }
    lethals = {
        gene: genotypes
        for gene, genotypes in settings.lethals.items()
        if gene in available_genes
    }
    x_linked_genes = [
        gene for gene in settings.x_linked_genes if gene in available_genes
    ]
    linked_pairs = [
        pair
        for pair in settings.linked_pairs
        if (not pair[0] or pair[0] in available_genes)
        and (not pair[1] or pair[1] in available_genes)
        and (
            not pair[3]
            or (
                pair[3] in fieldnames
                and pair[3] not in mapped_columns
                and pair[3] in available_phases
            )
        )
    ]
    return build_run_settings(
        targets=targets,
        lethals=lethals,
        x_linked_genes=x_linked_genes,
        linked_pairs=linked_pairs,
        exclude_related=settings.exclude_related,
        top=settings.top,
    )


def build_settings_from_widgets(
    target_genotypes: Mapping[str, str],
    lethal_genotypes: Mapping[str, str],
    x_linked_genes: Sequence[str],
    linked_pairs: Sequence[Tuple[str, str, str, str]],
    exclude_related: bool,
    top: Optional[int],
) -> RunSettings:
    """Convert widget text areas and linked-pair rows into RunSettings."""
    targets = {
        gene: text.splitlines()
        for gene, text in target_genotypes.items()
    }
    lethals = {
        gene: text.splitlines()
        for gene, text in lethal_genotypes.items()
    }
    return build_run_settings(
        targets=targets,
        lethals=lethals,
        x_linked_genes=x_linked_genes,
        linked_pairs=linked_pairs,
        exclude_related=exclude_related,
        top=top,
    )


def drop_stale_genes(
    selected: Sequence[str],
    options: Sequence[str],
) -> List[str]:
    """Keep selected genes that remain available, preserving their order."""
    available = set(options)
    return [gene for gene in selected if gene in available]


def input_signature(
    file_name: str,
    content_hash: str,
    sheet_name: Optional[str],
    header_row: int,
    column_mapping: Mapping[str, Optional[str]],
    sex_mapping: Mapping[str, str],
    gene_selection: Sequence[str],
) -> Tuple[object, ...]:
    """Return a stable, hashable key for the complete set of load inputs."""
    return (
        file_name,
        content_hash,
        sheet_name,
        header_row,
        tuple(sorted(column_mapping.items())),
        tuple(sorted(sex_mapping.items())),
        tuple(gene_selection),
    )


def preview_rows(
    upload: UploadResult,
    limit: int = 10,
) -> List[Dict[Optional[str], object]]:
    """Return a bounded deep copy of upload rows for safe presentation."""
    return deepcopy(upload.rows[:limit])


def format_load_result(result: UploadLoadResult) -> LoadResultView:
    """Convert a helper result into an immutable view without rewriting warnings."""
    grouped_warnings = tuple(
        (group, tuple(messages))
        for group, messages in result.warnings.items()
    )
    gene_value_summaries = tuple(
        (gene, tuple(values))
        for gene, values in result.distinct_values.items()
    )
    return LoadResultView(
        animal_count=len(result.animals),
        unknown_sex_count=result.unknown_sex_count,
        unknown_sex_ids=tuple(result.unknown_sex_ids),
        grouped_warnings=grouped_warnings,
        error_message=result.error_message,
        gene_value_summaries=gene_value_summaries,
    )


def colony_sex_warnings(animals: Sequence[Animal]) -> List[str]:
    """Warn when loaded animals do not include at least one male and one female."""
    has_males = any(animal.sex == "M" for animal in animals)
    has_females = any(animal.sex == "F" for animal in animals)
    warnings = []
    if not has_males:
        warnings.append("No males are present among the loaded animals.")
    if not has_females:
        warnings.append("No females are present among the loaded animals.")
    return warnings
