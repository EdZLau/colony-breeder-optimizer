"""Pure presentation logic for the Streamlit upload workflow."""

import csv
from copy import deepcopy
from dataclasses import dataclass
from fractions import Fraction
from io import StringIO
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

from colony_breeder.app.ui_helpers import (
    NONE_COLUMN,
    UploadLoadResult,
    build_loader_args,
    guess_column_mapping,
    guess_sex_value_mapping,
)
from colony_breeder.app.upload import UploadResult
from colony_breeder.app.run_settings import RunSettings, build_run_settings
from colony_breeder.app.run_ranking import (
    ExcludedAnimalView,
    RunResult,
    SkippedPairView,
    skipped_counts,
    visible_ranked,
)
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


def build_ranking_upload_args(
    upload: UploadResult,
    column_mapping: Mapping[str, Optional[str]],
    loaded_animals: Sequence[Animal],
) -> Dict[str, object]:
    """Build copied loader inputs with sex values normalized by the load step."""
    args = build_loader_args(
        upload.fieldnames,
        upload.rows,
        column_mapping,
        upload.first_data_row,
    )
    normalized_mapping = args["column_mapping"]
    animals_by_id = {animal.animal_id: animal for animal in loaded_animals}
    animal_id_column = normalized_mapping["animal_id"]
    sex_column = normalized_mapping["sex"]

    for row in args["rows"]:
        raw_animal_id = row.get(animal_id_column, "")
        animal_id = "" if raw_animal_id is None else str(raw_animal_id).strip()
        if not animal_id:
            continue
        animal = animals_by_id.get(animal_id)
        if animal is None:
            raise ValueError(
                f"Animal {animal_id!r} is missing from the loaded colony."
            )
        row[sex_column] = animal.sex
    return args


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


def stale_linked_pair_notices(
    previous: Optional[RunSettings],
    current: Optional[RunSettings],
) -> Tuple[str, ...]:
    """Describe linked-pair settings removed while cleaning stale selections."""
    if previous is None:
        return ()

    remaining = list(current.linked_pairs) if current is not None else []
    notices = []
    for index, pair in enumerate(previous.linked_pairs, start=1):
        if pair in remaining:
            remaining.remove(pair)
        else:
            notices.append(
                f"Linked pair {index} was removed because its gene or phase "
                "column is no longer available."
            )
    return tuple(notices)


def format_fraction_with_decimal(probability: Fraction) -> str:
    """Show an exact probability alongside a readable decimal approximation."""
    return f"{probability} ({float(probability):.4f})"


def split_error_message(message: str) -> Tuple[str, ...]:
    """Split semicolon-delimited validation/run errors into display bullets."""
    return tuple(part.strip() for part in message.split("; ") if part.strip())


def format_flags(flags: Sequence[str]) -> str:
    """Format core flag names without changing their wording."""
    return ", ".join(flags) if flags else "-"


def _descriptive_gene_names(settings: RunSettings) -> Set[str]:
    """Return genes whose per-gene probabilities are descriptive marginals."""
    genes = set(settings.x_linked_genes)
    genes.update(
        gene
        for gene1, gene2, _, _ in settings.linked_pairs
        for gene in (gene1, gene2)
    )
    return genes


def ranked_pair_rows(
    result: RunResult,
    settings: RunSettings,
) -> List[Dict[str, str]]:
    """Build display rows for visible pairs, preserving the core ranking order."""
    shown_pairs, _, _ = visible_ranked(result, settings.top)
    descriptive_genes = _descriptive_gene_names(settings)
    rows = []
    for pair in shown_pairs:
        per_gene = "; ".join(
            f"{gene}={probability}"
            + (" (descriptive)" if gene in descriptive_genes else "")
            for gene, probability in pair.per_gene
        )
        linked_factors = "; ".join(
            f"{gene1},{gene2}={format_fraction_with_decimal(probability)}"
            for (gene1, gene2), probability in pair.linked_pair_probabilities
        )
        rows.append(
            {
                "Male": pair.male_id,
                "Female": pair.female_id,
                "Probability (exact; decimal)": format_fraction_with_decimal(
                    pair.probability
                ),
                "Flags": format_flags(pair.flags),
                "Per-gene values (linked/X values marked descriptive; do not multiply those)": (
                    per_gene or "-"
                ),
                "Linked-pair factors": linked_factors or "-",
            }
        )
    return rows


_RANKED_PAIRS_CSV_HEADER = (
    "male",
    "female",
    "probability_numerator",
    "probability_denominator",
    "probability_decimal",
    "flags",
    "per_gene",
    "linked_pair_factors",
)


def _spreadsheet_safe_text(value: str) -> str:
    """Prefix formula-leading text so spreadsheet programs treat it as text."""
    if value.startswith(("=", "+", "-", "@", "\t", "\r")):
        return "'" + value
    return value


def ranked_pairs_csv(result: RunResult, settings: RunSettings) -> bytes:
    """Export every nonzero ranked pair, independent of the display top limit."""
    output = StringIO(newline="")
    writer = csv.writer(output, lineterminator="\r\n")
    writer.writerow(_RANKED_PAIRS_CSV_HEADER)

    if result.error_message is None:
        ranked_pairs, _, _ = visible_ranked(result, top=None)
        descriptive_genes = _descriptive_gene_names(settings)
        for pair in ranked_pairs:
            per_gene = "; ".join(
                f"{gene}={probability}"
                + (
                    " (descriptive) (do not multiply)"
                    if gene in descriptive_genes
                    else ""
                )
                for gene, probability in pair.per_gene
            )
            linked_factors = "; ".join(
                f"{gene1},{gene2}={format_fraction_with_decimal(probability)}"
                for (gene1, gene2), probability in pair.linked_pair_probabilities
            )
            flags = format_flags(pair.flags) if pair.flags else ""
            writer.writerow(
                (
                    _spreadsheet_safe_text(pair.male_id),
                    _spreadsheet_safe_text(pair.female_id),
                    pair.probability.numerator,
                    pair.probability.denominator,
                    float(pair.probability),
                    _spreadsheet_safe_text(flags) if flags else "",
                    _spreadsheet_safe_text(per_gene) if per_gene else "",
                    _spreadsheet_safe_text(linked_factors) if linked_factors else "",
                )
            )

    return output.getvalue().encode("utf-8-sig")


def run_result_summary(result: RunResult, top: Optional[int]) -> Tuple[str, ...]:
    """Format complete pair-count summaries without omitting hidden categories."""
    shown_pairs, zero_count, hidden_by_top = visible_ranked(result, top)
    counts = skipped_counts(result)
    return (
        f"Ranked pairs: {len(result.ranked)}",
        f"Visible nonzero pairs: {len(shown_pairs)}",
        f"Zero-probability pairs hidden: {zero_count}",
        f"Pairs hidden by top: {hidden_by_top}",
        f"Skipped pairs: {len(result.skipped)}",
        f"Excluded animals: {len(result.excluded_animals)}",
        "Skipped pairs by kind: "
        + ", ".join(f"{kind}={count}" for kind, count in counts.items()),
    )


def skipped_pair_rows(
    skipped: Sequence[SkippedPairView],
) -> List[Dict[str, str]]:
    """Return every skipped pair as a display-ready row."""
    return [
        {
            "Male": pair.male_id,
            "Female": pair.female_id,
            "Kind": pair.kind,
            "Reason": pair.reason,
            "Flags": format_flags(pair.flags),
        }
        for pair in skipped
    ]


def excluded_animal_rows(
    excluded_animals: Sequence[ExcludedAnimalView],
) -> List[Dict[str, str]]:
    """Return every excluded animal as a display-ready row."""
    return [
        {"Animal": animal.animal_id, "Reason": animal.reason}
        for animal in excluded_animals
    ]


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
