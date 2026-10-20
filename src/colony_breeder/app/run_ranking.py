"""Pure orchestration and immutable views for a colony ranking run."""

from dataclasses import dataclass
from copy import deepcopy
from fractions import Fraction
from typing import Mapping, Optional, Sequence, Tuple

from colony_breeder.app.run_settings import (
    RunSettings,
    build_reload_args,
    run_signature,
    validate_run_settings,
)
from colony_breeder.core.colony import load_colony_rows, pedigree_warnings
from colony_breeder.core.ranking import PairRankingResult, rank_pairs

_CORE_SKIPPED_KINDS = (
    "missing_genotype",
    "missing_phase",
    "excluded_related",
    "all_lethal",
)


@dataclass(frozen=True)
class RankedPairView:
    male_id: str
    female_id: str
    probability: Fraction
    per_gene: Tuple[Tuple[str, Fraction], ...]
    linked_pair_probabilities: Tuple[
        Tuple[Tuple[str, str], Fraction], ...
    ]
    flags: Tuple[str, ...]


@dataclass(frozen=True)
class SkippedPairView:
    male_id: str
    female_id: str
    kind: str
    reason: str
    flags: Tuple[str, ...]


@dataclass(frozen=True)
class ExcludedAnimalView:
    animal_id: str
    reason: str


@dataclass(frozen=True)
class RunResult:
    ranked: Tuple[RankedPairView, ...]
    skipped: Tuple[SkippedPairView, ...]
    excluded_animals: Tuple[ExcludedAnimalView, ...]
    notes: Tuple[str, ...]
    pedigree_warnings: Tuple[Tuple[str, str], ...]
    signature: Optional[tuple]
    error_message: Optional[str]


def _error_result(message: str, signature: Optional[tuple] = None) -> RunResult:
    return RunResult((), (), (), (), (), signature, message)


def _ranked_view(result: PairRankingResult) -> Tuple[RankedPairView, ...]:
    return tuple(
        RankedPairView(
            male_id=pair.male.animal_id,
            female_id=pair.female.animal_id,
            probability=pair.probability,
            per_gene=tuple(pair.per_gene.items()),
            linked_pair_probabilities=tuple(
                pair.linked_pair_probabilities.items()
            ),
            flags=tuple(pair.flags),
        )
        for pair in result.ranked
    )


def _skipped_view(result: PairRankingResult) -> Tuple[SkippedPairView, ...]:
    return tuple(
        SkippedPairView(
            male_id=pair.male.animal_id,
            female_id=pair.female.animal_id,
            kind=pair.kind,
            reason=pair.reason,
            flags=tuple(pair.flags),
        )
        for pair in result.skipped
    )


def run_ranking(
    upload_args: Mapping[str, object],
    settings: RunSettings,
    gene_pool: Sequence[str],
    fieldnames: Sequence[Optional[str]],
    column_mapping: Mapping[str, Optional[str]],
    base_signature: tuple,
) -> RunResult:
    """Validate settings, reload genotypes, rank, and return immutable views."""
    errors = validate_run_settings(
        settings,
        gene_pool,
        fieldnames,
        column_mapping,
    )
    if errors:
        return _error_result("; ".join(errors))

    signature = None
    try:
        signature = run_signature(base_signature, settings)
        reload_args = build_reload_args(upload_args, settings)
        reload_args["rows"] = deepcopy(reload_args["rows"])
        animals = load_colony_rows(**reload_args)
        ranking = rank_pairs(
            animals,
            target=settings.targets,
            lethal=settings.lethals,
            exclude_related=settings.exclude_related,
            x_linked_genes=settings.x_linked_genes,
            linked_pairs=[
                (gene1, gene2, recombination)
                for gene1, gene2, recombination, _ in reload_args["linked_pairs"]
            ],
        )
        warnings = tuple(pedigree_warnings(animals))
    except ValueError as exc:
        return _error_result(str(exc), signature)

    return RunResult(
        ranked=_ranked_view(ranking),
        skipped=_skipped_view(ranking),
        excluded_animals=tuple(
            ExcludedAnimalView(animal.animal_id, animal.reason)
            for animal in ranking.excluded_animals
        ),
        notes=tuple(ranking.notes),
        pedigree_warnings=warnings,
        signature=signature,
        error_message=None,
    )


def visible_ranked(
    result: RunResult,
    top: Optional[int],
) -> Tuple[Tuple[RankedPairView, ...], int, int]:
    """Return visible nonzero pairs and the CLI-style hidden-pair counts."""
    nonzero_pairs = tuple(
        pair for pair in result.ranked if pair.probability != Fraction(0)
    )
    zero_probability_count = len(result.ranked) - len(nonzero_pairs)
    shown_pairs = nonzero_pairs if top is None else nonzero_pairs[:top]
    hidden_by_top_count = len(nonzero_pairs) - len(shown_pairs)
    return shown_pairs, zero_probability_count, hidden_by_top_count


def skipped_counts(result: RunResult) -> dict:
    """Count skipped pairs using the core ranking categories."""
    counts = {kind: 0 for kind in _CORE_SKIPPED_KINDS}
    for pair in result.skipped:
        counts[pair.kind] = counts.get(pair.kind, 0) + 1
    return counts


def result_is_current(result: RunResult, current_signature: tuple) -> bool:
    """Return whether a successful result matches the current settings."""
    return result.error_message is None and result.signature == current_signature
