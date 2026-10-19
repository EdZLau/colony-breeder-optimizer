"""Generate and rank male x female breeding pairs."""

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from colony_breeder.core.colony import Animal
from colony_breeder.core.genetics import (
    AllOffspringLethalError,
    parse_genotype,
    probability_of_target,
)


@dataclass(frozen=True)
class RankedPair:
    """A male-female pair and its probability of producing the target."""

    male: Animal
    female: Animal
    probability: Fraction


@dataclass(frozen=True)
class SkippedPair:
    """A male-female pair omitted because required data is unknown or lethal."""

    male: Animal
    female: Animal
    reason: str


@dataclass(frozen=True)
class PairRankingResult:
    """Ranked eligible pairs alongside pairs that could not be evaluated."""

    ranked: List[RankedPair]
    skipped: List[SkippedPair]


def generate_pairs(animals: Sequence[Animal]) -> List[Tuple[Animal, Animal]]:
    """Return every male x female pair from the colony."""
    males = [animal for animal in animals if animal.sex == "M"]
    females = [animal for animal in animals if animal.sex == "F"]
    return [(male, female) for male in males for female in females]


def rank_pairs(
    animals: Sequence[Animal],
    target: Mapping[str, Union[str, Iterable[str]]],
    lethal: Optional[Mapping[str, Union[str, Iterable[str]]]] = None,
) -> PairRankingResult:
    """Rank pairs for one target gene, preserving unknown data as skipped pairs.

    ``target`` and ``lethal`` are keyed by gene name to support the planned
    multi-gene interface, while this step evaluates exactly one target gene.
    """
    colony_genes = {gene for animal in animals for gene in animal.genotypes}
    target_genes = set(target)
    lethal_by_gene = lethal or {}
    lethal_genes = set(lethal_by_gene)
    unknown_genes = (target_genes | lethal_genes) - colony_genes
    if unknown_genes:
        raise ValueError(
            "Target or lethal settings name gene(s) not present in the colony: "
            + ", ".join(sorted(unknown_genes))
        )

    if len(target) != 1:
        raise ValueError("target must contain exactly one gene for single-gene ranking.")

    gene, target_genotypes = next(iter(target.items()))
    if not gene.strip():
        raise ValueError("target gene name cannot be blank.")

    target_genotypes = _validate_genotypes(target_genotypes)
    unsupported_lethal_genes = lethal_genes - {gene}
    if unsupported_lethal_genes:
        raise ValueError(
            "Single-gene ranking only supports lethal genotypes for the target gene "
            f"{gene!r}; got {sorted(unsupported_lethal_genes)!r}."
        )

    lethal_genotypes = lethal_by_gene.get(gene, ())
    lethal_genotypes = _validate_genotypes(lethal_genotypes)
    ranked: List[RankedPair] = []
    skipped: List[SkippedPair] = []

    for male, female in generate_pairs(animals):
        missing = []
        if not male.genotypes.get(gene, "").strip():
            missing.append(f"male {male.animal_id}")
        if not female.genotypes.get(gene, "").strip():
            missing.append(f"female {female.animal_id}")
        if missing:
            skipped.append(
                SkippedPair(
                    male=male,
                    female=female,
                    reason=(
                        f"Missing genotype for gene {gene!r} in "
                        + " and ".join(missing)
                        + "."
                    ),
                )
            )
            continue

        try:
            probability = probability_of_target(
                male.genotypes[gene],
                female.genotypes[gene],
                target_genotypes,
                lethal=lethal_genotypes,
            )
        except AllOffspringLethalError:
            skipped.append(
                SkippedPair(
                    male=male,
                    female=female,
                    reason=(
                        f"All offspring are lethal for gene {gene!r} from "
                        f"male {male.animal_id} and female {female.animal_id}."
                    ),
                )
            )
            continue
        ranked.append(RankedPair(male=male, female=female, probability=probability))

    ranked.sort(
        key=lambda pair: (
            -pair.probability,
            pair.male.animal_id,
            pair.female.animal_id,
        )
    )
    return PairRankingResult(ranked=ranked, skipped=skipped)


def _validate_genotypes(
    genotypes: Union[str, Iterable[str]],
) -> Union[str, List[str]]:
    """Raise before ranking if a target or lethal genotype is malformed."""
    if isinstance(genotypes, str):
        values = [genotypes]
        normalized: Union[str, List[str]] = genotypes
    else:
        values = list(genotypes)
        normalized = values
    for genotype in values:
        parse_genotype(genotype)
    return normalized
