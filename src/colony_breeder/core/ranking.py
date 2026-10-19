"""Generate and rank male x female breeding pairs."""

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from colony_breeder.core.colony import Animal
from colony_breeder.core.genetics import (
    AllOffspringLethalError,
    parse_genotype,
    probability_of_target,
)
from colony_breeder.core.relatedness import pair_flags


@dataclass(frozen=True)
class RankedPair:
    """A male-female pair and its probability of producing the target."""

    male: Animal
    female: Animal
    probability: Fraction
    per_gene: Dict[str, Fraction]
    flags: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class SkippedPair:
    """A male-female pair omitted because required data is unknown or lethal."""

    male: Animal
    female: Animal
    reason: str
    flags: List[str] = field(default_factory=list)


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
    exclude_related: bool = False,
) -> PairRankingResult:
    """Rank pairs by the probability of producing all targeted genotypes.

    Target genes are assumed to be unlinked, so their exact per-gene
    probabilities are multiplied. A lethal list applies only to its gene.
    """
    if not target:
        raise ValueError("target must contain at least one gene.")

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

    targets = {}
    for gene, genotypes in target.items():
        if not gene.strip():
            raise ValueError("target gene name cannot be blank.")
        validated = _validate_genotypes(genotypes)
        if not validated:
            raise ValueError(f"target genotype list for gene {gene!r} cannot be empty.")
        targets[gene] = validated

    lethals = {
        gene: _validate_genotypes(genotypes)
        for gene, genotypes in lethal_by_gene.items()
    }
    ranked: List[RankedPair] = []
    skipped: List[SkippedPair] = []

    for male, female in generate_pairs(animals):
        flags = pair_flags(male, female)
        missing_by_gene = {}
        for gene in targets:
            missing = []
            if not male.genotypes.get(gene, "").strip():
                missing.append(f"male {male.animal_id}")
            if not female.genotypes.get(gene, "").strip():
                missing.append(f"female {female.animal_id}")
            if missing:
                missing_by_gene[gene] = missing

        if missing_by_gene:
            missing_details = [
                f"gene {gene!r} in " + " and ".join(missing)
                for gene, missing in missing_by_gene.items()
            ]
            skipped.append(
                SkippedPair(
                    male=male,
                    female=female,
                    reason=(
                        "Missing genotype for "
                        + "; ".join(missing_details)
                        + "."
                    ),
                    flags=flags,
                )
            )
            continue

        excluded_flags = [
            flag
            for flag in flags
            if flag
            in {
                "full_siblings",
                "half_siblings_mother",
                "half_siblings_father",
                "parent_offspring",
            }
        ]
        if exclude_related and excluded_flags:
            labels = {
                "full_siblings": "full siblings",
                "half_siblings_mother": "maternal half siblings",
                "half_siblings_father": "paternal half siblings",
                "parent_offspring": "parent-offspring",
            }
            skipped.append(
                SkippedPair(
                    male=male,
                    female=female,
                    reason="Excluded: " + ", ".join(labels[flag] for flag in excluded_flags),
                    flags=flags,
                )
            )
            continue

        per_gene = {}
        for gene, target_genotypes in targets.items():
            try:
                per_gene[gene] = probability_of_target(
                    male.genotypes[gene],
                    female.genotypes[gene],
                    target_genotypes,
                    lethal=lethals.get(gene, []),
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
                        flags=flags,
                    )
                )
                break
        else:
            probability = Fraction(1)
            for gene_probability in per_gene.values():
                probability *= gene_probability
            ranked.append(
                RankedPair(
                    male=male,
                    female=female,
                    probability=probability,
                    per_gene=per_gene,
                    flags=flags,
                )
            )

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
) -> List[str]:
    """Raise before ranking if a target or lethal genotype is malformed."""
    values = [genotypes] if isinstance(genotypes, str) else list(genotypes)
    for genotype in values:
        parse_genotype(genotype)
    return values
