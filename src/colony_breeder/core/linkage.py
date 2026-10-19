"""Exact inheritance calculations for one linked pair of autosomal genes.

A linked pair must use the joint calculation whenever either gene is targeted
or has a lethal list; an untargeted partner's lethal genotypes can change the
surviving target probability. The pair can be ignored only when neither gene
is targeted and neither has lethals.

This module assumes equal recombination in males and females, no crossover
interference, and a recombination frequency supplied by the user. Ten map
units are represented as r = 1/10. X-linked genes are not handled here and
remain unlinked from one another (r = 1/2) in this step.
"""
from fractions import Fraction
from typing import Dict, Iterable, List, Optional, Tuple, Union

from colony_breeder.core.genetics import (
    AllOffspringLethalError,
    Genotype,
    parse_genotype,
)

Haplotype = Tuple[str, str]
PhasedGenotype = Tuple[Haplotype, Haplotype]
LinkedGenotype = Tuple[Genotype, Genotype]
LinkedDistribution = Dict[LinkedGenotype, Fraction]
GenotypeInput = Union[str, Iterable[str]]
RecombinationInput = Union[Fraction, str, int]


def _as_list(value: GenotypeInput) -> List[str]:
    return [value] if isinstance(value, str) else list(value)


def _fraction(value: RecombinationInput) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError(
            "recombination frequency must be a Fraction, string, or integer, not float or bool."
        )
    if not isinstance(value, (Fraction, str, int)):
        raise ValueError("recombination frequency must be a Fraction, string, or integer.")
    try:
        frequency = value if isinstance(value, Fraction) else Fraction(str(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid recombination frequency {value!r}.") from exc
    if not Fraction(0) <= frequency <= Fraction(1, 2):
        raise ValueError("recombination frequency must be between 0 and 1/2 inclusive.")
    return frequency


def parse_phased_genotype(text: str) -> PhasedGenotype:
    """Parse e.g. ``A|B/a|b`` into a canonical, sorted haplotype pair."""
    if not isinstance(text, str):
        raise ValueError("phased genotype must be text in the form 'A|B/a|b'.")
    chromosomes = text.strip().split("/")
    if len(chromosomes) != 2:
        raise ValueError(
            f"Cannot read phased genotype {text!r}: expected two chromosomes "
            "separated by '/', each with two alleles separated by '|', like 'A|B/a|b'."
        )
    haplotypes = []
    for chromosome in chromosomes:
        alleles = [allele.strip() for allele in chromosome.split("|")]
        if len(alleles) != 2 or not all(alleles):
            raise ValueError(
                f"Cannot read phased genotype {text!r}: each chromosome must have "
                "two non-empty alleles separated by '|', like 'A|B/a|b'."
            )
        haplotypes.append((alleles[0], alleles[1]))
    return tuple(sorted(haplotypes))


def validate_phase_matches_genotypes(
    phased_genotype: Union[str, PhasedGenotype],
    gene1_genotype: str,
    gene2_genotype: str,
) -> None:
    """Raise ValueError unless the phase alleles match both per-gene genotypes."""
    phase = (
        parse_phased_genotype(phased_genotype)
        if isinstance(phased_genotype, str)
        else tuple(sorted(phased_genotype))
    )
    observed_gene1 = tuple(sorted(allele for haplotype in phase for allele in [haplotype[0]]))
    observed_gene2 = tuple(sorted(allele for haplotype in phase for allele in [haplotype[1]]))
    expected_gene1 = parse_genotype(gene1_genotype)
    expected_gene2 = parse_genotype(gene2_genotype)
    if observed_gene1 != expected_gene1 or observed_gene2 != expected_gene2:
        raise ValueError(
            "phased genotype does not agree with the per-gene genotypes: "
            f"phase gives {observed_gene1!r} and {observed_gene2!r}, expected "
            f"{expected_gene1!r} and {expected_gene2!r}."
        )


def phase_is_required(gene1_genotype: str, gene2_genotype: str) -> bool:
    """Return whether both genes are heterozygous and their phase matters."""
    return (
        len(set(parse_genotype(gene1_genotype))) == 2
        and len(set(parse_genotype(gene2_genotype))) == 2
    )


def phase_for_unambiguous_genotypes(
    gene1_genotype: str,
    gene2_genotype: str,
) -> PhasedGenotype:
    """Build canonical phase when at least one linked gene is homozygous."""
    gene1 = parse_genotype(gene1_genotype)
    gene2 = parse_genotype(gene2_genotype)
    if phase_is_required(gene1_genotype, gene2_genotype):
        raise ValueError("phase must be supplied when both linked genes are heterozygous.")
    return tuple(sorted(((gene1[0], gene2[0]), (gene1[1], gene2[1]))))


def _gametes(phase: PhasedGenotype, recombination_frequency: Fraction) -> Dict[Haplotype, Fraction]:
    first, second = phase
    parental_probability = (Fraction(1) - recombination_frequency) / 2
    recombinant_probability = recombination_frequency / 2
    result: Dict[Haplotype, Fraction] = {}
    for haplotype, probability in (
        (first, parental_probability),
        (second, parental_probability),
        ((first[0], second[1]), recombinant_probability),
        ((second[0], first[1]), recombinant_probability),
    ):
        if probability:
            result[haplotype] = result.get(haplotype, Fraction(0)) + probability
    return result


def _lethal_set(genotypes: Optional[GenotypeInput]) -> set:
    if genotypes is None:
        return set()
    return {parse_genotype(genotype) for genotype in _as_list(genotypes)}


def linked_offspring_distribution(
    parent1_phase: Union[str, PhasedGenotype],
    parent2_phase: Union[str, PhasedGenotype],
    recombination_frequency: RecombinationInput,
    lethal_gene1: Optional[GenotypeInput] = None,
    lethal_gene2: Optional[GenotypeInput] = None,
) -> LinkedDistribution:
    """Return the jointly lethal-filtered distribution for two linked genes.

    Each result key is ``(gene1_genotype, gene2_genotype)``. Lethals at either
    gene remove every joint outcome carrying that genotype; survivors are
    renormalized once across the complete two-gene distribution.
    """
    phase1 = (
        parse_phased_genotype(parent1_phase)
        if isinstance(parent1_phase, str)
        else tuple(sorted(parent1_phase))
    )
    phase2 = (
        parse_phased_genotype(parent2_phase)
        if isinstance(parent2_phase, str)
        else tuple(sorted(parent2_phase))
    )
    frequency = _fraction(recombination_frequency)
    lethal1 = _lethal_set(lethal_gene1)
    lethal2 = _lethal_set(lethal_gene2)
    gametes1 = _gametes(phase1, frequency)
    gametes2 = _gametes(phase2, frequency)

    distribution: LinkedDistribution = {}
    for haplotype1, probability1 in gametes1.items():
        for haplotype2, probability2 in gametes2.items():
            gene1 = tuple(sorted((haplotype1[0], haplotype2[0])))
            gene2 = tuple(sorted((haplotype1[1], haplotype2[1])))
            outcome = (gene1, gene2)
            distribution[outcome] = (
                distribution.get(outcome, Fraction(0)) + probability1 * probability2
            )

    survivors = {
        outcome: probability
        for outcome, probability in distribution.items()
        if outcome[0] not in lethal1 and outcome[1] not in lethal2
    }
    total = sum(survivors.values(), Fraction(0))
    if total == 0:
        raise AllOffspringLethalError(
            "Every possible linked offspring genotype is lethal; nothing survives."
        )
    return {outcome: probability / total for outcome, probability in survivors.items()}


def probability_of_linked_target(
    parent1_phase: Union[str, PhasedGenotype],
    parent2_phase: Union[str, PhasedGenotype],
    recombination_frequency: RecombinationInput,
    target_gene1: Optional[GenotypeInput] = None,
    target_gene2: Optional[GenotypeInput] = None,
    lethal_gene1: Optional[GenotypeInput] = None,
    lethal_gene2: Optional[GenotypeInput] = None,
) -> Fraction:
    """Return the probability of the specified per-gene targets among survivors.

    Either target may be ``None`` for an untargeted partner gene. Its lethal
    list, if supplied, still participates in joint survival conditioning.
    """
    if target_gene1 is None and target_gene2 is None:
        raise ValueError("at least one target genotype list is required.")
    targets1 = (
        {parse_genotype(genotype) for genotype in _as_list(target_gene1)}
        if target_gene1 is not None
        else None
    )
    targets2 = (
        {parse_genotype(genotype) for genotype in _as_list(target_gene2)}
        if target_gene2 is not None
        else None
    )
    if targets1 is not None and not targets1:
        raise ValueError("target genotype list for gene 1 cannot be empty.")
    if targets2 is not None and not targets2:
        raise ValueError("target genotype list for gene 2 cannot be empty.")

    distribution = linked_offspring_distribution(
        parent1_phase,
        parent2_phase,
        recombination_frequency,
        lethal_gene1=lethal_gene1,
        lethal_gene2=lethal_gene2,
    )
    return sum(
        (
            probability
            for (genotype1, genotype2), probability in distribution.items()
            if (targets1 is None or genotype1 in targets1)
            and (targets2 is None or genotype2 in targets2)
        ),
        Fraction(0),
    )
