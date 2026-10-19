"""Genetics for ONE autosomal gene.

Everything here is pure logic: no files, no printing, no command line.
Probabilities are fractions.Fraction, so 1/3 stays exactly 1/3.

Vocabulary used in this file
- allele:   one version of a gene, e.g. "A", "a", "Cre", "+".
- genotype: the two alleles an animal carries, stored as a sorted tuple of
            two strings, e.g. ("A", "a"). Sorting is what makes "aA" and
            "Aa" the same genotype.
- gamete:   an egg or sperm. It carries ONE allele from its parent.
"""
from fractions import Fraction
from typing import Dict, Iterable, Tuple, Union

Genotype = Tuple[str, str]
Distribution = Dict[Genotype, Fraction]


def parse_genotype(text: str) -> Genotype:
    """Turn text like "Aa", "aA", "Cre/+" or "+/Cre" into a sorted allele pair.

    Two accepted styles:
    - slash style: alleles separated by "/", each can be any length ("Cre/+").
    - compact style: exactly two characters ("Aa"), one per allele.
    """
    cleaned = text.strip()
    if "/" in cleaned:
        parts = [part.strip() for part in cleaned.split("/")]
    else:
        parts = list(cleaned)  # "Aa" -> ["A", "a"]

    if len(parts) != 2 or not all(parts):
        raise ValueError(
            f"Cannot read genotype {text!r}: expected two alleles, like 'Aa' "
            "or 'Cre/+' (use a slash when an allele name is longer than one character)."
        )

    first, second = sorted(parts)
    return (first, second)


def format_genotype(genotype: Genotype) -> str:
    """Display a genotype in slash style, e.g. ("A", "a") -> "A/a"."""
    return "/".join(genotype)


def gametes(genotype: Genotype) -> Dict[str, Fraction]:
    """Return {allele: probability} for the gametes this genotype produces.

    Mendel's first law: each of the two alleles goes into a gamete with
    probability 1/2. A homozygote ("AA") has the same allele twice, so the
    two halves add up to a single allele with probability 1.
    """
    result: Dict[str, Fraction] = {}
    for allele in genotype:
        result[allele] = result.get(allele, Fraction(0)) + Fraction(1, 2)
    return result


def _as_list(value: Union[str, Iterable[str]]) -> list:
    """Let callers pass one genotype string OR a collection of them."""
    if isinstance(value, str):
        return [value]
    return list(value)


def remove_lethal(distribution: Distribution, lethal: Iterable[str] = ()) -> Distribution:
    """Drop lethal genotypes, then rescale so the survivors sum to exactly 1.

    Example: Ay/a x Ay/a gives 1/4 Ay/Ay, 1/2 Ay/a, 1/4 a/a. If Ay/Ay is
    lethal, the survivors hold 3/4 of the probability, so we divide each by
    3/4 and get Ay/a = 2/3 and a/a = 1/3.
    """
    lethal_set = {parse_genotype(text) for text in _as_list(lethal)}
    survivors = {g: p for g, p in distribution.items() if g not in lethal_set}

    total = sum(survivors.values(), Fraction(0))
    if total == 0:
        raise ValueError("Every possible offspring genotype is lethal; nothing survives.")

    return {g: p / total for g, p in survivors.items()}


def offspring_distribution(
    parent1: str, parent2: str, lethal: Iterable[str] = ()
) -> Distribution:
    """Return {genotype: probability} for pups from parent1 x parent2.

    Every gamete from parent1 can meet every gamete from parent2 (random
    fertilisation), and the chance of a meeting is the product of the two
    gamete probabilities.
    """
    gametes1 = gametes(parse_genotype(parent1))
    gametes2 = gametes(parse_genotype(parent2))

    distribution: Distribution = {}
    for allele1, prob1 in gametes1.items():
        for allele2, prob2 in gametes2.items():
            first, second = sorted((allele1, allele2))
            child = (first, second)
            distribution[child] = distribution.get(child, Fraction(0)) + prob1 * prob2

    return remove_lethal(distribution, lethal)


def probability_of_target(
    parent1: str,
    parent2: str,
    target: Union[str, Iterable[str]],
    lethal: Iterable[str] = (),
) -> Fraction:
    """Probability that a (surviving) pup has the target genotype.

    `target` can be one genotype ("aa") or several that all count as a hit
    (["AA", "Aa"] for "shows the dominant trait").
    """
    distribution = offspring_distribution(parent1, parent2, lethal)
    targets = {parse_genotype(text) for text in _as_list(target)}
    return sum((p for g, p in distribution.items() if g in targets), Fraction(0))
