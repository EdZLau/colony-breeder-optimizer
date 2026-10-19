"""Genetics for one autosomal or X-linked gene.

Everything here is pure logic: no files, no printing, no command line.
Probabilities are fractions.Fraction, so 1/3 stays exactly 1/3.
X-linked crosses assume a 1:1 sex ratio at conception, no sex-specific
viability apart from listed lethals, no pseudoautosomal region or X-inactivation,
and no recombination.

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
XLinkedState = Tuple[str, Genotype]
XLinkedDistribution = Dict[XLinkedState, Fraction]


class AllOffspringLethalError(ValueError):
    """Raised when a cross has no surviving offspring."""


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


def format_x_linked_genotype(genotype: Genotype) -> str:
    """Display an X-linked genotype with the Y chromosome last, if present."""
    if "Y" in genotype:
        allele = next(allele for allele in genotype if allele != "Y")
        return f"{allele}/Y"
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
        raise AllOffspringLethalError(
            "Every possible offspring genotype is lethal; nothing survives."
        )

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


def parse_x_linked_genotype(text: str, sex: str) -> Genotype:
    """Parse slash-style X-linked genotypes for females ("A/a") or males ("A/Y").

    "y" is accepted as the male sex-chromosome marker and normalized to "Y".
    The returned tuple is sorted like other internal genotypes.
    """
    cleaned = text.strip()
    normalized_sex = sex.strip().upper()
    if normalized_sex not in {"F", "M"}:
        raise ValueError(f"Invalid sex {sex!r} for X-linked genotype; expected 'F' or 'M'.")
    if "/" not in cleaned:
        expected = "'A/a'" if normalized_sex == "F" else "'A/Y'"
        raise ValueError(
            f"Cannot read X-linked genotype {text!r}: expected slash style like {expected}."
        )

    parts = [part.strip() for part in cleaned.split("/")]
    if len(parts) != 2 or not all(parts):
        raise ValueError(
            f"Cannot read X-linked genotype {text!r}: expected slash style with two parts, "
            "like 'A/a' or 'A/Y'."
        )

    first, second = parts
    if normalized_sex == "M":
        if first.lower() == "y" or second.lower() != "y":
            raise ValueError(
                f"Invalid X-linked genotype {text!r} for male; expected one allele "
                "followed by '/Y', like 'A/Y'."
            )
        return tuple(sorted((first, "Y")))

    if any(part.lower() == "y" for part in parts):
        raise ValueError(
            f"Invalid X-linked genotype {text!r} for female; expected two X alleles "
            "in slash style, like 'A/a'."
        )
    return tuple(sorted(parts))


def _parse_x_linked_genotype_without_sex(text: str) -> Genotype:
    """Parse a target or lethal genotype whose form identifies its offspring sex."""
    cleaned = text.strip()
    if "/" not in cleaned:
        raise ValueError(
            f"Cannot read X-linked genotype {text!r}: expected slash style like 'A/a' or 'A/Y'."
        )
    parts = [part.strip() for part in cleaned.split("/")]
    if len(parts) != 2 or not all(parts):
        raise ValueError(
            f"Cannot read X-linked genotype {text!r}: expected slash style with two parts, "
            "like 'A/a' or 'A/Y'."
        )
    sex = "M" if parts[1].lower() == "y" else "F"
    return parse_x_linked_genotype(cleaned, sex)


def x_linked_offspring_distribution(
    mother: str,
    father: str,
    lethal: Iterable[str] = (),
) -> XLinkedDistribution:
    """Return a joint {(sex, genotype): probability} distribution for X-linked pups.

    Sex is "F" for daughters and "M" for sons. The 1:1 sex ratio applies at
    conception; lethal outcomes are removed from the joint distribution and
    surviving sex/genotype outcomes are then renormalized together.
    """
    maternal_genotype = parse_x_linked_genotype(mother, "F")
    paternal_genotype = parse_x_linked_genotype(father, "M")
    lethal_set = {
        _parse_x_linked_genotype_without_sex(text)
        for text in _as_list(lethal)
    }
    maternal_gametes = gametes(maternal_genotype)
    paternal_x = next(allele for allele in paternal_genotype if allele != "Y")

    distribution: XLinkedDistribution = {}
    for maternal_allele, maternal_probability in maternal_gametes.items():
        for paternal_allele, sex in ((paternal_x, "F"), ("Y", "M")):
            genotype = tuple(sorted((maternal_allele, paternal_allele)))
            state = (sex, genotype)
            distribution[state] = (
                distribution.get(state, Fraction(0))
                + maternal_probability * Fraction(1, 2)
            )

    survivors = {
        state: probability
        for state, probability in distribution.items()
        if state[1] not in lethal_set
    }
    total = sum(survivors.values(), Fraction(0))
    if total == 0:
        raise AllOffspringLethalError(
            "Every possible offspring genotype is lethal; nothing survives."
        )
    return {state: probability / total for state, probability in survivors.items()}


def x_linked_probability_of_target(
    mother: str,
    father: str,
    target: Union[str, Iterable[str]],
    lethal: Iterable[str] = (),
) -> Fraction:
    """Return the all-pup probability of any target X-linked genotype.

    Male targets use "A/Y"; female targets use "A/a". The probability
    includes the sex ratio and is conditional on survival after joint lethal
    renormalization.
    """
    targets = {
        _parse_x_linked_genotype_without_sex(text)
        for text in _as_list(target)
    }
    distribution = x_linked_offspring_distribution(mother, father, lethal)
    return sum(
        (
            probability
            for (_, genotype), probability in distribution.items()
            if genotype in targets
        ),
        Fraction(0),
    )
