"""Generate and rank male x female breeding pairs."""

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from colony_breeder.core.colony import Animal
from colony_breeder.core.genetics import (
    AllOffspringLethalError,
    parse_genotype,
    parse_x_linked_genotype,
    probability_of_target,
    x_linked_offspring_distribution,
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
    kind: str
    flags: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class PairRankingResult:
    """Ranked eligible pairs alongside pairs that could not be evaluated."""

    ranked: List[RankedPair]
    skipped: List[SkippedPair]
    excluded_animals: List["ExcludedAnimal"] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


@dataclass(frozen=True)
class ExcludedAnimal:
    """An animal excluded from pairing because its sex is unusable."""

    animal_id: str
    reason: str


def generate_pairs(animals: Sequence[Animal]) -> List[Tuple[Animal, Animal]]:
    """Return every male x female pair from the colony."""
    males = [animal for animal in animals if animal.sex.strip().upper() == "M"]
    females = [animal for animal in animals if animal.sex.strip().upper() == "F"]
    return [(male, female) for male in males for female in females]


def rank_pairs(
    animals: Sequence[Animal],
    target: Mapping[str, Union[str, Iterable[str]]],
    lethal: Optional[Mapping[str, Union[str, Iterable[str]]]] = None,
    exclude_related: bool = False,
    x_linked_genes: Optional[Iterable[str]] = None,
) -> PairRankingResult:
    """Rank pairs by the probability of producing all targeted genotypes.

    Autosomal target probabilities are multiplied as before. X-linked genes
    are combined jointly by offspring sex and conditioned on the surviving
    joint distribution. Their per_gene entries are descriptive marginal
    probabilities among survivors, not factors to multiply. Multiple
    X-linked genes are treated as unlinked (r = 1/2), an approximation until
    step 6 accounts for linkage and recombination.

    Animals whose normalized sex is not "M" or "F" are excluded from pairing
    and reported in excluded_animals, regardless of inheritance settings.
    Missing targeted genotype data takes precedence over relatedness exclusion
    for pairs that can be formed.
    A blank genotype in an untargeted X-linked gene with a lethal list is also
    skipped as missing_genotype because that gene affects survival.
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

    x_linked_gene_order = list(dict.fromkeys(x_linked_genes or ()))
    x_linked_gene_names = set(x_linked_gene_order)
    unknown_x_linked_genes = x_linked_gene_names - colony_genes
    if unknown_x_linked_genes:
        raise ValueError(
            "X-linked gene name(s) not present in the colony: "
            + ", ".join(sorted(unknown_x_linked_genes))
        )
    for animal in animals:
        sex = animal.sex.strip().upper()
        if sex not in {"M", "F"}:
            continue
        for gene in x_linked_gene_order:
            genotype = animal.genotypes.get(gene, "").strip()
            if not genotype:
                continue
            try:
                parse_x_linked_genotype(genotype, sex)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid X-linked genotype for animal {animal.animal_id!r} "
                    f"in gene {gene!r}: {genotype!r}. {exc}"
                ) from exc

    targets = {}
    for gene, genotypes in target.items():
        if not gene.strip():
            raise ValueError("target gene name cannot be blank.")
        validated = (
            _validate_x_linked_genotypes(genotypes)
            if gene in x_linked_gene_names
            else _validate_genotypes(genotypes)
        )
        if not validated:
            raise ValueError(f"target genotype list for gene {gene!r} cannot be empty.")
        targets[gene] = validated

    lethals = {
        gene: (
            _validate_x_linked_genotypes(genotypes)
            if gene in x_linked_gene_names
            else _validate_genotypes(genotypes)
        )
        for gene, genotypes in lethal_by_gene.items()
    }
    x_linked_active_genes = list(dict.fromkeys(
        [
            gene
            for gene in targets
            if gene in x_linked_gene_names
        ]
        + [
            gene
            for gene in lethals
            if gene in x_linked_gene_names and lethals[gene]
        ]
    ))
    excluded_animals = [
        ExcludedAnimal(
            animal_id=animal.animal_id,
            reason=(
                "sex is blank or unrecognized, so it can't be assigned "
                "a male or female role"
            ),
        )
        for animal in animals
        if animal.sex.strip().upper() not in {"M", "F"}
    ]
    notes = (
        [
            "Some animals were excluded because sex is blank or unrecognized, "
            "so they can't be assigned a male or female role."
        ]
        if excluded_animals
        else []
    )
    ranked: List[RankedPair] = []
    skipped: List[SkippedPair] = []

    for male, female in generate_pairs(animals):
        flags = pair_flags(male, female)
        missing_by_gene = {}
        required_genotypes = list(dict.fromkeys([*targets, *x_linked_active_genes]))
        for gene in required_genotypes:
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
                    kind="missing_genotype",
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
                    kind="excluded_related",
                    flags=flags,
                )
            )
            continue

        per_gene = {}
        autosomal_probability = Fraction(1)
        gene_masses = {}
        if x_linked_active_genes:
            for gene in x_linked_active_genes:
                raw_distribution = x_linked_offspring_distribution(
                    female.genotypes[gene],
                    male.genotypes[gene],
                )
                target_set = (
                    {_x_linked_genotype_key(genotype) for genotype in targets[gene]}
                    if gene in targets
                    else None
                )
                lethal_set = {
                    _x_linked_genotype_key(genotype)
                    for genotype in lethals.get(gene, [])
                }
                gene_masses[gene] = {}
                for sex in ("F", "M"):
                    conditional_distribution = {
                        genotype: probability * 2
                        for (offspring_sex, genotype), probability in raw_distribution.items()
                        if offspring_sex == sex
                    }
                    surviving = {
                        genotype: probability
                        for genotype, probability in conditional_distribution.items()
                        if genotype not in lethal_set
                    }
                    gene_masses[gene][sex] = (
                        sum(surviving.values(), Fraction(0)),
                        sum(
                            (
                                probability
                                for genotype, probability in surviving.items()
                                if target_set is not None and genotype in target_set
                            ),
                            Fraction(0),
                        ),
                    )

            weight_by_sex = {}
            target_by_sex = {}
            for sex in ("F", "M"):
                weight = Fraction(1, 2)
                target_mass = Fraction(1, 2)
                for gene in x_linked_active_genes:
                    surviving_mass, gene_target_mass = gene_masses[gene][sex]
                    weight *= surviving_mass
                    target_mass *= (
                        gene_target_mass
                        if gene in targets
                        else surviving_mass
                    )
                weight_by_sex[sex] = weight
                target_by_sex[sex] = target_mass

            total_surviving_mass = sum(weight_by_sex.values(), Fraction(0))
            if total_surviving_mass == 0:
                lethal_genes_text = ", ".join(x_linked_active_genes)
                skipped.append(
                    SkippedPair(
                        male=male,
                        female=female,
                        reason=(
                            f"All offspring are lethal for X-linked gene(s) "
                            f"{lethal_genes_text!r} from male {male.animal_id} "
                            f"and female {female.animal_id}."
                        ),
                        kind="all_lethal",
                        flags=flags,
                    )
                )
                continue
            x_linked_probability = (
                sum(target_by_sex.values(), Fraction(0))
                / total_surviving_mass
            )
        else:
            total_surviving_mass = Fraction(1)
            x_linked_probability = Fraction(1)

        for gene, target_genotypes in targets.items():
            if gene in x_linked_gene_names:
                marginal_target_mass = Fraction(0)
                for sex in ("F", "M"):
                    branch_mass = Fraction(1, 2)
                    for other_gene in x_linked_active_genes:
                        surviving_mass, gene_target_mass = gene_masses[other_gene][sex]
                        branch_mass *= (
                            gene_target_mass
                            if other_gene == gene
                            else surviving_mass
                        )
                    marginal_target_mass += branch_mass
                per_gene[gene] = marginal_target_mass / total_surviving_mass
            else:
                try:
                    gene_probability = probability_of_target(
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
                            kind="all_lethal",
                            flags=flags,
                        )
                    )
                    break
                per_gene[gene] = gene_probability
                autosomal_probability *= gene_probability
        else:
            probability = autosomal_probability * x_linked_probability
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
    return PairRankingResult(
        ranked=ranked,
        skipped=skipped,
        excluded_animals=excluded_animals,
        notes=notes,
    )


def _validate_genotypes(
    genotypes: Union[str, Iterable[str]],
) -> List[str]:
    """Raise before ranking if a target or lethal genotype is malformed."""
    values = [genotypes] if isinstance(genotypes, str) else list(genotypes)
    for genotype in values:
        parse_genotype(genotype)
    return values


def _validate_x_linked_genotypes(
    genotypes: Union[str, Iterable[str]],
) -> List[str]:
    """Validate target/lethal forms whose notation identifies offspring sex."""
    values = [genotypes] if isinstance(genotypes, str) else list(genotypes)
    for genotype in values:
        parts = genotype.strip().split("/")
        sex = "M" if len(parts) == 2 and parts[1].strip().lower() == "y" else "F"
        parse_x_linked_genotype(genotype, sex)
    return values


def _x_linked_genotype_key(genotype: str) -> Tuple[str, str]:
    parts = genotype.strip().split("/")
    sex = "M" if len(parts) == 2 and parts[1].strip().lower() == "y" else "F"
    return parse_x_linked_genotype(genotype, sex)
