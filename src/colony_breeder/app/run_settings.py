"""Pure validation and loader-argument helpers for ranking settings."""

from dataclasses import dataclass
from fractions import Fraction
import re
from types import MappingProxyType
from typing import Iterable, Mapping, Optional, Sequence, Tuple, Union

from colony_breeder.app.ui_helpers import NONE_COLUMN
from colony_breeder.core.genetics import parse_genotype, parse_x_linked_genotype
from colony_breeder.core.linkage import linked_offspring_distribution

GenotypeValues = Union[str, Iterable[str]]
LinkedPairInput = Tuple[str, str, str, str]
_RECOMBINATION_TEXT = re.compile(r"(?:[0-9]+(?:\.[0-9]+)?|[0-9]+/[0-9]+)\Z")


@dataclass(frozen=True)
class RunSettings:
    """Normalized, immutable settings for loading and ranking a colony."""

    targets: Mapping[str, Tuple[str, ...]]
    lethals: Mapping[str, Tuple[str, ...]]
    x_linked_genes: Tuple[str, ...]
    linked_pairs: Tuple[LinkedPairInput, ...]
    exclude_related: bool
    top: Optional[int]


def _normalise_genotypes(
    values_by_gene: Optional[Mapping[str, GenotypeValues]],
) -> Mapping[str, Tuple[str, ...]]:
    normalized = {}
    for raw_gene, raw_values in (values_by_gene or {}).items():
        gene = raw_gene.strip() if isinstance(raw_gene, str) else str(raw_gene).strip()
        if isinstance(raw_values, str):
            values = (raw_values,)
        else:
            values = tuple(raw_values)
        seen = set()
        normalized_values = []
        for raw_value in values:
            if not isinstance(raw_value, str):
                raise ValueError(
                    f"Genotype values for gene {gene!r} must be text."
                )
            value = raw_value.strip()
            if not value:
                continue
            if value not in seen:
                seen.add(value)
                normalized_values.append(value)
        if gene not in normalized:
            normalized[gene] = normalized_values
        else:
            for value in normalized_values:
                if value not in normalized[gene]:
                    normalized[gene].append(value)
    return MappingProxyType(
        {gene: tuple(values) for gene, values in normalized.items()}
    )


def _normalise_names(names: Iterable[str]) -> Tuple[str, ...]:
    if isinstance(names, str):
        names = (names,)
    return tuple(
        dict.fromkeys(
            name.strip() if isinstance(name, str) else str(name).strip()
            for name in names
        )
    )


def build_run_settings(
    targets: Mapping[str, GenotypeValues],
    lethals: Optional[Mapping[str, GenotypeValues]] = None,
    x_linked_genes: Iterable[str] = (),
    linked_pairs: Iterable[LinkedPairInput] = (),
    exclude_related: bool = False,
    top: Optional[int] = None,
) -> RunSettings:
    """Strip text, remove duplicate values, and preserve input order."""
    normalized_pairs = []
    for pair in linked_pairs:
        try:
            gene1, gene2, r_text, phase_column = pair
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Each linked pair must contain gene1, gene2, r text, and phase column."
            ) from exc
        if not all(isinstance(value, str) for value in (gene1, gene2, r_text, phase_column)):
            raise ValueError("Linked-pair genes, r, and phase column must be text.")
        normalized_pairs.append(
            (
                gene1.strip(),
                gene2.strip(),
                r_text.strip(),
                phase_column.strip(),
            )
        )
    return RunSettings(
        targets=_normalise_genotypes(targets),
        lethals=_normalise_genotypes(lethals),
        x_linked_genes=_normalise_names(x_linked_genes),
        linked_pairs=tuple(normalized_pairs),
        exclude_related=bool(exclude_related),
        top=top,
    )


def _validate_genotypes(
    gene: str,
    genotypes: Sequence[str],
    x_linked: bool,
    category: str,
    errors: list,
) -> None:
    for genotype in genotypes:
        try:
            if x_linked:
                parts = genotype.strip().split("/")
                sex = (
                    "M"
                    if len(parts) == 2 and parts[1].strip().lower() == "y"
                    else "F"
                )
                parse_x_linked_genotype(genotype, sex)
            else:
                parse_genotype(genotype)
        except ValueError as exc:
            errors.append(
                f"Invalid {category} genotype {genotype!r} for gene {gene!r}: {exc}"
            )


def validate_run_settings(
    settings: RunSettings,
    gene_pool: Sequence[str],
    fieldnames: Sequence[Optional[str]],
    column_mapping: Mapping[str, Optional[str]],
) -> list:
    """Return user-readable errors; an empty list means settings are valid."""
    errors = []
    available_genes = set(gene_pool)
    available_columns = set(fieldnames)
    mapped_columns = {
        column
        for column in column_mapping.values()
        if column is not None and column != NONE_COLUMN
    }

    if not settings.targets:
        errors.append("Choose at least one target gene and genotype.")

    configured_genes = dict.fromkeys(
        [
            *settings.targets,
            *settings.lethals,
            *settings.x_linked_genes,
            *(
                gene
                for gene1, gene2, _, _ in settings.linked_pairs
                for gene in (gene1, gene2)
            ),
        ]
    )
    for gene in configured_genes:
        if gene in available_genes and gene not in available_columns:
            errors.append(
                f"Selected gene {gene!r} is not present in the uploaded headers."
            )
        if gene in mapped_columns:
            errors.append(
                f"Selected gene {gene!r} is mapped to a colony field."
            )

    for category, settings_by_gene in (
        ("target", settings.targets),
        ("lethal", settings.lethals),
    ):
        for gene, genotypes in settings_by_gene.items():
            if not gene:
                errors.append(f"A {category} setting has a blank gene name.")
            elif gene not in available_genes:
                errors.append(
                    f"{category.title()} gene {gene!r} is not in the selected gene columns."
                )
            if category == "target" and not genotypes:
                errors.append(f"Target genotype list for gene {gene!r} cannot be empty.")
            _validate_genotypes(
                gene,
                genotypes,
                gene in settings.x_linked_genes,
                category,
                errors,
            )

    for gene in settings.x_linked_genes:
        if not gene:
            errors.append("An X-linked gene name cannot be blank.")
        elif gene not in available_genes:
            errors.append(
                f"X-linked gene {gene!r} is not in the selected gene columns."
            )

    linked_genes = set()
    used_phase_columns = set()
    for linked_pair in settings.linked_pairs:
        if len(linked_pair) != 4:
            errors.append(
                "Each linked pair must contain gene1, gene2, r text, and phase column."
            )
            continue
        gene1, gene2, r_text, phase_column = linked_pair
        label = f"({gene1!r}, {gene2!r})"
        if not gene1 or not gene2:
            errors.append(f"Linked pair {label} must name two genes.")
        if gene1 == gene2:
            errors.append(f"Linked pair {label} cannot use the same gene twice.")
        for gene in dict.fromkeys((gene1, gene2)):
            if gene and gene not in available_genes:
                errors.append(
                    f"Linked gene {gene!r} is not in the selected gene columns."
                )
            if gene in linked_genes:
                errors.append(
                    f"Linked gene {gene!r} is repeated within or across linked pairs."
                )
            linked_genes.add(gene)
            if gene in settings.x_linked_genes:
                errors.append(
                    f"Linked pair {label} cannot include X-linked gene {gene!r}."
                )

        if not phase_column:
            errors.append(f"Linked pair {label} needs a phase column.")
        else:
            if phase_column not in available_columns:
                errors.append(
                    f"Phase column {phase_column!r} for linked pair {label} "
                    "is not present in the uploaded headers."
                )
            if phase_column in mapped_columns:
                errors.append(
                    f"Phase column {phase_column!r} cannot be mapped to a colony field."
                )
            if phase_column in available_genes:
                errors.append(
                    f"Phase column {phase_column!r} cannot also be a selected gene column."
                )
            if phase_column in used_phase_columns:
                errors.append(
                    f"Phase column {phase_column!r} cannot be shared between linked pairs."
                )
            used_phase_columns.add(phase_column)

        try:
            parse_recombination_fraction(r_text)
        except ValueError as exc:
            errors.append(f"Invalid recombination fraction for linked pair {label}: {exc}")

    if settings.top is not None and (
        isinstance(settings.top, bool)
        or not isinstance(settings.top, int)
        or settings.top < 1
    ):
        errors.append("Top N must be an integer greater than or equal to 1.")

    return errors


def build_reload_args(
    upload_args: Mapping[str, object],
    settings: RunSettings,
) -> dict:
    """Add configured genotype columns and linked options to upload arguments."""
    gene_columns = list(
        dict.fromkeys(
            [
                *settings.targets,
                *settings.lethals,
                *settings.x_linked_genes,
                *(
                    gene
                    for gene1, gene2, _, _ in settings.linked_pairs
                    for gene in (gene1, gene2)
                ),
            ]
        )
    )
    phase_columns = {phase for _, _, _, phase in settings.linked_pairs}
    gene_columns = [gene for gene in gene_columns if gene not in phase_columns]

    reload_args = dict(upload_args)
    reload_args.update(
        gene_columns=gene_columns,
        x_linked_genes=list(settings.x_linked_genes),
        linked_pairs=[
            (gene1, gene2, parse_recombination_fraction(r_text), phase_column)
            for gene1, gene2, r_text, phase_column in settings.linked_pairs
        ],
        allow_unknown_sex=True,
    )
    return reload_args


def parse_recombination_fraction(r_text: str) -> Fraction:
    """Validate recombination text with the core and return its exact fraction."""
    if not isinstance(r_text, str):
        raise ValueError("r must be entered as text, not a float or boolean.")
    cleaned = r_text.strip()
    if not _RECOMBINATION_TEXT.fullmatch(cleaned):
        raise ValueError(
            f"Invalid recombination fraction {r_text!r}: enter digits, a decimal, "
            "or a fraction such as '1/10'."
        )
    try:
        linked_offspring_distribution("A|B/A|B", "A|B/A|B", cleaned)
        return Fraction(cleaned)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"Invalid recombination fraction {r_text!r}: {exc}") from exc


def _canonical_genotype(text: str, x_linked: bool) -> Tuple[str, str]:
    if not x_linked:
        return parse_genotype(text)
    parts = text.strip().split("/")
    sex = (
        "M"
        if len(parts) == 2 and parts[1].strip().lower() == "y"
        else "F"
    )
    return parse_x_linked_genotype(text, sex)


def _canonical_genotype_settings(
    values_by_gene: Mapping[str, Tuple[str, ...]],
    x_linked_genes: set,
) -> Tuple[Tuple[str, Tuple[Tuple[str, str], ...]], ...]:
    return tuple(
        (
            gene,
            tuple(
                sorted(
                    {
                        _canonical_genotype(genotype, gene in x_linked_genes)
                        for genotype in genotypes
                    }
                )
            ),
        )
        for gene, genotypes in sorted(values_by_gene.items())
    )


def run_signature(base_signature: tuple, settings: RunSettings) -> tuple:
    """Build a hashable signature for all settings that affect ranking output."""
    x_linked_genes = set(settings.x_linked_genes)
    linked_pairs = tuple(
        sorted(
            (
                gene1,
                gene2,
                parse_recombination_fraction(r_text),
                phase_column,
            )
            for gene1, gene2, r_text, phase_column in settings.linked_pairs
        )
    )
    return (
        base_signature,
        _canonical_genotype_settings(settings.targets, x_linked_genes),
        _canonical_genotype_settings(settings.lethals, x_linked_genes),
        tuple(sorted(x_linked_genes)),
        linked_pairs,
        settings.exclude_related,
    )
