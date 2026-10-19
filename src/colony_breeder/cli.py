"""Command-line interface for ranking colony breeding pairs."""

import argparse
import csv
import sys
from typing import List, Optional, Sequence, Tuple

from colony_breeder.core.colony import load_colony, pedigree_warnings
from colony_breeder.core.ranking import PairRankingResult, rank_pairs


def _gene_genotype(value: str) -> Tuple[str, str]:
    """Parse one GENE=GENOTYPE argument, preserving later equals signs."""
    if "=" not in value:
        raise argparse.ArgumentTypeError("expected GENE=GENOTYPE")
    gene, genotype = value.split("=", 1)
    gene = gene.strip()
    genotype = genotype.strip()
    if not gene:
        raise argparse.ArgumentTypeError("gene name cannot be empty")
    if not genotype:
        raise argparse.ArgumentTypeError("genotype cannot be empty")
    return gene, genotype


def _positive_integer(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer greater than or equal to 1") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be an integer greater than or equal to 1")
    return parsed


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Rank male x female colony pairs by target genotype probability."
    )
    parser.add_argument("csv_path", help="path to a colony CSV file")
    parser.add_argument(
        "--target",
        action="append",
        type=_gene_genotype,
        required=True,
        metavar="GENE=GENOTYPE",
        help="target genotype; repeat for multiple targets or genes",
    )
    parser.add_argument(
        "--lethal",
        action="append",
        type=_gene_genotype,
        default=[],
        metavar="GENE=GENOTYPE",
        help="lethal genotype; repeat as needed",
    )
    parser.add_argument(
        "--top",
        type=_positive_integer,
        metavar="N",
        help="show at most the top N nonzero-probability pairs",
    )
    parser.add_argument(
        "--exclude-related",
        action="store_true",
        help="skip pairs with sibling or parent-offspring relationship flags",
    )
    parser.add_argument(
        "--show-skipped",
        action="store_true",
        help="show every skipped pair reason instead of the first five",
    )
    return parser


def _settings(entries: Sequence[Tuple[str, str]]) -> dict:
    settings = {}
    for gene, genotype in entries:
        settings.setdefault(gene, []).append(genotype)
    return settings


def _ascii(text: str) -> str:
    return text.encode("ascii", "backslashreplace").decode("ascii")


def _print_results(result: PairRankingResult, top: Optional[int], show_skipped: bool) -> None:
    zero_probability_count = sum(pair.probability == 0 for pair in result.ranked)
    nonzero_pairs = [pair for pair in result.ranked if pair.probability != 0]
    visible_pairs = nonzero_pairs
    if top is not None:
        visible_pairs = visible_pairs[:top]

    print("Assumptions: autosomal genes; target genes unlinked; all animals fertile and available.")
    if not nonzero_pairs:
        print("No pair can produce the target.")
    print("Ranked pairs:")
    print(
        f"{'Rank':>4}  {'Male':<12}  {'Female':<12}  "
        f"{'Probability (exact)':<20}  {'Decimal':>7}  Flags"
    )
    for rank, pair in enumerate(visible_pairs, start=1):
        flags = ",".join(pair.flags) if pair.flags else "-"
        print(_ascii(
            f"{rank:>4}  {pair.male.animal_id:<12}  {pair.female.animal_id:<12}  "
            f"{str(pair.probability):<20}  {float(pair.probability):>7.4f}  {flags}"
        ))
    print(f"Showing {len(visible_pairs)} of {len(nonzero_pairs)} pairs with nonzero probability")
    print(f"Hidden zero-probability pairs: {zero_probability_count}")

    categories = {
        "missing_genotype": 0,
        "excluded_related": 0,
        "all_lethal": 0,
    }
    for pair in result.skipped:
        if pair.kind in categories:
            categories[pair.kind] += 1
    print(f"Skipped pairs: {len(result.skipped)}")
    print(
        "Skipped categories: "
        f"missing genotype={categories['missing_genotype']}, "
        f"excluded related={categories['excluded_related']}, "
        f"all offspring lethal={categories['all_lethal']}"
    )
    skipped_to_show = result.skipped if show_skipped else result.skipped[:5]
    if skipped_to_show:
        print("Skipped pair examples:")
        for pair in skipped_to_show:
            print(_ascii(
                f"{pair.male.animal_id} x {pair.female.animal_id}: "
                f"{pair.reason}"
            ))


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Parse arguments, call core functions, and print the CLI result."""
    args = _parser().parse_args(argv)
    targets = _settings(args.target)
    lethals = _settings(args.lethal)
    genes = list(dict.fromkeys([*targets, *lethals]))
    ignored_lethal_genes = [gene for gene in lethals if gene not in targets]

    try:
        colony = load_colony(args.csv_path, gene_columns=genes)
        warnings = pedigree_warnings(colony)
        result = rank_pairs(
            colony,
            target=targets,
            lethal=lethals,
            exclude_related=args.exclude_related,
        )
    except UnicodeDecodeError:
        print(
            "Error: CSV file is not valid UTF-8; save it as a UTF-8 CSV file.",
            file=sys.stderr,
        )
        return 1
    except (ValueError, OSError, csv.Error) as exc:
        print(_ascii(f"Error: {exc}"), file=sys.stderr)
        return 1

    for severity, message in warnings:
        if severity == "warning":
            print(_ascii(f"Pedigree warning: {message}"))
    info_count = sum(severity == "info" for severity, _ in warnings)
    print(f"Pedigree info messages: {info_count}")
    for gene in ignored_lethal_genes:
        print(_ascii(f"Note: lethal setting for '{gene}' ignored (not a target gene)"))
    _print_results(result, args.top, args.show_skipped)
    return 0


if __name__ == "__main__":
    sys.exit(main())
