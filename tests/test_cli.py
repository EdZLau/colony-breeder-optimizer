from pathlib import Path

import pytest

from colony_breeder.cli import _gene_genotype, main

SAMPLE = Path(__file__).resolve().parent.parent / "data" / "sample_mice.csv"


def test_gene_genotype_splits_only_at_first_equals():
    assert _gene_genotype("Cre=Cre/=+") == ("Cre", "Cre/=+")


def test_main_ranks_pairs_with_repeated_targets_and_fraction_decimals(tmp_path, capsys):
    csv_path = tmp_path / "fraction.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        "F01,F,C2,,,+/+\n",
        encoding="utf-8",
    )

    exit_code = main(
        [
            str(csv_path),
            "--target",
            "Cre=Cre/+",
            "--target",
            "Cre=Cre/Cre",
            "--top",
            "1",
        ]
    )

    output = capsys.readouterr()

    assert exit_code == 0
    assert "1/2" in output.out
    assert "0.5000" in output.out
    assert "Hidden zero-probability pairs:" in output.out
    assert "Showing 1 of 1 pairs with nonzero probability" in output.out
    assert output.err == ""


def test_main_splits_target_once_and_loads_gene_also_named_lethal(capsys):
    exit_code = main(
        [
            str(SAMPLE),
            "--target",
            "Cre=Cre/=+",
            "--lethal",
            "Cre=Cre/+",
        ]
    )

    output = capsys.readouterr()
    assert exit_code == 0
    assert output.err == ""
    assert "Ranked pairs:" in output.out


@pytest.mark.parametrize(
    "option_value",
    [
        ("--target", "CreCre/+"),
        ("--target", "=Cre/+"),
        ("--target", "Cre="),
        ("--lethal", "CreCre/+"),
        ("--lethal", "=Cre/+"),
        ("--lethal", "Cre="),
    ],
)
def test_main_rejects_malformed_gene_arguments(option_value, capsys):
    option, value = option_value
    with pytest.raises(SystemExit) as exc_info:
        main([str(SAMPLE), "--target", "Cre=Cre/+", option, value])

    output = capsys.readouterr()
    assert exc_info.value.code == 2
    assert output.err.count("error:") == 1
    assert "Traceback" not in output.err


def test_main_requires_at_least_one_target(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main([str(SAMPLE)])

    output = capsys.readouterr()
    assert exc_info.value.code == 2
    assert "the following arguments are required: --target" in output.err


def test_main_rejects_top_less_than_one(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main([str(SAMPLE), "--target", "Cre=Cre/+", "--top", "0"])

    output = capsys.readouterr()
    assert exc_info.value.code == 2
    assert "greater than or equal to 1" in output.err


@pytest.mark.parametrize(
    "gene_argument",
    [
        ("--target", "Missing=Aa"),
        ("--lethal", "Missing=aa"),
    ],
)
def test_main_reports_unknown_genes_without_traceback(gene_argument, capsys):
    option, value = gene_argument
    assert main([str(SAMPLE), "--target", "Cre=Cre/+", option, value]) == 1

    output = capsys.readouterr()
    assert output.out == ""
    assert "Error:" in output.err
    assert "Traceback" not in output.err


def test_main_reports_missing_file_without_traceback(tmp_path, capsys):
    missing_file = tmp_path / "missing.csv"

    assert main([str(missing_file), "--target", "Cre=Cre/+"]) == 1

    output = capsys.readouterr()
    assert "Error:" in output.err
    assert "Traceback" not in output.err
    assert output.out == ""


def test_main_suggests_utf8_for_non_utf8_csv(tmp_path, capsys):
    csv_path = tmp_path / "not_utf8.csv"
    csv_path.write_bytes(b"animal_id,sex,cage,mother_id,father_id,Cre\nM01,M,C1,,,\xff\n")

    assert main([str(csv_path), "--target", "Cre=Cre/+"]) == 1

    output = capsys.readouterr()
    assert "save it as a UTF-8 CSV file" in output.err
    assert "Traceback" not in output.err
    assert output.out == ""


def test_main_reports_csv_parser_error_without_traceback(tmp_path, capsys):
    csv_path = tmp_path / "nul_byte.csv"
    csv_path.write_bytes(
        b"animal_id,sex,cage,mother_id,father_id,Cre\n"
        b"M01,M,C1,,,\x00\n"
    )

    assert main([str(csv_path), "--target", "Cre=Cre/+"]) == 1

    output = capsys.readouterr()
    assert output.out == ""
    assert output.err.startswith("Error:")
    assert output.err.count("\n") == 1
    assert "Traceback" not in output.err


def test_main_prints_pedigree_warnings_before_ranked_table(tmp_path, capsys):
    csv_path = tmp_path / "pedigree_warning.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,\n"
        "F01,F,C2,M01,,Cre/+\n",
        encoding="utf-8",
    )

    exit_code = main([str(csv_path), "--target", "Cre=Cre/+"])

    output = capsys.readouterr()
    assert exit_code == 0
    assert output.out.index("Pedigree warning:") < output.out.index("Ranked pairs:")
    assert output.out.index("Pedigree info messages:") < output.out.index("Ranked pairs:")


def test_main_prints_skipped_counts_and_limits_examples_by_default(tmp_path, capsys):
    csv_path = tmp_path / "skipped.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,\n"
        "F01,F,C2,,,\n"
        "F02,F,C3,,,\n"
        "F03,F,C4,,,\n"
        "F04,F,C5,,,\n"
        "F05,F,C6,,,\n",
        encoding="utf-8",
    )

    assert main([str(csv_path), "--target", "Cre=Cre/+"]) == 0

    output = capsys.readouterr()
    assert "missing genotype=5" in output.out
    assert output.out.count("Missing genotype for gene") == 5
    assert output.out.count("x F0") == 5


def test_show_skipped_prints_all_skips_and_empty_rankable_result(tmp_path, capsys):
    csv_path = tmp_path / "all_skipped.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/Cre\n"
        "F01,F,C2,,,Cre/Cre\n"
        "F02,F,C3,,,Cre/Cre\n"
        "F03,F,C4,,,Cre/Cre\n"
        "F04,F,C5,,,Cre/Cre\n"
        "F05,F,C6,,,Cre/Cre\n",
        encoding="utf-8",
    )

    assert main(
        [
            str(csv_path),
            "--target",
            "Cre=Cre/+",
            "--lethal",
            "Cre=Cre/Cre",
            "--show-skipped",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "Skipped pairs: 5" in output.out
    assert "all offspring lethal=5" in output.out
    assert output.out.count("All offspring are lethal") == 5
    assert "No pairs could be evaluated (see skipped pairs below)." in output.out
    assert "Showing 0 of 0 pairs with nonzero probability" in output.out


def test_exclude_related_and_show_skipped_report_related_count(capsys):
    assert main(
        [
            str(SAMPLE),
            "--target",
            "Cre=Cre/Cre",
            "--exclude-related",
            "--show-skipped",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "excluded related=12" in output.out
    assert "Excluded:" in output.out


def test_main_reports_errors_as_one_stderr_line(tmp_path, capsys):
    csv_path = tmp_path / "unknown_gene.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        "F01,F,C2,,,+/+\n",
        encoding="utf-8",
    )

    assert main([str(csv_path), "--target", "Missing=Aa"]) == 1

    output = capsys.readouterr()
    assert output.out == ""
    assert output.err.startswith("Error:")
    assert output.err.count("\n") == 1
    assert "Traceback" not in output.err


def test_main_reports_invalid_genotype_before_printing_results(capsys):
    assert main([str(SAMPLE), "--target", "Cre=Cre+"]) == 1

    output = capsys.readouterr()
    assert output.out == ""
    assert output.err.startswith("Error:")
    assert output.err.count("\n") == 1
    assert "Cannot read genotype" in output.err


def test_top_limits_number_of_printed_pair_rows(tmp_path, capsys):
    csv_path = tmp_path / "several_pairs.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Cre\n"
        "M01,M,C1,,,Cre/+\n"
        "F01,F,C2,,,+/+\n"
        "F02,F,C3,,,+/+\n"
        "F03,F,C4,,,+/+\n",
        encoding="utf-8",
    )

    assert main([str(csv_path), "--target", "Cre=Cre/+", "--top", "2"]) == 0

    output = capsys.readouterr()
    pair_rows = [
        line for line in output.out.splitlines()
        if line.startswith(("   1  ", "   2  ", "   3  "))
    ]
    assert len(pair_rows) == 2
    assert "Showing 2 of 3 pairs with nonzero probability" in output.out


def test_lethal_only_gene_prints_ignored_setting_note(capsys):
    assert main(
        [
            str(SAMPLE),
            "--target",
            "Cre=Cre/Cre",
            "--lethal",
            "Ob=ob/ob",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "Note: lethal setting for 'Ob' ignored (not a target gene)" in output.out


def test_untargeted_x_linked_lethal_affects_survival_without_ignored_note(
    tmp_path, capsys
):
    csv_path = tmp_path / "xlinked_lethal.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,X1,X2\n"
        "M01,M,C1,,,A/Y,A/Y\n"
        "F01,F,C2,,,A/a,A/a\n",
        encoding="utf-8",
    )

    assert main(
        [
            str(csv_path),
            "--target",
            "X2=a/Y",
            "--lethal",
            "X1=a/Y",
            "--xlinked",
            "X1",
            "--xlinked",
            "X2",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "1/6" in output.out
    assert "ignored (not a target gene)" not in output.out


def test_x_linked_and_autosomal_untargeted_lethals_have_distinct_notes(
    capsys,
):
    assert main(
        [
            str(SAMPLE),
            "--target",
            "Cre=Cre/Cre",
            "--lethal",
            "Ob=ob/ob",
        ]
    ) == 0
    autosomal_output = capsys.readouterr()
    assert "Note: lethal setting for 'Ob' ignored (not a target gene)" in autosomal_output.out


@pytest.mark.parametrize(
    ("animal_id", "sex", "genotype"),
    [
        ("M01", "M", "A/a"),
        ("F01", "F", "A/Y"),
    ],
)
def test_x_linked_sex_form_error_names_animal_and_gene_and_suggests_option(
    tmp_path, capsys, animal_id, sex, genotype
):
    csv_path = tmp_path / "invalid_xlinked_sex_form.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,Xgene\n"
        f"{animal_id},{sex},C1,,,{genotype}\n",
        encoding="utf-8",
    )

    assert main(
        [
            str(csv_path),
            "--target",
            "Xgene=A/a",
            "--xlinked",
            "Xgene",
        ]
    ) == 1

    output = capsys.readouterr()
    assert output.out == ""
    assert f"animal '{animal_id}'" in output.err
    assert "'Xgene'" in output.err
    assert output.err.rstrip().endswith(
        "If this gene is autosomal, remove it from --xlinked."
    )


def test_xlinked_cli_joint_probability_and_descriptive_footnote(tmp_path, capsys):
    csv_path = tmp_path / "two_xlinked_genes.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,X1,X2\n"
        "M01,M,C1,,,A/Y,A/Y\n"
        "F01,F,C2,,,A/a,A/a\n",
        encoding="utf-8",
    )

    assert main(
        [
            str(csv_path),
            "--target",
            "X2=a/Y",
            "--lethal",
            "X1=a/Y",
            "--xlinked",
            "X1",
            "--xlinked",
            "X2",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "1/6" in output.out
    assert "X-linked per_gene" in output.out
    assert "X2=1/6" in output.out
    assert (
        "Note: X-linked per-gene probabilities are descriptive marginal values "
        "among survivors and must not be multiplied. Multiple X-linked genes "
        "are treated as unlinked (r = 1/2), an approximation until step 6."
    ) in output.out


def test_unmarked_xlinked_gene_warning_does_not_change_xlinked_probability(
    tmp_path, capsys
):
    csv_path = tmp_path / "partially_marked_xlinked_genes.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,X1,X2\n"
        "M01,M,C1,,,A/Y,A/Y\n"
        "F01,F,C2,,,A/a,A/a\n",
        encoding="utf-8",
    )

    assert main(
        [
            str(csv_path),
            "--target",
            "X2=a/Y",
            "--lethal",
            "X1=a/Y",
            "--xlinked",
            "X2",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "Warning: gene 'X1' has genotypes like 'A/Y' but is not marked --xlinked." in output.out
    assert "If it is X-linked, add --xlinked X1; results for it are wrong otherwise." in output.out
    assert "1/2" in output.out


def test_marked_xlinked_genes_do_not_get_unmarked_warning(tmp_path, capsys):
    csv_path = tmp_path / "fully_marked_xlinked_genes.csv"
    csv_path.write_text(
        "animal_id,sex,cage,mother_id,father_id,X1,X2\n"
        "M01,M,C1,,,A/Y,A/Y\n"
        "F01,F,C2,,,A/a,A/a\n",
        encoding="utf-8",
    )

    assert main(
        [
            str(csv_path),
            "--target",
            "X2=a/Y",
            "--lethal",
            "X1=a/Y",
            "--xlinked",
            "X1",
            "--xlinked",
            "X2",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "Warning: gene " not in output.out


def test_autosomal_sample_run_has_no_unmarked_xlinked_warning(capsys):
    assert main([str(SAMPLE), "--target", "Cre=Cre/Cre"]) == 0

    output = capsys.readouterr()
    assert "Warning: gene " not in output.out


def test_sample_xgene_cli_command_ranks_expected_top_pair(capsys):
    assert main(
        [
            str(SAMPLE),
            "--xlinked",
            "Xgene",
            "--target",
            "Xgene=a/a",
            "--lethal",
            "Xgene=a/Y",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "M02" in output.out
    ranked_lines = [
        line for line in output.out.splitlines()
        if line.startswith("   1  ")
    ]
    assert len(ranked_lines) == 1
    assert ranked_lines[0].split() == [
        "1",
        "M02",
        "F03",
        "1",
        "1.0000",
        "-",
        "Xgene=1",
    ]
    assert "X-linked per-gene probabilities are descriptive marginal" in output.out


def test_sample_autosomal_only_run_has_no_xlinked_output_or_warning(capsys):
    assert main([str(SAMPLE), "--target", "Cre=Cre/Cre"]) == 0

    output = capsys.readouterr()
    assert "X-linked per_gene" not in output.out
    assert "X-linked per-gene probabilities" not in output.out
    assert "Warning: gene " not in output.out


def test_cli_prints_excluded_animals_with_five_example_limit(tmp_path, capsys):
    csv_path = tmp_path / "unknown_sexes.csv"
    rows = [
        "animal_id,sex,cage,mother_id,father_id,Cre,X\n",
        "M01,M,C1,,,Cre/+,A/Y\n",
        "F01,F,C2,,,+/+,A/a\n",
    ]
    rows.extend(f"U{i:02d},?,C{i + 3},,,Cre/+, \n" for i in range(1, 7))
    csv_path.write_text("".join(rows), encoding="utf-8")

    assert main(
        [
            str(csv_path),
            "--target",
            "Cre=Cre/+",
            "--xlinked",
            "X",
        ]
    ) == 0

    output = capsys.readouterr()
    assert "Some animals were excluded because sex is blank or unrecognized" in output.out
    assert "Excluded animals:" in output.out
    assert "Showing 1 of 1 pairs with nonzero probability" in output.out
    assert "0.5000" in output.out
    assert all(f"U{i:02d}:" in output.out for i in range(1, 6))
    assert "U06:" not in output.out
    assert "and 1 more" in output.out


def test_autosomal_only_output_has_no_x_linked_additions(capsys):
    assert main([str(SAMPLE), "--target", "Cre=Cre/Cre"]) == 0

    output = capsys.readouterr()
    assert "X-linked per_gene" not in output.out
    assert "unlinked (r = 1/2)" not in output.out
    assert "Excluded animals:" not in output.out
