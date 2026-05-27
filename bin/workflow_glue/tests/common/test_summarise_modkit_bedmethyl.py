"""Tests for modkit bedMethyl summarisation."""

import gzip

from workflow_glue import summarise_modkit_bedmethyl


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def test_summarise_bedmethyl_aggregates_per_mod_code(tmp_path):
    """Rows should be summed into one row per modification code."""
    label_tsv = _write(
        tmp_path / "mod_code_labels.tsv",
        "mod_code\tlabel\nA:a\tm6A\nT:17802\tpseU\n",
    )
    bedmethyl = _write(
        tmp_path / "sample.mods.bedmethyl",
        (
            "chr1\t10\t11\ta\t1\t+\t10\t11\t255,0,0\t5\t40.00\t2\t3\t0\t0\t1\t0\t0\n"
            "chr1\t11\t12\ta\t1\t+\t11\t12\t255,0,0\t7\t57.14\t4\t3\t0\t0\t2\t0\t0\n"
            "chr1\t20\t21\t17802\t1\t+\t20\t21\t255,0,0\t4\t25.00\t1\t3\t0\t0\t0\t0\t"
            "0\n"
        ),
    )

    rows = summarise_modkit_bedmethyl.summarise_bedmethyl(
        bedmethyl,
        "sample1",
        "A:a,T:17802",
        summarise_modkit_bedmethyl.load_mod_code_labels(label_tsv),
    )

    assert rows == [
        {
            "sample": "sample1",
            "full_mod_code": "T:17802",
            "mod_code": "17802",
            "mod_label": "pseU",
            "valid_coverage": 4,
            "modified_calls": 1,
            "canonical_calls": 3,
            "other_calls": 0,
            "delete_calls": 0,
            "fail_calls": 0,
            "diff_calls": 0,
            "nocall_calls": 0,
            "modification_percent": "25.00",
        },
        {
            "sample": "sample1",
            "full_mod_code": "A:a",
            "mod_code": "a",
            "mod_label": "m6A",
            "valid_coverage": 12,
            "modified_calls": 6,
            "canonical_calls": 6,
            "other_calls": 0,
            "delete_calls": 0,
            "fail_calls": 3,
            "diff_calls": 0,
            "nocall_calls": 0,
            "modification_percent": "50.00",
        },
    ]


def test_summarise_bedmethyl_reads_gzipped_input(tmp_path):
    """The helper should aggregate multiple rows from gzipped bedMethyl input."""
    bedmethyl = tmp_path / "sample.mods.bedmethyl.gz"
    with gzip.open(bedmethyl, "wt", encoding="utf-8") as handle:
        handle.write(
            "chr1\t10\t11\ta\t1\t+\t10\t11\t255,0,0\t5\t40.00\t2\t3\t0\t0\t1\t0\t0\n"
            "chr1\t11\t12\ta\t1\t+\t11\t12\t255,0,0\t7\t57.14\t4\t3\t0\t0\t2\t0\t0\n"
        )

    rows = summarise_modkit_bedmethyl.summarise_bedmethyl(
        bedmethyl,
        "sample1",
        "A:a",
        {"A:a": "m6A"},
    )

    assert rows == [
        {
            "sample": "sample1",
            "full_mod_code": "A:a",
            "mod_code": "a",
            "mod_label": "m6A",
            "valid_coverage": 12,
            "modified_calls": 6,
            "canonical_calls": 6,
            "other_calls": 0,
            "delete_calls": 0,
            "fail_calls": 3,
            "diff_calls": 0,
            "nocall_calls": 0,
            "modification_percent": "50.00",
        }
    ]
