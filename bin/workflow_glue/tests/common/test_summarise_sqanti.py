"""Tests for SQANTI summary extraction."""

import csv

import pytest
from workflow_glue import summarise_sqanti


def _args(*argv):
    return summarise_sqanti.argparser().parse_args(list(argv))


def test_main_writes_summary_from_structural_category(tmp_path):
    """The standard SQANTI structural_category column should be summarised."""
    sqanti_dir = tmp_path / "sqanti"
    sqanti_dir.mkdir()
    (sqanti_dir / "sample_classification.txt").write_text(
        (
            "isoform\tstructural_category\n"
            "tx1\tFSM\n"
            "tx2\tNIC\n"
            "tx3\tFSM\n"
        ),
        encoding="utf-8",
    )
    output = tmp_path / "classification_summary.tsv"

    summarise_sqanti.main(
        _args("--sqanti_dir", str(sqanti_dir), "--output", str(output))
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {"structural_category": "FSM", "count": "2"},
        {"structural_category": "NIC", "count": "1"},
    ]


def test_main_accepts_category_fallback_column(tmp_path):
    """Older SQANTI-style `category` columns should still be accepted."""
    sqanti_dir = tmp_path / "sqanti"
    sqanti_dir.mkdir()
    (sqanti_dir / "sample_classification.tsv").write_text(
        (
            "isoform\tcategory\n"
            "tx1\tnovel\n"
            "tx2\tknown\n"
            "tx3\tnovel\n"
        ),
        encoding="utf-8",
    )
    output = tmp_path / "classification_summary.tsv"

    summarise_sqanti.main(
        _args("--sqanti_dir", str(sqanti_dir), "--output", str(output))
    )

    with output.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows == [
        {"structural_category": "novel", "count": "2"},
        {"structural_category": "known", "count": "1"},
    ]


def test_main_errors_when_no_classification_table_exists(tmp_path):
    """A missing SQANTI classification table should fail clearly."""
    sqanti_dir = tmp_path / "sqanti"
    sqanti_dir.mkdir()
    output = tmp_path / "classification_summary.tsv"

    with pytest.raises(
        SystemExit, match="Could not find a SQANTI classification table"
    ):
        summarise_sqanti.main(
            _args("--sqanti_dir", str(sqanti_dir), "--output", str(output))
        )
