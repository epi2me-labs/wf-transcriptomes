"""Tests for DE/DTU design validation."""

import pytest
from workflow_glue import check_experiment_design


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def _args(*argv):
    return check_experiment_design.argparser().parse_args(list(argv))


def test_split_covariates_trims_and_drops_empty_values():
    """Covariates should be normalised into a clean list."""
    assert check_experiment_design._split_covariates(
        " batch, sex ,, site "
    ) == ["batch", "sex", "site"]


def test_main_accepts_valid_design(tmp_path):
    """A balanced two-condition sample sheet should validate cleanly."""
    sample_sheet = _write(
        tmp_path / "sample_sheet.csv",
        (
            "alias,condition,batch\n"
            "control_rep1,control,b1\n"
            "control_rep2,control,b2\n"
            "treated_rep1,treated,b1\n"
            "treated_rep2,treated,b2\n"
        ),
    )

    check_experiment_design.main(
        _args("--sample_sheet", str(sample_sheet), "--covariates", "batch")
    )


def test_main_rejects_duplicate_aliases(tmp_path):
    """Duplicate aliases should fail validation."""
    sample_sheet = _write(
        tmp_path / "sample_sheet.csv",
        (
            "alias,condition\n"
            "rep1,control\n"
            "rep1,control\n"
            "rep3,treated\n"
            "rep4,treated\n"
        ),
    )

    with pytest.raises(SystemExit, match="aliases must be unique"):
        check_experiment_design.main(_args("--sample_sheet", str(sample_sheet)))


def test_main_requires_reference_level_when_control_is_absent(tmp_path):
    """Non-control condition labels require an explicit reference level."""
    sample_sheet = _write(
        tmp_path / "sample_sheet.csv",
        (
            "alias,condition\n"
            "rep1,baseline\n"
            "rep2,baseline\n"
            "rep3,treated\n"
            "rep4,treated\n"
        ),
    )

    with pytest.raises(SystemExit, match="Provide --reference_level"):
        check_experiment_design.main(_args("--sample_sheet", str(sample_sheet)))
