"""Test check_sample_sheet.py."""

import json

import pytest
from workflow_glue.wfg_helpers import check_sample_sheet
from workflow_glue.wfg_helpers.validators.wf import _validate_r_formula_names


def _run_check_sample_sheet(sample_sheet_path, params_json, capsys):
    """Run sample-sheet validation and return captured stdout."""
    args = [str(sample_sheet_path), str(params_json), "--no_barcode"]
    parsed_args = check_sample_sheet.argparser().parse_args(args)

    try:
        check_sample_sheet.main(parsed_args)
    except SystemExit:
        pass

    out, _ = capsys.readouterr()
    return out


@pytest.mark.parametrize(
    "names",
    [
        [
            "condition", "batch_1", "site.2",
            "S01", "a1", "sample_group"
        ],
    ],
)
def test_validate_r_formula_names_accepts_valid_names(names):
    """R formula names accept letter-led names with safe punctuation."""
    _validate_r_formula_names(names)


@pytest.mark.parametrize(
    "names",
    [
        ["_condition"],
        ["condition 1"],
        ["size;batch"],
        ["size:batch"],
        ["size|batch"],
        ["sueño"],
        ["01"],
        ["1batch"],
        ["2.site-3"],
        [""]
    ],
)
def test_validate_r_formula_names_rejects_invalid_names(names):
    """Invalid R formula names should raise a validation error."""
    with pytest.raises(ValueError):
        _validate_r_formula_names(names)


def test_check_sample_sheet_skips_condition_check_without_condition_column_param(
        tmp_path, capsys
):
    """Test that condition checks are skipped without the param."""
    # This would fail (as there is no enough replicates)
    # but in this case is ok as de_analysis is not enabled
    sample_sheet_path = tmp_path / "sample_sheet_condition_counts.csv"
    sample_sheet_path.write_text(
        "sample_name,condition\n"
        "test_name1,condition1\n"
        "test_name2,condition2\n"
        "test_name3,condition2\n"
    )

    params_json = tmp_path / "params.json"
    # condition has a default from the nextflow schema
    params_json.write_text(
        json.dumps({"de_analysis": False, "condition_column": "condition"})
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out == ""


def test_check_sample_sheet_rejects_condition_with_single_sample_when_param_present(
        tmp_path, capsys
):
    """Test that each condition must have at least two samples when DE is enabled."""
    sample_sheet_path = tmp_path / "sample_sheet_condition_counts.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition\n"
        "test_name1,test_sample,condition1\n"
        "test_name2,positive_control,condition2\n"
        "test_name3,negative_control,condition2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({"de_analysis": True, "condition_column": "condition"})
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith(
        "Condition must have at least 2 samples: condition1 "
        "(column: condition)"
    )


def test_check_sample_sheet_requires_at_least_two_condition_levels(
        tmp_path, capsys
):
    """DE-enabled sample sheets must contain at least two condition levels."""
    sample_sheet_path = tmp_path / "sample_sheet_single_condition_level.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition\n"
        "test_name1,test_sample,condition1\n"
        "test_name2,positive_control,condition1\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({"de_analysis": True, "condition_column": "condition"})
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith(
        "Condition column must contain at least 2 condition levels "
        "(column: condition)"
    )


def test_check_sample_sheet_rejects_missing_condition_values(
        tmp_path, capsys
):
    """Sample sheets must not contain empty condition cells."""
    sample_sheet_path = tmp_path / "sample_sheet_missing_condition_value.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition\n"
        "test_name1,test_sample,\n"
        "test_name2,positive_control,treated\n"
        "test_name3,negative_control,treated\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({"de_analysis": True, "condition_column": "condition"})
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith(
        "Condition column must not contain missing or empty values "
        "(column: condition, line: 1)"
    )


def test_check_sample_sheet_requires_condition_column_when_param_present(
        tmp_path, capsys
):
    """Configured condition columns must exist when DE is enabled."""
    sample_sheet_path = tmp_path / "sample_sheet_missing_condition.csv"
    sample_sheet_path.write_text(
        "sample_name,batch\n"
        "test_name1,b1\n"
        "test_name2,b2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({"de_analysis": True, "condition_column": "condition"})
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith(
        "Sample sheet must contain the 'condition' column. "
        "(column: condition)"
    )


def test_check_sample_sheet_requires_covariate_columns_when_present(
        tmp_path, capsys
):
    """Check covariate columns exist when provided and DE is enabled."""
    sample_sheet_path = tmp_path / "sample_sheet_missing_covariate.csv"
    sample_sheet_path.write_text(
        "sample_name,condition\n"
        "test_name1,control\n"
        "test_name2,control\n"
        "test_name3,treated\n"
        "test_name4,treated\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith("Missing covariate columns: batch, site")


def test_check_sample_sheet_rejects_missing_covariate_values(
        tmp_path, capsys
):
    """Configured covariates must have a value in each row when DE is enabled."""
    sample_sheet_path = tmp_path / "sample_sheet_missing_covariate_value.csv"
    sample_sheet_path.write_text(
        "sample_name,condition,batch,site\n"
        "test_name1,control,,s1\n"
        "test_name2,control,b1,s2\n"
        "test_name3,treated,b2,s1\n"
        "test_name4,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith(
        "Covariate column must not contain missing or empty values "
        "(column: batch, line: 1)"
    )


def test_check_sample_sheet_rejects_unsafe_covariate_values(
        tmp_path, capsys
):
    """Configured covariate values must use safe contrast/design characters."""
    sample_sheet_path = tmp_path / "sample_sheet_invalid_covariate_value.csv"
    sample_sheet_path.write_text(
        "sample_name,condition,good_batch,site\n"
        "test_name1,control,b1,s1\n"
        "test_name2,control,b2,s2\n"
        "test_name3,treated,b 1,s1\n"
        "test_name4,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "good_batch,site",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith(
        "Covariate value names must be safe for R formulas. "
        "Invalid names: b 1. Names must start with a letter and "
        "contain only letters, numbers, underscores, and dots. "
        "(column: good_batch, line: 3)"
    )


def test_check_sample_sheet_rejects_unsafe_column_names(
        tmp_path, capsys
):
    """Configured DE design names may not include hyphens or spaces."""
    sample_sheet_path = tmp_path / "sample_sheet_invalid_design_name.csv"
    sample_sheet_path.write_text(
        "sample_name,condition-column,bad covariate\n"
        "test_name1,control,b1\n"
        "test_name2,control,b2\n"
        "test_name3,treated,b1\n"
        "test_name4,treated,b2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition-column",
            "covariates": "bad covariate",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)
    assert out.splitlines() == [
        "Design column names must be safe for R formulas. "
        "Invalid names: condition-column. Names must start with a letter and "
        "contain only letters, numbers, underscores, and dots.",
        "Covariate column names must be safe for R formulas. "
        "Invalid names: bad covariate. Names must start with a letter and "
        "contain only letters, numbers, underscores, and dots.",
    ]


def test_check_sample_sheet_rejects_unsafe_condition_values(
        tmp_path, capsys
):
    """Condition values must be safe for contrast output paths."""
    sample_sheet_path = tmp_path / "sample_sheet_invalid_condition_value.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition,batch,site\n"
        "test_name1,test_sample,control/group,b1,s1\n"
        "test_name2,positive_control,control/group,b1,s2\n"
        "test_name3,test_sample,treated,b2,s1\n"
        "test_name4,negative_control,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.splitlines() == [
        "Condition value names must be safe for R formulas. "
        "Invalid names: control/group. Names must start with a letter and "
        "contain only letters, numbers, underscores, and dots. "
        "(column: condition, line: 1)",
        "Condition value names must be safe for R formulas. "
        "Invalid names: control/group. Names must start with a letter and "
        "contain only letters, numbers, underscores, and dots. "
        "(column: condition, line: 2)",
    ]


def test_check_sample_sheet_accepts_control_as_default_reference_level(
        tmp_path, capsys
):
    """Control is accepted as the default reference level when present."""
    sample_sheet_path = tmp_path / "sample_sheet_reference_control.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition,batch,site\n"
        "test_name1,test_sample,control,b1,s1\n"
        "test_name2,positive_control,control,b1,s2\n"
        "test_name3,test_sample,treated,b2,s1\n"
        "test_name4,negative_control,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out == ""


def test_check_sample_sheet_requires_reference_level_without_control(
        tmp_path, capsys
):
    """A non-control design requires an explicit reference level."""
    sample_sheet_path = tmp_path / "sample_sheet_missing_reference.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition,batch,site\n"
        "test_name1,test_sample,baseline,b1,s1\n"
        "test_name2,positive_control,baseline,b1,s2\n"
        "test_name3,test_sample,treated,b2,s1\n"
        "test_name4,negative_control,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)
    assert out == (
        "Provide --reference_level when the condition column "
        "does not match the default 'control'.\n"
    )


def test_check_sample_sheet_accepts_requested_reference_level(
        tmp_path, capsys
):
    """A configured reference level is accepted when present."""
    sample_sheet_path = tmp_path / "sample_sheet_valid_reference.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition,batch,site\n"
        "test_name1,test_sample,baseline,b1,s1\n"
        "test_name2,positive_control,baseline,b1,s2\n"
        "test_name3,test_sample,treated,b2,s1\n"
        "test_name4,negative_control,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
            "reference_level": "baseline",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out == ""


def test_check_sample_sheet_rejects_missing_requested_reference_level(
        tmp_path, capsys
):
    """Configured reference levels must exist in the condition column."""
    sample_sheet_path = tmp_path / "sample_sheet_invalid_reference.csv"
    sample_sheet_path.write_text(
        "sample_name,type,condition,batch,site\n"
        "test_name1,test_sample,baseline,b1,s1\n"
        "test_name2,positive_control,baseline,b1,s2\n"
        "test_name3,test_sample,treated,b2,s1\n"
        "test_name4,negative_control,treated,b2,s2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": True,
            "condition_column": "condition",
            "covariates": "batch,site",
            "reference_level": "control",
        })
    )

    out = _run_check_sample_sheet(sample_sheet_path, params_json, capsys)

    assert out.startswith((
        "The requested reference level 'control' is not present in "
        "the condition column.")
    )


def test_check_sample_sheet_skips_de_checks_when_de_analysis_is_disabled(
        tmp_path, capsys
):
    """DE design columns are optional unless DE analysis is enabled."""
    sample_sheet_path = tmp_path / "sample_sheet_no_condition.csv"
    sample_sheet_path.write_text(
        "barcode,alias\n"
        "barcode01,test_name1\n"
        "barcode02,test_name2\n"
    )

    params_json = tmp_path / "params.json"
    params_json.write_text(
        json.dumps({
            "de_analysis": False,
            "condition_column": "condition",
            "covariates": "batch",
        })
    )

    args = [str(sample_sheet_path), str(params_json)]
    parsed_args = check_sample_sheet.argparser().parse_args(args)

    try:
        check_sample_sheet.main(parsed_args)
    except SystemExit:
        pass

    out, _ = capsys.readouterr()

    assert out == ""
