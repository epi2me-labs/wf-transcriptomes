"""Workflow-specific sample sheet validators."""

from collections import Counter
import re

from .default import SampleSheetValidator  # noqa: ABS101

_R_FORMULA_NAME_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_.]*$")


def _split_covariates(value):
    """Normalise workflow covariates into a clean list."""
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]


def _validate_r_formula_names(names, label="Column"):
    """Validate names against the shared design-name pattern."""
    empty = names.count("")
    invalid = [
        name for name in names
        if name != "" and (
            not isinstance(name, str) or _R_FORMULA_NAME_RE.match(name) is None
        )
    ]
    if empty or invalid:
        details = []
        if empty:
            details.append('Empty value: ""')
        if invalid:
            details.append(f"Invalid names: {', '.join(invalid)}")
        raise ValueError(
            f"{label} names must be safe for R formulas. {'; '.join(details)}. "
            "Names must start with a letter and contain "
            "only letters, numbers, underscores, and dots."
        )


class Covariates(SampleSheetValidator):
    """Validate optional DE/DTU covariate columns from workflow params."""

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        if wf_params.get("covariates"):
            self.covariates = _split_covariates(wf_params.get("covariates"))
        else:
            self.covariates = []

    def on_header(self, header):
        """Require configured covariate columns when they are requested."""
        if not self.wf_params.get("de_analysis"):
            return
        # Check values of condition are safe for R
        try:
            _validate_r_formula_names(
                self.covariates,
                label="Covariate column",
            )
        except ValueError as exc:
            self.log_error(str(exc))

        missing_covariates = [
            covariate for covariate in self.covariates
            if covariate not in header
        ]
        if missing_covariates:
            self.log_error(
                f"Missing covariate columns: {', '.join(missing_covariates)}"
            )

    def add_sheet_row(self, row, lineno):
        """Require values in requested design columns for each row."""
        if not self.wf_params.get("de_analysis"):
            return

        for covariate in self.covariates:
            if covariate in row and not row.get(covariate):
                self.log_error(
                    "Covariate column must not contain missing or empty values",
                    column=covariate,
                    lineno=lineno,
                )
                continue

            if covariate in row:
                try:
                    _validate_r_formula_names(
                        [row.get(covariate)],
                        label="Covariate value",
                    )
                except ValueError as exc:
                    self.log_error(
                        str(exc),
                        column=covariate,
                        lineno=lineno,
                    )


class ConditionReplicates(SampleSheetValidator):
    """Validate the condition column and check that each level has replicates."""

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        self.condition_column = wf_params.get("condition_column")
        self.has_condition = False
        self.condition_counts = Counter()

    def on_header(self, header):
        """Validate and track whether the sample sheet contains a condition column."""
        if not self.wf_params.get("de_analysis"):
            return

        try:
            _validate_r_formula_names(
                [self.condition_column],
                label="Design column",
            )
        except ValueError as exc:
            self.log_error(str(exc))

        self.has_condition = (
            self.condition_column is not None and self.condition_column in header
        )

        if self.condition_column and self.condition_column not in header:
            self.log_error(
                f"Sample sheet must contain the '{self.condition_column}' column.",
                column=self.condition_column,
            )

    def add_sheet_row(self, row, lineno):
        """Validate condition values and count samples for each condition."""
        if not self.has_condition:
            return

        condition = row.get(self.condition_column)
        if not condition:
            self.log_error(
                "Condition column must not contain missing or empty values",
                column=self.condition_column,
                lineno=lineno,
            )
            return
        # Check values of condition are safe for R
        try:
            _validate_r_formula_names(
                [condition],
                label="Condition value",
            )
        except ValueError as exc:
            self.log_error(str(exc), column=self.condition_column, lineno=lineno)

        self.condition_counts[condition] += 1

    @property
    def is_valid(self):
        """Check condition levels and replicates are sufficient for DE."""
        valid = True
        if self.has_condition and len(self.condition_counts) < 2:
            self.log_error(
                "Condition column must contain at least 2 condition levels",
                column=self.condition_column,
            )
            valid = False
        for condition, count in self.condition_counts.items():
            if count < 2:
                self.log_error(
                    f"Condition must have at least 2 samples: {condition}",
                    column=self.condition_column,
                )
                valid = False
        return valid and super().is_valid


class ReferenceLevel(SampleSheetValidator):
    """Validate the DE reference level against observed condition values."""

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        self.condition_column = wf_params.get("condition_column")
        self.has_condition = False
        self.condition_values = set()

    def on_header(self, header):
        """Track whether reference-level validation should run."""
        self.has_condition = (
            self.wf_params.get("de_analysis")
            and self.condition_column is not None
            and self.condition_column in header
        )

    def add_sheet_row(self, row, lineno):
        """Collect condition values for reference-level validation."""
        if not self.has_condition:
            return

        condition = row.get(self.condition_column)
        if not condition:
            return

        self.condition_values.add(condition)

    @property
    def is_valid(self):
        """Validate the configured or inferred reference level."""
        if not self.has_condition:
            return super().is_valid

        reference_level = self.wf_params.get("reference_level")
        if reference_level is None:
            if "control" in self.condition_values:
                reference_level = "control"
            else:
                self.log_error(
                    (
                        "Provide --reference_level when the condition column "
                        "does not match the default 'control'."
                    )
                )
                return False

        if reference_level not in self.condition_values:
            self.log_error((
                f"The requested reference level '{reference_level}' is not present "
                "in the condition column."))
            return False

        return super().is_valid
