"""Required sample sheet validators.

These should be applied to sample sheets from all workflows.
"""

from abc import ABC
import re


# Define the base class for sample sheet validators
class SampleSheetValidator(ABC):
    """Base class for sample sheet validators."""

    def __init__(self, wf_params, options=None):
        """Initialize the validator with workflow parameters."""
        self.wf_params = wf_params
        self.options = options or {}
        self.alias_field = "sample_name" if self.options.get("no_barcode") else "alias"
        self.errors = []

    def log_error(self, msg, column=None, lineno=None):
        """Log an error message with optional column and line context."""
        context = []
        if column is not None:
            context.append(f"column: {column}")
        if lineno is not None:
            context.append(f"line: {lineno}")
        if context:
            msg = f"{msg} ({', '.join(context)})"
        self.errors.append(msg)

    def on_header(self, header):
        """Handle header line."""
        pass

    def add_sheet_row(self, row, lineno):
        """Handle a single row from the sample sheet."""
        pass

    @property
    def is_valid(self):
        """Check if the sample sheet is valid according to this validator."""
        if self.errors:
            return False
        return True


class ProhibitedColumns(SampleSheetValidator):
    """Check for prohibited columns."""

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)

    def on_header(self, header):
        """Don't allow `barcode` and `alias` when `--no_barcode` set."""
        if self.options.get("no_barcode"):
            for field in ["alias", "barcode"]:
                if field in header:
                    self.log_error(
                        "Column must not be present with --no_barcode",
                        column=field,
                    )


class RequiredColumns(SampleSheetValidator):
    """Check for required columns."""

    def on_header(self, header):
        """Check for required columns."""
        self.header = header

        required = [self.alias_field]
        if not self.options.get("no_barcode"):
            # in barcode mode required are alias_field and barcode
            required.insert(0, "barcode")

        for req in required:
            if req not in header:
                self.log_error("Column missing", column=req)

    def add_sheet_row(self, row, lineno):
        """Check for consistent number of columns."""
        if len(row) != len(self.header):
            self.log_error("Unexpected number of cells in row", lineno=lineno)


class SampleTypeRules(SampleSheetValidator):
    """Check sample types."""

    ALLOWED = {
        "test_sample", "positive_control", "negative_control", "no_template_control"
    }

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        self.types = []
        self.unexpected_types = []
        self.required = wf_params.get("required_sample_types", [])

    def add_sheet_row(self, row, lineno):
        """Collect sample types."""
        t = row.get("type")
        if t:
            self.types.append(t)
            if t not in self.ALLOWED:
                self.unexpected_types.append((t, lineno))

    @property
    def is_valid(self):
        """Check if the sample types are valid."""
        valid = True
        if self.unexpected_types:
            for t, lineno in self.unexpected_types:
                self.log_error(
                    f"Unexpected value: {t}", column="type", lineno=lineno)
            valid = False
        for required in self.required:
            if required not in self.ALLOWED:
                self.log_error(f"Not an allowed sample type: {required}")
                valid = False
            elif required not in self.types:
                self.log_error(f"Sample sheet requires at least 1 of '{required}'")
                valid = False
        return valid


class BarcodeRules(SampleSheetValidator):
    """Check barcode rules."""

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        self.first_len = None
        self.barcodes = []

    def add_sheet_row(self, row, lineno):
        """Check barcode rules."""
        if self.options.get("no_barcode"):
            return

        bc = row.get("barcode")
        if bc is None:
            self.log_error("Column missing", column="barcode")
            return
        if bc in self.barcodes:
            self.log_error(
                f"Value not unique: {bc}", column="barcode", lineno=lineno
            )
        self.barcodes.append(bc)
        if not re.match(r"^barcode\d\d+$", bc):
            self.log_error(
                f"Value has incorrect format: {bc}",
                column="barcode",
                lineno=lineno,
            )
        if self.first_len is None:
            self.first_len = len(bc)
        elif len(bc) != self.first_len:
            self.log_error(
                f"Values are different lengths: {bc}",
                column="barcode",
                lineno=lineno,
            )


class AliasRules(SampleSheetValidator):
    """Alias rules."""

    ALIAS_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        self.aliases = set()

    def add_sheet_row(self, row, lineno):
        """Check alias rules."""
        alias = row.get(self.alias_field)
        if alias is None:
            self.log_error("Column missing", column=self.alias_field)
            return
        if alias in self.aliases:
            self.log_error(
                f"Value not unique: {alias}", column=self.alias_field, lineno=lineno,
            )
        self.aliases.add(alias)
        # Specific error message if empty "" to improve user error message
        if not alias:
            self.log_error(
                f"Empty {self.alias_field}. "
                "Allowed values start with letters or numbers "
                "and may contain only letters, numbers, '.', '_' or '-'",
                column=self.alias_field,
                lineno=lineno,
            )
        if not self.ALIAS_PATTERN.match(alias):
            self.log_error(
                f"Invalid value {alias}. Allowed values start with letters or numbers "
                "and may contain only letters, numbers, '.', '_' or '-'",
                column=self.alias_field,
                lineno=lineno,
            )
        if alias.startswith("barcode"):
            self.log_error(
                f"Value must not begin with 'barcode': {alias}",
                column=self.alias_field,
                lineno=lineno,
            )


class AnalysisGroupCompleteness(SampleSheetValidator):
    """Analysis groups."""

    def __init__(self, wf_params, options=None):
        """Initialize with workflow parameters."""
        super().__init__(wf_params, options)
        self.has_analysis_groups = False

    def on_header(self, header):
        """Check for required columns."""
        if "analysis_group" in header:
            self.has_analysis_groups = True

    def add_sheet_row(self, row, lineno):
        """Check analysis group completeness."""
        if self.has_analysis_groups and not row.get("analysis_group"):
            self.log_error(
                "Column exists but needs values in each row",
                column="analysis_group",
                lineno=lineno,
            )
