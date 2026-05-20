"""Sample sheet validators.

Add workflow-specific validators by adding modules to this package. Each module
can define one or more concrete subclasses of `SampleSheetValidator`; these will be
discovered automatically by `check_sample_sheet`.
"""

from .default import SampleSheetValidator  # noqa: ABS101

__all__ = ["SampleSheetValidator"]
