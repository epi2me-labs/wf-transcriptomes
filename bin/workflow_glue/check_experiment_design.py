"""Validate DE/DTU sample sheet settings."""

from collections import Counter
import csv
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


def _split_covariates(value):
    if not value:
        return []
    return [part.strip() for part in value.split(",") if part.strip()]


def main(args):
    """Validate sample sheet content for DE/DTU."""
    logger = get_named_logger("checkDesign")
    covariates = _split_covariates(args.covariates)
    with open(args.sample_sheet, "r", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            sys.exit("Sample sheet is empty.")
        fieldnames = set(reader.fieldnames)
        required = {"alias", args.condition_column}
        missing = sorted(required - fieldnames)
        if missing:
            sys.exit(
                "Sample sheet is missing required columns: "
                + ", ".join(missing)
            )

        missing_covariates = [name for name in covariates if name not in fieldnames]
        if missing_covariates:
            sys.exit(
                "Sample sheet is missing requested covariate columns: "
                + ", ".join(missing_covariates)
            )

        aliases = []
        levels = Counter()
        for row in reader:
            alias = row.get("alias", "").strip()
            if not alias:
                sys.exit("Sample sheet contains a row with an empty alias value.")
            aliases.append(alias)

            value = row.get(args.condition_column, "").strip()
            if not value:
                sys.exit(
                    f"Sample sheet contains a row with an empty "
                    f"'{args.condition_column}' value."
                )
            levels[value] += 1

            for covariate in covariates:
                if not row.get(covariate, "").strip():
                    sys.exit(
                        f"Sample sheet contains an empty value in covariate column "
                        f"'{covariate}'."
                    )

    duplicate_aliases = [
        alias for alias, count in Counter(aliases).items() if count > 1
    ]
    if duplicate_aliases:
        sys.exit(
            "Sample sheet aliases must be unique. Duplicates: "
            + ", ".join(sorted(duplicate_aliases))
        )

    if len(levels) < 2:
        sys.exit(
            f"The condition column '{args.condition_column}' must contain at "
            "least two levels."
        )

    reference_level = args.reference_level
    if not reference_level:
        if "control" in levels:
            reference_level = "control"
        else:
            sys.exit(
                "Provide --reference_level when the condition column does not "
                "contain 'control'."
            )

    if reference_level not in levels:
        sys.exit(
            f"Reference level '{reference_level}' was absent from "
            f"'{args.condition_column}'."
        )

    underpowered = [level for level, count in levels.items() if count < 2]
    if underpowered:
        sys.exit(
            "Each condition level must contain at least two samples. "
            "Levels with too few samples: "
            + ", ".join(sorted(underpowered))
        )

    logger.info(
        "Validated sample sheet %s using condition column '%s' with reference '%s'.",
        args.sample_sheet,
        args.condition_column,
        reference_level,
    )


def argparser():
    """Argument parser for the validation entry point."""
    parser = wf_parser("check_experiment_design")
    parser.add_argument("--sample_sheet", required=True, help="Sample sheet CSV.")
    parser.add_argument(
        "--condition_column", default="condition",
        help="Primary biological variable column."
    )
    parser.add_argument(
        "--covariates", default=None,
        help="Comma-separated nuisance covariates."
    )
    parser.add_argument(
        "--reference_level", default=None,
        help="Reference level of the primary condition column."
    )
    return parser
