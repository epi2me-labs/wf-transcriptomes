"""Summarise SQANTI3 classification output into stable TSVs."""

from pathlib import Path

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def _find_one(directory, patterns):
    for pattern in patterns:
        matches = sorted(directory.glob(pattern))
        if matches:
            return matches[0]
    return None


def main(args):
    """Extract a simple classification summary from SQANTI output."""
    logger = get_named_logger("summSQANTI")
    sqanti_dir = Path(args.sqanti_dir)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    classification = _find_one(
        sqanti_dir,
        ["*_classification.txt", "*classification.txt", "*_classification.tsv"]
    )
    if classification is None:
        raise SystemExit(
            f"Could not find a SQANTI classification table in {sqanti_dir}."
        )

    df = pd.read_csv(classification, sep="\t")
    structural_column = None
    for name in ("structural_category", "category"):
        if name in df.columns:
            structural_column = name
            break
    if structural_column is None:
        raise SystemExit(
            f"Could not find a structural category column in {classification}."
        )

    summary = (
        df[structural_column]
        .value_counts(dropna=False)
        .rename_axis("structural_category")
        .reset_index(name="count")
        .sort_values(["count", "structural_category"], ascending=[False, True])
    )
    summary.to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote SQANTI summary to %s.", out_path)


def argparser():
    """Argument parser for the SQANTI summariser."""
    parser = wf_parser("summarise_sqanti")
    parser.add_argument("--sqanti_dir", required=True, help="SQANTI run directory.")
    parser.add_argument("--output", required=True, help="Summary TSV output path.")
    return parser
