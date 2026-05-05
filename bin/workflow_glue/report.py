"""Create workflow report for wf-transcriptomes."""

import json
from pathlib import Path

from dominate.tags import div, h3, p, pre, strong
from dominate.util import raw
from ezcharts.components import fastcat
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def _find_table(directory, pattern):
    matches = sorted(Path(directory).glob(pattern))
    return matches[0] if matches else None


def _read_table(path, **kwargs):
    if path is None or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t", **kwargs)


def _cohort_summary(cohort_dir):
    tx_meta = _read_table(Path(cohort_dir) / "transcript_metadata.tsv")
    if tx_meta is None:
        return None, None

    summary = {
        "Transcripts": len(tx_meta),
        "Genes": (
            tx_meta["GENEID"].nunique() if "GENEID" in tx_meta.columns else "N/A"
        ),
    }
    for column in ("newTxClass", "newGeneClass"):
        if column in tx_meta.columns:
            counts = tx_meta[column].fillna("NA").value_counts().head(10)
            class_df = counts.rename_axis(column).reset_index(name="count")
            return (
                pd.DataFrame(summary.items(), columns=["Metric", "Value"]),
                class_df,
            )

    return pd.DataFrame(summary.items(), columns=["Metric", "Value"]), None


def _sample_summaries(samples_dir):
    summaries = {}
    for sample_dir in sorted(Path(samples_dir).iterdir()):
        if not sample_dir.is_dir():
            continue
        tx_meta = _read_table(sample_dir / "transcript_metadata.tsv")
        if tx_meta is None:
            continue
        summaries[sample_dir.name] = pd.DataFrame(
            [
                ("Transcripts", len(tx_meta)),
                (
                    "Genes",
                    (
                        tx_meta["GENEID"].nunique()
                        if "GENEID" in tx_meta.columns
                        else "N/A"
                    ),
                ),
            ],
            columns=["Metric", "Value"],
        )
    return summaries


def _sqanti_tables(sqanti_dir):
    tables = {}
    for summary in sorted(Path(sqanti_dir).rglob("classification_summary.tsv")):
        label = summary.parent.name
        tables[label] = _read_table(summary)
    return tables


def _pychopper_tables(pychopper_dir):
    tables = {}
    for summary in sorted(Path(pychopper_dir).rglob("pychopper_summary.tsv")):
        table = _read_table(summary)
        if table is None or table.empty:
            continue
        label = summary.parent.name.replace("_pychopper_output", "")
        tables[label] = table
    return tables


def _top_results(de_dir, filename, n=20):
    tables = {}
    for contrast_dir in sorted(Path(de_dir).iterdir()):
        if not contrast_dir.is_dir():
            continue
        table = _read_table(contrast_dir / filename)
        if table is None or table.empty:
            continue
        tables[contrast_dir.name] = table.head(n)
    return tables


def _load_bambu_qc(cohort_dir):
    """Load bambu QC statistics JSON."""
    qc_file = Path(cohort_dir) / "bambu_qc_stats.json"
    if qc_file.exists():
        with open(qc_file) as f:
            return json.load(f)
    return None


def _load_de_qc(de_dir):
    """Load DE/DTU QC statistics JSON."""
    qc_file = Path(de_dir) / "de_qc_stats.json"
    if qc_file.exists():
        with open(qc_file) as f:
            return json.load(f)
    return None


def _create_warning_banner(message, level="warning"):
    """Create a styled warning banner."""
    colors = {
        "warning": "#fff3cd",
        "danger": "#f8d7da",
        "info": "#d1ecf1",
    }
    border_colors = {
        "warning": "#ffc107",
        "danger": "#dc3545",
        "info": "#0dcaf0",
    }
    style = (
        "padding: 15px; margin: 10px 0; "
        "background-color: {}; "
        "border-left: 4px solid {}; "
        "border-radius: 4px;".format(
            colors.get(level, colors["warning"]),
            border_colors.get(level, border_colors["warning"]),
        )
    )
    with div(style=style):
        with strong():
            raw(
                "⚠️ "
                if level == "warning"
                else "❌ "
                if level == "danger"
                else "ℹ️ "
            )
        raw(message)


def main(args):
    """Run the report entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Transcriptomes Sequencing report",
        "wf-transcriptomes",
        args.params,
        args.versions,
        args.wf_version,
    )

    with open(args.metadata, "r") as handle:
        metadata = json.load(handle)

    with report.add_section("Workflow overview", "Overview"):
        p(
            "This report summarises joint bambu transcript discovery "
            "and quantification, per-sample transcriptomes, optional "
            "SQANTI3 structural classification, and optional "
            "differential analysis outputs."
        )

    if args.stats:
        with report.add_section("Read summary", "Reads"):
            stats = tuple(args.stats)
            sample_names = tuple(
                item["alias"] for item in metadata if item.get("has_stats")
            )
            if len(stats) == 1:
                stats = stats[0]
                sample_names = sample_names[0] if sample_names else None
            fastcat.SeqSummary(stats, sample_names=sample_names)

    with report.add_section("Sample metadata", "Samples"):
        tabs = Tabs()
        for item in sorted(metadata, key=lambda value: value["alias"]):
            with tabs.add_tab(item["alias"]):
                DataTable.from_pandas(
                    pd.DataFrame.from_dict(item, orient="index", columns=["Value"])
                    .reset_index()
                    .rename(columns={"index": "Field"})
                )

    # Load bambu QC statistics
    bambu_qc = _load_bambu_qc(args.cohort_dir)

    # Add Bambu QC section with warnings
    if bambu_qc:
        with report.add_section("Bambu Quality Control", "Bambu QC"):
            # Check for warnings
            if bambu_qc.get("library_size_warning"):
                _create_warning_banner(
                    f"Library Size Variation: {bambu_qc['library_size_warning']}. "
                    "Large variation (>3x) may affect CPM normalization. "
                    "Consider reviewing per-sample library sizes.",
                    level="warning",
                )

            # Library size statistics
            with h3("Library Size Statistics"):
                lib_stats = pd.DataFrame(
                    [
                        ("Samples analyzed", bambu_qc.get("samples", "N/A")),
                        (
                            "Median library size",
                            "{} reads".format(
                                format(
                                    bambu_qc.get("median_library_size", 0),
                                    ",",
                                )
                            ),
                        ),
                        (
                            "Min library size",
                            "{} reads".format(
                                format(
                                    bambu_qc.get("min_library_size", 0),
                                    ",",
                                )
                            ),
                        ),
                        (
                            "Max library size",
                            "{} reads".format(
                                format(
                                    bambu_qc.get("max_library_size", 0),
                                    ",",
                                )
                            ),
                        ),
                        (
                            "Library size ratio (max/min)",
                            "{:.2f}x".format(
                                bambu_qc.get("library_size_ratio", 1.0)
                            ),
                        ),
                    ],
                    columns=["Metric", "Value"],
                )
                DataTable.from_pandas(lib_stats, paging=False, searchable=False)

            # Transcript discovery statistics
            with h3("Transcript Discovery"):
                discovery_stats = pd.DataFrame(
                    [
                        (
                            "Transcriptome mode",
                            bambu_qc.get("transcriptome_mode", "N/A"),
                        ),
                        ("NDR used", str(bambu_qc.get("ndr_used", "N/A"))),
                        (
                            "Transcripts before filtering",
                            bambu_qc.get("total_transcripts_before_filter", 0),
                        ),
                        (
                            "Transcripts after filtering",
                            bambu_qc.get("total_transcripts_after_filter", 0),
                        ),
                        (
                            "Transcripts removed",
                            bambu_qc.get("transcripts_filtered", 0),
                        ),
                        (
                            "Median transcripts per sample",
                            bambu_qc.get("median_transcripts_detected", 0),
                        ),
                        (
                            "Unique genes (after filter)",
                            bambu_qc.get("total_genes_after_filter", 0),
                        ),
                    ],
                    columns=["Metric", "Value"],
                )
                DataTable.from_pandas(discovery_stats, paging=False, searchable=False)

            # Per-sample library sizes
            if "library_sizes" in bambu_qc and bambu_qc["library_sizes"]:
                with h3("Per-Sample Library Sizes"):
                    lib_size_data = []
                    for sample, size in bambu_qc["library_sizes"].items():
                        lib_size_data.append(
                            {
                                "Sample": sample,
                                "Library Size": format(size, ","),
                                "Reads": size,
                            }
                        )
                    lib_df = pd.DataFrame(lib_size_data).sort_values(
                        "Reads", ascending=False
                    )
                    DataTable.from_pandas(
                        lib_df[["Sample", "Library Size"]],
                        paging=False,
                        use_index=False,
                    )

    with report.add_section("Cohort transcriptome", "Cohort"):
        cohort_metrics, cohort_classes = _cohort_summary(args.cohort_dir)
        if cohort_metrics is not None:
            DataTable.from_pandas(cohort_metrics, paging=False, searchable=False)
        if cohort_classes is not None:
            DataTable.from_pandas(cohort_classes, paging=False, searchable=False)

        tx_counts = _read_table(Path(args.cohort_dir) / "transcript_counts.tsv")
        if tx_counts is not None and not tx_counts.empty:
            p("Top transcript rows from the cohort abundance table.")
            DataTable.from_pandas(tx_counts.head(20), use_index=False)

    with report.add_section("Per-sample transcriptomes", "Per sample"):
        tabs = Tabs()
        for sample, summary_df in _sample_summaries(args.samples_dir).items():
            with tabs.add_tab(sample):
                DataTable.from_pandas(summary_df, paging=False, searchable=False)

    if args.alignment_stats_dir and Path(args.alignment_stats_dir).exists():
        with report.add_section("Alignment statistics", "Alignments"):
            tabs = Tabs()
            for stats_file in sorted(
                Path(args.alignment_stats_dir).glob("*.flagstat.txt")
            ):
                with tabs.add_tab(stats_file.stem.replace(".flagstat", "")):
                    pre(stats_file.read_text())

    pychopper_tables = {}
    if args.pychopper_dir and Path(args.pychopper_dir).exists():
        pychopper_tables = _pychopper_tables(args.pychopper_dir)
    if pychopper_tables:
        with report.add_section("Pychopper preprocessing", "Pychopper"):
            tabs = Tabs()
            for label, table in pychopper_tables.items():
                with tabs.add_tab(label):
                    DataTable.from_pandas(table, use_index=False)

    sqanti_tables = _sqanti_tables(args.sqanti_dir)
    if sqanti_tables:
        with report.add_section("SQANTI3 classification", "SQANTI3"):
            tabs = Tabs()
            for label, table in sqanti_tables.items():
                with tabs.add_tab(label):
                    DataTable.from_pandas(table, use_index=False)

    if args.de_dir and Path(args.de_dir).exists():
        # Load DE QC statistics
        de_qc = _load_de_qc(args.de_dir)

        # Add DE/DTU QC section with warnings
        if de_qc:
            with report.add_section(
                "Differential Analysis Quality Control",
                "DE/DTU QC",
            ):
                # Check for critical warnings
                has_warnings = False

                if (
                    de_qc.get("sample_size_warnings")
                    and de_qc["sample_size_warnings"] != "none"
                ):
                    _create_warning_banner(
                        f"Sample Size Warning: {de_qc['sample_size_warnings']}. "
                        "Underpowered designs may have reduced statistical "
                        "power and increased false negative rate.",
                        level="warning",
                    )
                    has_warnings = True

                if de_qc.get("multiple_testing_note"):
                    _create_warning_banner(
                        de_qc["multiple_testing_note"]
                        + ". See MULTIPLE_TESTING_WARNING.txt for details.",
                        level="info",
                    )

                # Check for dispersion fallbacks
                dispersion_fallbacks = []
                for contrast_name, contrast_data in de_qc.get(
                    "contrasts", {}
                ).items():
                    dispersion_file = (
                        Path(args.de_dir)
                        / f"DESeq2_dispersion_fallback_{contrast_name}.txt"
                    )
                    if dispersion_file.exists():
                        dispersion_fallbacks.append(contrast_name)

                if dispersion_fallbacks:
                    _create_warning_banner(
                        "Dispersion Estimation Fallback: "
                        f"{len(dispersion_fallbacks)} contrast(s) used "
                        "gene-wise dispersion (reduced power). "
                        f"Affected: {', '.join(dispersion_fallbacks)}",
                        level="warning",
                    )
                    has_warnings = True

                # Check for failed DTU analyses
                failed_dtu = []
                for contrast_name, contrast_data in de_qc.get(
                    "contrasts", {}
                ).items():
                    if contrast_data.get("dtu_status") == "FAILED":
                        failed_dtu.append(contrast_name)

                if failed_dtu:
                    _create_warning_banner(
                        "DTU Analysis Failed: "
                        f"{len(failed_dtu)} contrast(s) could not perform "
                        "DTU testing. "
                        f"Affected: {', '.join(failed_dtu)}. "
                        "See DTU_ANALYSIS_FAILED.txt files for details.",
                        level="danger",
                    )
                    has_warnings = True

                # Experimental design summary
                with h3("Experimental Design"):
                    covariates = de_qc.get("covariates", [])
                    covariates_value = (
                        ", ".join(covariates)
                        if de_qc.get("covariates") != "none"
                        else "none"
                    )
                    design_stats = pd.DataFrame(
                        [
                            ("Total samples", de_qc.get("total_samples", 0)),
                            (
                                "Condition column",
                                de_qc.get("condition_column", "N/A"),
                            ),
                            (
                                "Reference level",
                                de_qc.get("reference_level", "N/A"),
                            ),
                            ("Covariates", covariates_value),
                            (
                                "Number of contrasts",
                                de_qc.get("num_contrasts", 0),
                            ),
                        ],
                        columns=["Parameter", "Value"],
                    )
                    DataTable.from_pandas(design_stats, paging=False, searchable=False)

                # Sample sizes per group
                if "samples_per_group" in de_qc:
                    with h3("Sample Sizes per Group"):
                        sample_size_data = []
                        for group, count in de_qc["samples_per_group"].items():
                            status = (
                                "✓" if count >= 3 else "⚠️" if count >= 2 else "❌"
                            )
                            note = (
                                "OK"
                                if count >= 3
                                else "Low power"
                                if count >= 2
                                else "Too few"
                            )
                            sample_size_data.append(
                                {
                                    "Group": group,
                                    "Samples": count,
                                    "Status": status,
                                    "Note": note,
                                }
                            )
                        sample_df = pd.DataFrame(sample_size_data)
                        DataTable.from_pandas(
                            sample_df,
                            paging=False,
                            use_index=False,
                        )

                # Per-contrast summary
                if "contrasts" in de_qc:
                    with h3("Results Summary by Contrast"):
                        contrast_summary_data = []
                        for contrast_name, contrast_data in de_qc["contrasts"].items():
                            dtu_genes = (
                                contrast_data.get("dtu_significant_genes", 0)
                                if contrast_data.get("dtu_status") == "SUCCESS"
                                else "N/A"
                            )
                            contrast_summary_data.append(
                                {
                                    "Contrast": contrast_name,
                                    "Samples": (
                                        f"{contrast_data.get('n_target', 0)} "
                                        f"vs "
                                        f"{contrast_data.get('n_reference', 0)}"
                                    ),
                                    "DGE Significant (FDR<0.05)": (
                                        contrast_data.get(
                                            "dge_significant_fdr05", 0
                                        )
                                    ),
                                    "DGE Up": contrast_data.get(
                                        "dge_upregulated", 0
                                    ),
                                    "DGE Down": contrast_data.get(
                                        "dge_downregulated", 0
                                    ),
                                    "DTU Status": contrast_data.get(
                                        "dtu_status", "N/A"
                                    ),
                                    "DTU Genes (q<0.05)": dtu_genes,
                                }
                            )
                        contrast_summary_df = pd.DataFrame(contrast_summary_data)
                        DataTable.from_pandas(
                            contrast_summary_df,
                            paging=False,
                            use_index=False,
                        )

                # Warnings summary table
                if has_warnings:
                    with h3("Quality Warnings Summary"):
                        warnings_data = []
                        if (
                            de_qc.get("sample_size_warnings")
                            and de_qc["sample_size_warnings"] != "none"
                        ):
                            warnings_data.append(
                                {
                                    "Warning Type": "Sample Size",
                                    "Details": de_qc["sample_size_warnings"],
                                }
                            )
                        if dispersion_fallbacks:
                            warnings_data.append(
                                {
                                    "Warning Type": "Dispersion Estimation",
                                    "Details": (
                                        f"{len(dispersion_fallbacks)} "
                                        "contrasts affected"
                                    ),
                                }
                            )
                        if failed_dtu:
                            warnings_data.append(
                                {
                                    "Warning Type": "DTU Failure",
                                    "Details": (
                                        f"{len(failed_dtu)} contrasts failed"
                                    ),
                                }
                            )

                        warnings_df = pd.DataFrame(warnings_data)
                        DataTable.from_pandas(
                            warnings_df,
                            paging=False,
                            use_index=False,
                        )

        with report.add_section("Differential gene expression", "DGE"):
            tabs = Tabs()
            for contrast, table in _top_results(args.de_dir, "results_dge.tsv").items():
                with tabs.add_tab(contrast):
                    # Check for contrast-specific warnings
                    if de_qc and contrast in de_qc.get("contrasts", {}):
                        contrast_data = de_qc["contrasts"][contrast]
                        if contrast_data.get("dtu_power_warning"):
                            with div(
                                style=(
                                    "padding: 10px; margin-bottom: 10px; "
                                    "background-color: #fff3cd; "
                                    "border-radius: 4px;"
                                )
                            ):
                                with p():
                                    strong("Note: ")
                                    raw(contrast_data["dtu_power_warning"])

                    DataTable.from_pandas(table, use_index=False)

        with report.add_section("Differential transcript usage", "DTU"):
            tabs = Tabs()
            dtu_tables = _top_results(args.de_dir, "results_dtu_transcript.tsv")

            for contrast in sorted(Path(args.de_dir).iterdir()):
                if not contrast.is_dir():
                    continue
                contrast_name = contrast.name

                with tabs.add_tab(contrast_name):
                    # Check if DTU failed for this contrast
                    if de_qc and contrast_name in de_qc.get("contrasts", {}):
                        contrast_data = de_qc["contrasts"][contrast_name]
                        if contrast_data.get("dtu_status") == "FAILED":
                            _create_warning_banner(
                                "DTU analysis failed for this contrast. "
                                f"See {contrast_name}/"
                                "DTU_ANALYSIS_FAILED.txt for detailed "
                                "explanation.",
                                level="danger",
                            )
                            p(
                                "Empty results indicate analysis failure, "
                                "not 'no DTU detected'."
                            )
                        elif contrast_data.get("dtu_power_warning"):
                            _create_warning_banner(
                                contrast_data["dtu_power_warning"],
                                level="warning",
                            )

                    # Show table if available
                    if contrast_name in dtu_tables:
                        DataTable.from_pandas(
                            dtu_tables[contrast_name],
                            use_index=False,
                        )
                    else:
                        p("No DTU results available for this contrast.")

    report.write(args.report)
    logger.info("Report written to %s.", args.report)


def argparser():
    """Argument parser for the report entry point."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file.")
    parser.add_argument("--metadata", required=True, help="Sample metadata JSON.")
    parser.add_argument("--stats", nargs="+", help="Per-read stats paths.")
    parser.add_argument(
        "--alignment_stats_dir",
        default=None,
        help="Alignment stats directory.",
    )
    parser.add_argument(
        "--cohort_dir",
        required=True,
        help="Cohort output directory.",
    )
    parser.add_argument(
        "--samples_dir",
        required=True,
        help="Per-sample output directory.",
    )
    parser.add_argument(
        "--pychopper_dir",
        default=None,
        help="Pychopper output directory.",
    )
    parser.add_argument(
        "--sqanti_dir",
        required=True,
        help="SQANTI output directory.",
    )
    parser.add_argument(
        "--de_dir",
        default=None,
        help="Differential analysis directory.",
    )
    parser.add_argument("--versions", required=True, help="Versions directory.")
    parser.add_argument("--params", required=True, help="Workflow params JSON.")
    parser.add_argument(
        "--wf_version",
        default="unknown",
        help="Workflow version.",
    )
    return parser
