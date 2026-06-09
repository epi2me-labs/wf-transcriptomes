"""Create workflow report for wf-transcriptomes."""

from html import escape
import json
import math
import os
from pathlib import Path
import warnings

from bokeh.resources import INLINE as BOKEH_INLINE
from dominate.tags import (
    br, div, h3, h4, p, pre, script, small, strong, style as dom_style
)
from dominate.util import raw
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.components.theme import LAB_head_resources
from ezcharts.layout.resource import Resource as EZC_Resource
from ezcharts.layout.snippets import Tabs
from ezcharts.layout.snippets.table import DataTable
import pandas as pd
from .hierarchical_clustering import hierarchical, clustering_info   # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101
from .volcano import volcano  # noqa: ABS101

logger = get_named_logger("Report")


# Suppress asyncio deprecation warning triggered by dominate on Python 3.10+.
# dominate calls asyncio.get_event_loop() outside a running async context.
warnings.filterwarnings(
    "ignore",
    message="There is no current event loop",
    category=DeprecationWarning,
)

classification_categories = {
    "Full splice match": (
        "Reference and query isoforms have the same number of exons and "
        "all internal junctions agree."
    ),
    "Incomplete splice match": (
        "Query isoform has fewer 5&prime; exons than the reference, with "
        "matching internal junctions."
    ),
    "Novel in catalog": (
        "No full or incomplete splice match, but uses a "
        "combination of known donor/acceptor splice sites."
    ),
    "Novel not in catalog": (
        "No full or incomplete splice match, with at least "
        "one unannotated donor or acceptor splice site."
    ),
    "Antisense": (
        "No same-strand reference overlap, but antisense to an "
        "annotated gene."
    ),
    "Genic intron": (
        "Query isoform is fully contained within an annotated intron."
    ),
    "Genic": "Query isoform overlaps introns and exons.",
    "Intergenic": "Query isoform lies in an intergenic region.",
}


def get_bokeh_widgets_js():
    """Return the inline Bokeh widgets JavaScript bundle."""
    widgets_index = BOKEH_INLINE.components_for("js").index("bokeh-widgets")
    return raw(BOKEH_INLINE.js_raw[widgets_index])


def get_bokeh_tables_js():
    """Return the inline Bokeh tables JavaScript bundle."""
    tables_index = BOKEH_INLINE.components_for("js").index("bokeh-tables")
    return raw(BOKEH_INLINE.js_raw[tables_index])


def _read_table(path, **kwargs):
    """Read a TSV file into a DataFrame, returning None if path is absent."""
    if path is None or not Path(path).exists():
        return None
    return pd.read_csv(path, sep="\t", **kwargs)


def _coerce_float(value):
    """Return a finite float when possible, otherwise None."""
    if value in (None, "", "N/A", "NA", "nan", "NaN"):
        return None
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(numeric):
        return None
    return numeric


def _format_count_value(value):
    """Format count-like values for the report, tolerating NA-like strings."""
    numeric = _coerce_float(value)
    if numeric is None:
        return "N/A"
    return format(round(numeric), ",")


def _format_ratio_value(value):
    """Format ratio values for the report, tolerating NA-like strings."""
    numeric = _coerce_float(value)
    if numeric is None:
        return "N/A"
    return f"{numeric:.2f}x"


def _sorted_transcript_abundance_table(tx_counts_file, sample_aliases):
    """Get transcript abundance rows by descending total abundance of sample columns."""
    tx_counts = _read_table(tx_counts_file)
    if tx_counts is None or tx_counts.empty:
        return tx_counts

    count_columns = [
        column for column in tx_counts.columns if column in set(sample_aliases)
    ]
    if not count_columns:
        return tx_counts

    numeric_counts = tx_counts[count_columns].apply(pd.to_numeric, errors="coerce")
    abundance = numeric_counts.sum(axis=1, min_count=1)
    order = abundance.fillna(float("-inf")).sort_values(ascending=False).index
    return tx_counts.loc[order].reset_index(drop=True)


def _format_classification_label(name):
    """Return canonical report label for a summary classification value."""
    return str(name).strip().replace("-", "_").replace("_", " ").capitalize()


def _transcriptome_summary(transcriptome_dir):
    """Return transcriptome metrics and transcript class counts DataFrames."""
    tx_meta = _read_table(Path(transcriptome_dir) / "transcript_metadata.tsv")
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
    """Return a dict of per-sample metrics DataFrames keyed by sample name."""
    summaries = {}
    samples_path = Path(samples_dir)
    if not samples_path.exists() or not samples_path.is_dir():
        return summaries
    for sample_dir in sorted(samples_path.iterdir()):
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


def _sqanti_table(sqanti_dir):
    """Return a SQANTI3 classification summary DataFrame."""
    rows = []
    summaries = []
    for root, _, files in os.walk(sqanti_dir, followlinks=True):
        if "classification_summary.tsv" in files:
            summaries.append(Path(root) / "classification_summary.tsv")
    for summary in sorted(summaries):
        table = _read_table(summary)
        if table is None or table.empty:
            continue

        sample = summary.parent.name

        sample_counts = {"Sample": sample}
        for _, row in table.iterrows():
            feature = _format_classification_label(row["structural_category"])
            count = _coerce_float(row["count"])
            sample_counts[feature] = int(round(count)) if count is not None else 0
        rows.append(sample_counts)

    if not rows:
        return None

    sqanti_df = pd.DataFrame(rows).fillna(0)
    feature_cols = list(classification_categories.keys())
    for col in feature_cols:
        if col not in sqanti_df.columns:
            sqanti_df[col] = 0
        sqanti_df[col] = sqanti_df[col].astype(int)
    sqanti_df = sqanti_df[["Sample"] + feature_cols]
    is_cohort = sqanti_df["Sample"].str.lower().eq("cohort")
    sqanti_df = sqanti_df.assign(_is_cohort=is_cohort)
    return sqanti_df.sort_values(["_is_cohort", "Sample"]).drop(columns="_is_cohort")


def _sample_mod_summaries(summary_dir):
    """Return a combined modified base summary table across all samples."""
    if not summary_dir:
        return None
    summaries_path = Path(summary_dir)
    if not summaries_path.exists() or not summaries_path.is_dir():
        return None
    summaries = []
    for summary_file in sorted(summaries_path.glob("*.mods.summary.tsv")):
        if not summary_file.is_file():
            continue
        summary = _read_table(summary_file)
        if summary is None or summary.empty:
            continue
        if "sample" not in summary.columns:
            raise ValueError(
                f"Modified base summary is missing required 'sample' column: "
                f"{summary_file}"
            )

        sample_names = summary["sample"].unique()
        if len(sample_names) != 1:
            raise ValueError(
                f"Modified base summary must contain exactly one sample name: "
                f"{summary_file}"
            )

        summary["mod"] = summary["mod_label"].where(
            summary["mod_label"] != "",
            summary["full_mod_code"],
        )
        summary = summary[
            [
                "sample",
                "mod",
                "full_mod_code",
                "modification_percent",
                "modified_calls",
                "valid_coverage",
            ]
        ]
        summaries.append(summary)

    if not summaries:
        return None

    summary = pd.concat(summaries, ignore_index=True)
    return summary.sort_values(["mod", "sample"]).reset_index(drop=True)


def _render_mod_summary_matrix(summary):
    """Render modified base summaries as a cross-sample comparison matrix."""
    if summary.duplicated(["sample", "mod"]).any():
        raise ValueError(
            "Modified base summary contains duplicate sample/mod combinations."
        )

    samples = sorted(summary["sample"].astype(str).unique())
    mods = sorted(summary["mod"].astype(str).unique())
    lookup = {
        (str(row.sample), str(row.mod)): row
        for row in summary.itertuples(index=False)
    }

    header_html = "".join(
        f"<th>{escape(mod)}</th>" for mod in mods
    )
    body_rows = []
    for sample in samples:
        cells = []
        for mod in mods:
            row = lookup.get((sample, mod))
            if row is None:
                cells.append(
                    "<td class='mod-summary-cell mod-summary-empty'>"
                    "<div class='mod-summary-empty-mark'>&mdash;</div>"
                    "</td>"
                )
                continue
            cells.append(
                "<td class='mod-summary-cell'>"
                f"<div class='mod-summary-percent'>"
                f"{_format_mod_summary_percent(row.modification_percent)}%"
                f"</div>"
                f"<div class='mod-summary-meta'>"
                f"{_format_mod_summary_count(row.modified_calls)} modified"
                f"</div>"
                f"<div class='mod-summary-meta'>"
                f"{_format_mod_summary_count(row.valid_coverage)} valid"
                f"</div>"
                "</td>"
            )
        body_rows.append(
            "<tr>"
            f"<th class='mod-summary-sample'>{escape(sample)}</th>"
            + "".join(cells)
            + "</tr>"
        )

    return (
        "<div class='mod-summary-matrix-wrap'>"
        "<table class='mod-summary-matrix'>"
        "<thead><tr><th>Sample</th>"
        f"{header_html}</tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table>"
        "</div>"
    )


def _format_mod_summary_percent(value):
    """Format a modified base percentage for display."""
    return f"{float(value):.2f}"


def _format_mod_summary_count(value):
    """Format a modified base count for display."""
    return format(int(round(float(value))), ",")


def _mod_summary_matrix_style():
    """Return CSS for the modified base summary comparison matrix."""
    return """
    .mod-summary-matrix-wrap { overflow-x: auto; margin-bottom: 0.75rem; }
    .mod-summary-matrix { width: 100%; border-collapse: collapse; }
    .mod-summary-matrix th,
    .mod-summary-matrix td {
        border-bottom: 1px solid #e5e7eb;
        padding: 0.75rem 0.875rem;
        vertical-align: top;
        text-align: left;
    }
    .mod-summary-matrix thead th {
        font-weight: 600;
        white-space: nowrap;
    }
    .mod-summary-sample {
        white-space: nowrap;
        font-weight: 600;
    }
    .mod-summary-percent {
        font-size: 1.05rem;
        font-weight: 700;
        line-height: 1.2;
    }
    .mod-summary-meta {
        margin-top: 0.15rem;
        font-size: 0.8rem;
        color: #6b7280;
        line-height: 1.25;
        white-space: nowrap;
    }
    .mod-summary-empty {
        color: #9ca3af;
    }
    .mod-summary-empty-mark {
        font-size: 1rem;
        line-height: 1.2;
    }
    """


def _contrast_results(de_dir, filename):
    """Return a dict of per-contrast result DataFrames read from filename."""
    tables = {}
    for contrast_dir in sorted(Path(de_dir).iterdir()):
        if not contrast_dir.is_dir():
            continue
        # Enforce str dtype in case of all Nan values.
        table = _read_table(
            contrast_dir / filename,
            dtype={'gene_name': str, 'transcript_name': str}
        )
        if table is None or table.empty:
            continue

        if 'gene_name' in table.columns:
            table["gene_name"] = table["gene_name"].fillna("-")
        if 'transcript_name' in table.columns:
            table["transcript_name"] = table["transcript_name"].fillna("-")

        if "padj" in table.columns:
            table.sort_values("padj", ascending=True, inplace=True)
        tables[contrast_dir.name] = table

    return tables


def _round_de_table(table):
    """Return a display-formatted copy of a DE/DTU result table."""
    rounded_columns = [
        "baseMean", "exonBaseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    ]

    def _format_numeric(value):
        numeric_value = pd.to_numeric(value, errors="coerce")
        if pd.isna(numeric_value):
            return value
        fixed_decimal = f"{numeric_value:.3f}"
        if numeric_value != 0 and (
            abs(numeric_value) < 0.0001 or fixed_decimal in {"0.000", "-0.000"}
        ):
            return f"{numeric_value:.3e}"
        return fixed_decimal

    rounded = table.copy()
    for column in rounded_columns:
        if column not in rounded.columns:
            continue
        rounded[column] = rounded[column].map(_format_numeric)
    return rounded


def _load_bambu_qc(bambu_dir):
    """Load bambu QC statistics JSON."""
    qc_file = Path(bambu_dir) / "bambu_qc_stats.json"
    if qc_file.exists():
        with open(qc_file) as f:
            return json.load(f)
    return None


def _load_annotation_reference_summary(summary_file):
    """Load reference/annotation preparation summary JSON from a file path."""
    if summary_file is None:
        return None
    summary_path = Path(summary_file)
    if summary_path.exists():
        with open(summary_path) as f:
            return json.load(f)
    return None


def _load_de_qc(de_dir):
    """Load DE/DTU QC statistics JSON."""
    qc_file = Path(de_dir) / "de_qc_stats.json"
    if qc_file.exists():
        with open(qc_file) as f:
            return json.load(f)
    return None


def _load_cpm_tables(cohort_dir):
    """Load cohort-level gene and transcript CPM tables."""
    cohort_dir = Path(cohort_dir)
    gene_cpm = _read_table(cohort_dir / "gene_cpm.tsv")

    if gene_cpm is None or gene_cpm.empty:
        gene_cpm = None

    transcript_cpm = _read_table(cohort_dir / "transcript_cpm.tsv")
    if transcript_cpm is None or transcript_cpm.empty:
        transcript_cpm = None

    return {
        "gene": gene_cpm,
        "transcript": transcript_cpm,
    }


def _load_cohort_samples(cohort_dir):
    """Load cohort sample metadata CSV."""
    sample_file = Path(cohort_dir) / "samples.csv"
    if not sample_file.exists():
        return None
    return pd.read_csv(sample_file)


def _format_hint_values(hints):
    """Format provenance hints for a compact table cell."""
    if not hints:
        return "None detected"
    return ", ".join(hints)


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


def _heatmap_style():
    return """
        .heatmap-table-grid {
            display: grid;
            grid-template-columns: repeat(3, minmax(0, 1fr));
            gap: 20px 10px;
            align-items: start;
        }
        .heatmap-table-grid > * {
            min-width: 0;
        }
        @media screen and (max-width: 1000px) {
            .heatmap-table-grid {
                grid-template-columns: 1fr;
            }
        }
        .clustering-info {
            font-size: 11px;
        }"""


def _volcano_style():
    return """
        .volcano-plot-grid {
            display: grid;
            grid-template-columns: minmax(0, 1fr) minmax(430px, 34%);
            gap: 18px;
            align-items: start;
        }
        .volcano-plot-grid > *,
        .volcano-side-panel > * {
            min-width: 0;
        }
        .volcano-side-panel {
            display: grid;
            gap: 10px;
            align-items: start;
        }
        @media screen and (max-width: 1100px) {
            .volcano-plot-grid {
                grid-template-columns: 1fr;
            }
        }
        """


def _as_string_list(value):
    """Normalize optional values to a compact list of strings."""
    if value is None or value == "none":
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value if item not in (None, "")]
    if isinstance(value, str):
        return [value] if value else []
    return [str(value)]


def _collect_de_method_rows(de_qc):
    """Build per-contrast method rows and warning metadata."""
    rows = []
    deseq2_gene_wise = []
    dexseq_gene_wise = []
    dexseq_covariate_drops = []
    failed_dge = []

    for contrast_name, contrast_data in de_qc.get("contrasts", {}).items():
        fallback = contrast_data.get("deseq2_dispersion_fallback") or {}
        fallback_applied = bool(fallback.get("applied", False))
        deseq2_method = fallback.get("method_used")
        deseq2_size_factors = contrast_data.get("deseq2_size_factor_method") or "ratio"

        if not deseq2_method:
            deseq2_method = "gene-wise" if fallback_applied else "parametric"

        if deseq2_method == "gene-wise":
            deseq2_gene_wise.append(contrast_name)

        dexseq_method = contrast_data.get("dexseq_dispersion_method") or "parametric"
        dexseq_size_factors = contrast_data.get("dexseq_size_factor_method") or "ratio"
        if dexseq_method == "gene-wise":
            dexseq_gene_wise.append(contrast_name)

        dropped_covariates = _as_string_list(
            contrast_data.get("dexseq_covariates_dropped")
        )
        if dropped_covariates:
            dexseq_covariate_drops.append((contrast_name, dropped_covariates))

        if contrast_data.get("dge_status") == "FAILED":
            failed_dge.append(contrast_name)

        rows.append(
            {
                "Contrast": contrast_name,
                "DESeq2 size factors": deseq2_size_factors,
                "DESeq2 dispersion": (
                    f"{deseq2_method} (fallback)"
                    if fallback_applied
                    else deseq2_method
                ),
                "DEXSeq size factors": dexseq_size_factors,
                "DEXSeq dispersion": dexseq_method,
                "DEXSeq covariates dropped": (
                    ", ".join(dropped_covariates) if dropped_covariates else "none"
                ),
                "DGE status": contrast_data.get("dge_status", "N/A"),
                "DTU status": contrast_data.get("dtu_status", "N/A"),
            }
        )

    return rows, deseq2_gene_wise, dexseq_gene_wise, dexseq_covariate_drops, failed_dge


def main(args):
    """Run the report entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "Transcriptomes Sequencing report",
        "wf-transcriptomes",
        args.params,
        args.versions,
        args.wf_version,
        head_resources=[
            *LAB_head_resources,
            EZC_Resource(func=get_bokeh_widgets_js, tag=script),
            EZC_Resource(func=get_bokeh_tables_js, tag=script)]
    )

    with open(args.metadata, "r") as handle:
        metadata = json.load(handle)

    if args.stats:
        with report.add_section("Read summary", "Reads"):
            p(
                "Read quality and alignment statistics generated by fastcat bamstats "
                "for each input sample. Shown are distributions of read "
                "length and quality score, and a summary of alignment "
                "outcomes (mapped, unmapped, primary, supplementary) "
                "against the reference genome."
            )
            stats = tuple(args.stats)
            sample_names = tuple(
                item["alias"] for item in metadata if item.get("has_stats")
            )
            flagstats = tuple(
                Path(stats_dir) / "bamstats.flagstat.tsv" for stats_dir in stats
            )
            if len(stats) == 1:
                stats = stats[0]
                flagstats = flagstats[0]
                sample_names = sample_names[0] if sample_names else None
            try:
                fastcat.SeqSummary(
                    stats,
                    flagstat=flagstats,
                    sample_names=sample_names,
                    alignment_stats=True,
                )
            except Exception as exc:  # pragma: no cover - defensive
                logger.warning("Skipping read summary plot: %s", exc)
                _create_warning_banner(
                    (
                        "Read summary plots could not be rendered for this run. "
                        "This can happen for degenerate or extremely "
                        "small input statistics."
                    ),
                    level="info",
                )

    with report.add_section("Sample metadata", "Samples"):
        p(
            "Metadata for each input sample as parsed from the sample sheet "
            "and workflow parameters. Includes the sample alias, barcode, "
            "type, and any experimental design columns such as condition or "
            "batch that were supplied for differential analysis."
        )
        tabs = Tabs()
        for item in sorted(metadata, key=lambda value: value["alias"]):
            with tabs.add_tab(item["alias"]):
                DataTable.from_pandas(
                    pd.DataFrame.from_dict(item, orient="index", columns=["Value"])
                    .reset_index()
                    .rename(columns={"index": "Field"}),
                    use_index=False,
                )
    mod_summaries = _sample_mod_summaries(args.mod_summary_dir)
    if mod_summaries is not None and not mod_summaries.empty:
        with report.add_section("Modified base summaries", "Modifications"):
            p(
                "Summary of base modification calls from modkit pileup for each "
                "combination of sample and modification type."
            )
            dom_style(raw(_mod_summary_matrix_style()))
            raw(_render_mod_summary_matrix(mod_summaries))
            small(raw(
                "<b>Valid coverage</b> Sum of the valid residues: all the "
                "modified, canonical and other mod (where the modification "
                "is different from the listed base) bedMethyl columns "
                "counts for this combination of sample and mod. "
                "<b>Modified calls</b> Number of calls passing filters "
                "that were classified as a residue with a specified base "
                "modification. "
                "<b>Modification percent</b> (Modified calls / Valid "
                "coverage) * 100"
            ))

    annotation_reference_summary = _load_annotation_reference_summary(args.ref_summary)

    if annotation_reference_summary:
        with report.add_section("Reference and Annotation Checks", "Reference"):
            p(
                "Results of compatibility checks between the supplied reference "
                "genome and annotation performed before bambu transcript "
                "modelling. "
                "Build and provider hints are inferred from sequence names and "
                "file content to help identify mismatched genome/annotation "
                "combinations."
            )
            for warning in annotation_reference_summary.get("warnings", []):
                _create_warning_banner(warning, level="warning")

            seqname_rows = [
                (
                    "Overlapping seqnames",
                    len(annotation_reference_summary.get("seqname_overlap", [])),
                ),
                (
                    "Seqnames only in annotation",
                    len(annotation_reference_summary.get("only_in_annotation", [])),
                ),
                (
                    "Seqnames only in reference",
                    len(annotation_reference_summary.get("only_in_reference", [])),
                ),
            ]
            annotation_summary = annotation_reference_summary.get("annotation", {})
            seqname_rows.extend([
                (
                    "Annotation records retained",
                    annotation_summary.get("kept_records", "N/A"),
                ),
                (
                    "Unstranded records excluded",
                    annotation_summary.get("excluded_unstranded_records", "N/A"),
                ),
                (
                    "Annotation attributes sanitised",
                    annotation_summary.get("sanitised_attribute_records", "N/A"),
                ),
            ])
            DataTable.from_pandas(
                pd.DataFrame(seqname_rows, columns=["Check", "Value"]),
                paging=False,
                searchable=False,
                use_index=False,
            )

            h4("Build and Provider Hints")
            hint_rows = [
                (
                    "Reference build",
                    _format_hint_values(
                        annotation_reference_summary.get(
                            "reference_build_hints", []
                        )
                    ),
                ),
                (
                    "Annotation build",
                    _format_hint_values(
                        annotation_reference_summary.get(
                            "annotation_build_hints", []
                        )
                    ),
                ),
                (
                    "Reference provider",
                    _format_hint_values(
                        annotation_reference_summary.get(
                            "reference_provider_hints", []
                        )
                    ),
                ),
                (
                    "Annotation provider",
                    _format_hint_values(
                        annotation_reference_summary.get(
                            "annotation_provider_hints", []
                        )
                    ),
                ),
            ]
            DataTable.from_pandas(
                pd.DataFrame(hint_rows, columns=["Evidence", "Hints"]),
                paging=False,
                searchable=False,
                use_index=False,
            )

            examples = annotation_summary.get("unstranded_examples") or []
            if examples:
                h4("Unstranded Annotation Examples")
                pre("\n".join(examples))

    # Setup for using cohort or single sample bambu results
    is_single_sample = len(metadata) == 1
    primary_label = metadata[0]["alias"] if is_single_sample else "Cohort"

    # Load bambu QC statistics
    bambu_dir = (
        Path(args.samples_dir) / metadata[0]["alias"] if is_single_sample
        else args.cohort_dir
    )
    bambu_qc = _load_bambu_qc(bambu_dir)
    # Add Bambu QC section with warnings
    if bambu_qc:
        with report.add_section("Bambu Quality Control", "Bambu QC"):
            p(
                "Quality metrics from the bambu transcript discovery and "
                "quantification run. Library size statistics show cohort-level "
                "summary values to "
                "flag large inter-sample variation that could affect CPM "
                "normalisation, with the individual per-sample counts listed "
                "below. The transcript discovery table shows the transcriptome "
                "mode used, the novel discovery rate (NDR) threshold applied, "
                "and how many transcripts were present before and after "
                "low-count filtering."
            )
            # Check for warnings
            if bambu_qc.get("library_size_warning"):
                _create_warning_banner(
                    f"Library Size Variation: {bambu_qc['library_size_warning']}. "
                    "Large variation (>3x) may affect CPM normalization. "
                    "Consider reviewing per-sample library sizes.",
                    level="warning",
                )

            h4("Library Size Statistics")
            lib_stats = pd.DataFrame(
                [
                    ("Samples analyzed", bambu_qc.get("samples", "N/A")),
                    (
                        "Median library size",
                        "{} reads".format(
                            _format_count_value(
                                bambu_qc.get("median_library_size", 0)
                            )
                        ),
                    ),
                    (
                        "Min library size",
                        "{} reads".format(
                            _format_count_value(
                                bambu_qc.get("min_library_size", 0)
                            )
                        ),
                    ),
                    (
                        "Max library size",
                        "{} reads".format(
                            _format_count_value(
                                bambu_qc.get("max_library_size", 0)
                            )
                        ),
                    ),
                    (
                        "Library size ratio (max/min)",
                        _format_ratio_value(
                            bambu_qc.get("library_size_ratio", 1.0)
                        ),
                    ),
                ],
                columns=["Metric", "Value"],
            )
            DataTable.from_pandas(
                lib_stats,
                paging=False,
                searchable=False,
                use_index=False,
            )

            # Transcript discovery statistics
            h4("Transcript Discovery")
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
            DataTable.from_pandas(
                discovery_stats,
                paging=False,
                searchable=False,
                use_index=False,
            )

            # Per-sample library sizes
            if "library_sizes" in bambu_qc and bambu_qc["library_sizes"]:
                h4("Per-Sample Library Sizes")
                lib_size_data = []
                for sample, size in bambu_qc["library_sizes"].items():
                    numeric_size = _coerce_float(size)
                    lib_size_data.append(
                        {
                            "Sample": sample,
                            "Library Size": _format_count_value(size),
                            "Reads": (
                                numeric_size if numeric_size is not None else -1
                            ),
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

    with report.add_section(
        f"{primary_label} transcriptome",
        f"{primary_label} transcriptome"
    ):
        p(
            f"Summary of the "
            f"{'per-sample' if is_single_sample else 'joint cohort'} bambu "
            "transcriptome model. Shows the total number of transcripts and "
            "genes in the final model, and "
            f"the top 500 rows of the transcript abundance count table."
        )
        transcriptome_metrics, transcriptome_classes = _transcriptome_summary(bambu_dir)
        if transcriptome_metrics is not None:
            DataTable.from_pandas(
                transcriptome_metrics,
                paging=False,
                searchable=False,
                use_index=False,
            )
        if transcriptome_classes is not None:
            DataTable.from_pandas(
                transcriptome_classes,
                paging=False,
                searchable=False,
                use_index=False,
            )

        tx_counts = _sorted_transcript_abundance_table(
            Path(bambu_dir) / "transcript_counts.tsv",
            [item["alias"] for item in metadata]
        )
        if tx_counts is not None and not tx_counts.empty:
            p(
                f"Top {args.de_table_size} most abundant transcripts"
            )
            DataTable.from_pandas(tx_counts.head(args.de_table_size), use_index=False)
        else:
            _create_warning_banner("No trancrips discovered")

    if not is_single_sample:
        with report.add_section(
            "Per-sample transcriptomes",
            "Per-sample transcriptomes"
        ):
            p(
                "Per-sample transcript and gene count summaries derived from "
                "the individual bambu quantification runs. Each tab shows the "
                "number of transcripts and genes detected in that sample after "
                "filtering. These can be used to spot samples with unusually "
                "low transcript detection compared to the rest of the cohort."
            )
            tabs = Tabs()
            for sample, summary_df in _sample_summaries(args.samples_dir).items():
                with tabs.add_tab(sample):
                    DataTable.from_pandas(
                        summary_df,
                        paging=False,
                        searchable=False,
                        use_index=False,
                    )

    sqanti_table = _sqanti_table(args.sqanti_dir)
    if sqanti_table is not None and not sqanti_table.empty:
        with report.add_section("SQANTI3 classification", "SQANTI3"):
            p(
                "Structural classification of transcript isoforms by SQANTI3. "
                "Each transcript is assigned a category based on how its "
                "splice junctions and exon structure compared to the reference "
                "annotation. The table shows the count of "
                "transcripts in each category per sample and for "
                "the whole cohort."
            )
            DataTable.from_pandas(sqanti_table, use_index=False)
            with p():
                for category, description in classification_categories.items():
                    small(strong(f"{category}: "))
                    small(raw(f"{description}<br>"))

    if args.de_dir and Path(args.de_dir).exists():
        # Load DE QC statistics
        de_qc = _load_de_qc(args.de_dir)

        # Add DE/DTU QC section with warnings
        if de_qc:
            with report.add_section(
                "Differential Analysis Quality Control",
                "DE/DTU QC",
            ):
                p(
                    "Quality control summary for the differential expression "
                    "analyses. Shows the experimental design "
                    "and the number of samples per group. For each contrast, "
                    "the statistical methods chosen by DESeq2 and DEXSeq are "
                    "reported , including the dispersion estimation strategy "
                    "(parametric or gene-wise fallback) and size factor "
                    "normalisation method, alongside the analysis status and "
                    "count of significant hits. Any warnings about low sample "
                    "numbers, dispersion fallbacks, dropped covariates, or "
                    "analysis failures are highlighted here."
                )
                # Check for critical warnings
                has_warnings = False
                sample_size_warnings = _as_string_list(
                    de_qc.get("sample_size_warnings")
                )
                (
                    method_rows,
                    deseq2_gene_wise,
                    dexseq_gene_wise,
                    dexseq_covariate_drops,
                    failed_dge,
                ) = _collect_de_method_rows(de_qc)

                if sample_size_warnings:
                    _create_warning_banner(
                        "Sample Size Warning: "
                        + "; ".join(sample_size_warnings)
                        + ". "
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

                if deseq2_gene_wise or dexseq_gene_wise:
                    gene_wise_details = []
                    if deseq2_gene_wise:
                        gene_wise_details.append(
                            "DESeq2: " + ", ".join(sorted(deseq2_gene_wise))
                        )
                    if dexseq_gene_wise:
                        gene_wise_details.append(
                            "DEXSeq: " + ", ".join(sorted(dexseq_gene_wise))
                        )
                    _create_warning_banner(
                        "Gene-wise dispersion fallback used (reduced power). "
                        + " ".join(gene_wise_details),
                        level="warning",
                    )
                    has_warnings = True

                if dexseq_covariate_drops:
                    drop_details = [
                        f"{contrast} ({', '.join(columns)})"
                        for contrast, columns in sorted(dexseq_covariate_drops)
                    ]
                    _create_warning_banner(
                        "DEXSeq covariates dropped due to rank-deficient design. "
                        f"Affected: {'; '.join(drop_details)}",
                        level="warning",
                    )
                    has_warnings = True

                if failed_dge:
                    _create_warning_banner(
                        "DGE Analysis Failed: "
                        f"{len(failed_dge)} contrast(s) could not complete "
                        "DGE testing. "
                        f"Affected: {', '.join(failed_dge)}. "
                        "See DGE_ANALYSIS_FAILED.txt files for details.",
                        level="danger",
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
                h4("Experimental Design")
                covariates = _as_string_list(de_qc.get("covariates"))
                covariates_value = ", ".join(covariates) if covariates else "none"
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
                DataTable.from_pandas(
                    design_stats,
                    paging=False,
                    searchable=False,
                    use_index=False,
                )

                # Sample sizes per group
                if "samples_per_group" in de_qc:
                    h4("Sample Sizes per Group")
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

                h4("Statistical Methods & Warnings")
                if method_rows:
                    method_df = pd.DataFrame(method_rows)
                    DataTable.from_pandas(
                        method_df,
                        paging=False,
                        use_index=False,
                    )
                else:
                    p("No contrast-level QC metadata was found.")

                # Per-contrast summary
                if "contrasts" in de_qc:
                    h4("Results Summary by Contrast")
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
                                "DGE Status": contrast_data.get(
                                    "dge_status", "N/A"
                                ),
                                "DGE Significant (FDR<0.05)": (
                                    contrast_data.get(
                                        "dge_significant_fdr05", 0
                                    )
                                    if contrast_data.get("dge_status") == "SUCCESS"
                                    else "N/A"
                                ),
                                "DGE Up": (
                                    contrast_data.get("dge_upregulated", 0)
                                    if contrast_data.get("dge_status") == "SUCCESS"
                                    else "N/A"
                                ),
                                "DGE Down": (
                                    contrast_data.get("dge_downregulated", 0)
                                    if contrast_data.get("dge_status") == "SUCCESS"
                                    else "N/A"
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
                    h4("Quality Warnings Summary")
                    warnings_data = []
                    if sample_size_warnings:
                        warnings_data.append(
                            {
                                "Warning Type": "Sample Size",
                                "Details": "; ".join(sample_size_warnings),
                            }
                        )
                    if deseq2_gene_wise or dexseq_gene_wise:
                        engines = []
                        if deseq2_gene_wise:
                            engines.append(
                                f"DESeq2 ({len(deseq2_gene_wise)} contrasts)"
                            )
                        if dexseq_gene_wise:
                            engines.append(
                                f"DEXSeq ({len(dexseq_gene_wise)} contrasts)"
                            )
                        warnings_data.append(
                            {
                                "Warning Type": "Gene-wise Dispersion Fallback",
                                "Details": "; ".join(engines),
                            }
                        )
                    if dexseq_covariate_drops:
                        warnings_data.append(
                            {
                                "Warning Type": "DEXSeq Covariates Dropped",
                                "Details": (
                                    f"{len(dexseq_covariate_drops)} "
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
                    if failed_dge:
                        warnings_data.append(
                            {
                                "Warning Type": "DGE Failure",
                                "Details": (
                                    f"{len(failed_dge)} contrasts failed"
                                ),
                            }
                        )

                    warnings_df = pd.DataFrame(warnings_data)
                    DataTable.from_pandas(warnings_df, paging=False, use_index=False)

        if de_qc:
            condition_column = de_qc.get("condition_column")
            cohort_cpm = _load_cpm_tables(args.cohort_dir)
            cohort_samples = _load_cohort_samples(args.cohort_dir)

        with report.add_section("Differential gene expression", "DGE"):
            p(
                "Differential gene expression results from DESeq2. The "
                "heatmap, PCA plot, and sample distance matrix are derived "
                "from CPM-normalised counts "
                "and give an overview of sample clustering relative to "
                "the experimental conditions. Each contrast tab shows a "
                "results table of genes ranked by adjusted p-value and a "
                "volcano plot highlighting significantly up- and "
                "down-regulated genes."
            )
            dom_style(raw(_heatmap_style() + _volcano_style()))
            if condition_column:
                if cohort_cpm['gene'] is None:
                    _create_warning_banner(
                        "Cohort gene CPM table is missing or empty. ")
                else:
                    hierarchical_result = hierarchical(
                        cohort_cpm["gene"],
                        id_column="GENEID",
                        samples=cohort_samples,
                        condition_column=condition_column,
                        top_n=150,
                    )
                    if hierarchical_result.error is not None:
                        _create_warning_banner(
                            hierarchical_result.error, level='warning')
                    else:
                        with div(cls="heatmap-table-grid"):
                            EZChart(hierarchical_result.heatmap, width="100%")
                            EZChart(hierarchical_result.pca, width="100%")
                            EZChart(hierarchical_result.distance, width="100%")
                        with div(cls="clustering-info"):
                            br()
                            clustering_info('gene')
            tabs = Tabs()
            for contrast, table in _contrast_results(
                args.de_dir, "results_dge.tsv"
            ).items():
                with tabs.add_tab(contrast):
                    # Check for contrast-specific warnings
                    if de_qc and contrast in de_qc.get("contrasts", {}):
                        contrast_data = de_qc["contrasts"][contrast]
                        if contrast_data.get("dge_status") == "FAILED":
                            _create_warning_banner(
                                "DGE analysis failed for this contrast. "
                                f"See {contrast}/"
                                "DGE_ANALYSIS_FAILED.txt for detailed "
                                "explanation.",
                                level="danger",
                            )
                            p(
                                "Empty results indicate analysis failure, "
                                "not 'no DGE detected'."
                            )
                        elif contrast_data.get("dtu_power_warning"):
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
                    DataTable.from_pandas(
                        _round_de_table(table.head(args.de_table_size)),
                        use_index=False
                    )
                    with div(cls="clustering-info"):
                        raw(
                            f"Table showing the top {args.de_table_size} genes sorted "
                            "by adjusted p-value. <br><br><br>"
                        )

                    h3("Gene expression volcano Plot")
                    gn_vol, gn_class_table, gn_selected_table = volcano(table)
                    with div(_class="volcano-plot-grid"):
                        EZChart(gn_vol, width="100%", height="550")
                        with div(_class="volcano-side-panel"):
                            EZChart(gn_class_table, width="100%", height="auto")
                            EZChart(gn_selected_table, width="100%", height="auto")

        with report.add_section("Differential transcript usage", "DTU"):
            p(
                "Differential transcript usage results from DEXSeq. Unlike "
                "DGE, DTU tests whether individual transcripts change their "
                "proportional contribution to total gene expression between "
                "conditions , i.e. isoform switching , rather than testing "
                "for changes in total gene abundance. The heatmap, PCA, and "
                "sample distance matrix use CPM-normalised counts. "
                "Each contrast tab shows a "
                "results table ranked by adjusted p-value and a volcano plot."
            )
            if condition_column:
                if cohort_cpm["transcript"] is None:
                    _create_warning_banner(
                        "Cohort transcript CPM table is missing or empty.")
                else:
                    hierarchical_result = hierarchical(
                        cohort_cpm["transcript"],
                        id_column="TXNAME",
                        top_n=150,
                        samples=cohort_samples,
                        condition_column=condition_column
                    )
                    if hierarchical_result.error is not None:
                        _create_warning_banner(
                            hierarchical_result.error, level='warning')
                    else:
                        with div(cls="heatmap-table-grid"):
                            EZChart(hierarchical_result.heatmap, width="100%")
                            EZChart(hierarchical_result.pca, width="100%")
                            EZChart(hierarchical_result.distance, width="100%")
                        with div(cls="clustering-info"):
                            br()
                            clustering_info('transcript')

            tabs = Tabs()
            dtu_tables = _contrast_results(
                args.de_dir, "results_dtu_transcript.tsv")

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
                            failure_hint = contrast_data.get("dtu_failure_hint")
                            if failure_hint:
                                p(f"Probable cause: {failure_hint}")
                            p(
                                "Empty results indicate analysis failure, "
                                "not 'no DTU detected'."
                            )
                        elif contrast_data.get("dtu_power_warning"):
                            _create_warning_banner(
                                contrast_data["dtu_power_warning"],
                                level="warning",
                            )

                    if contrast_name in dtu_tables:
                        dtu_table = dtu_tables[contrast_name]
                        DataTable.from_pandas(
                            _round_de_table(dtu_table.head(args.de_table_size)),
                            use_index=False
                        )
                        with div(cls="clustering-info"):
                            raw(
                                f"Table showing the top {args.de_table_size} "
                                "transcripts sorted by adjusted p-value. <br><br><br>"
                            )

                        h3("Transcript expression volcano Plot")
                        tr_vol, tr_class_table, tr_selected_table = volcano(dtu_table)
                        with div(_class="volcano-plot-grid"):
                            EZChart(tr_vol, width="100%", height="550")
                            with div(_class="volcano-side-panel"):
                                EZChart(tr_class_table, width="100%", height="auto")
                                EZChart(tr_selected_table, width="100%", height="auto")
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
        "--mod_summary_dir",
        default=None,
        help="Modified base summary directory.",
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
    parser.add_argument(
        "--ref_summary",
        default=None,
        help="Annotation reference summary TSV.",
    )
    parser.add_argument(
        "--de_table_size",
        default=500,
        type=int,
        help="Number of rows to show in DE/DTU result tables.",
    )
    parser.add_argument("--versions", required=True, help="Versions directory.")
    parser.add_argument("--params", required=True, help="Workflow params JSON.")
    parser.add_argument(
        "--wf_version",
        default="unknown",
        help="Workflow version.",
    )
    return parser
