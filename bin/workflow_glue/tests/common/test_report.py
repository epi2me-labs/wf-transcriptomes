"""Tests for the workflow report entry point."""

import json
from pathlib import Path

import pandas as pd
from workflow_glue import report


class _NullContext:
    """Minimal context manager used by report section and tab stubs."""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False


class _FakeReport:
    """Small stand-in for the ezcharts report wrapper."""

    def __init__(self, *args, **kwargs):
        self.sections = []

    def add_section(self, title, key):
        self.sections.append((title, key))
        return _NullContext()

    def write(self, path):
        Path(path).write_text("report ok\n", encoding="utf-8")


class _FakeTabs:
    """Small stand-in for the tab layout helper."""

    def add_tab(self, label):
        return _NullContext()


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def _build_report_args(tmp_path, de_qc=None):
    """Create minimal report inputs, optionally including DE QC JSON."""
    metadata = _write(
        tmp_path / "metadata.json",
        json.dumps([{"alias": "sampleA", "has_stats": False}]),
    )
    params = _write(tmp_path / "params.json", "{}")
    versions = tmp_path / "versions"
    versions.mkdir()
    _write(versions / "versions.txt", "tool,1.0\n")

    cohort = tmp_path / "cohort"
    cohort.mkdir()
    reference = cohort / "reference"
    reference.mkdir()
    _write(
        reference / "annotation_reference_summary.json",
        json.dumps(
            {
                "seqname_overlap": ["chr1"],
                "only_in_annotation": [],
                "only_in_reference": [],
                "annotation": {
                    "kept_records": 10,
                    "excluded_unstranded_records": 0,
                    "sanitised_attribute_records": 0,
                },
                "warnings": [],
            }
        ),
    )

    if de_qc is not None and "contrasts" in de_qc:
        # Create samples.csv with multisample metadata for heatmaps
        samples_csv = "barcode,sample_id,alias,condition\n"
        for i in range(3):
            samples_csv += (
                f"BC{i:03d},sample_control_{i},sample_control_{i},control\n"
            )
        for i in range(3):
            samples_csv += (
                f"BC{i+3:03d},sample_treated_{i},sample_treated_{i},treated\n"
            )
        _write(cohort / "samples.csv", samples_csv)

        # Create minimal CPM tables for hierarchical clustering
        sample_cols = (
            "\tsample_control_0\tsample_control_1\tsample_control_2"
            "\tsample_treated_0\tsample_treated_1\tsample_treated_2\n"
        )
        gene_cpm = f"GENEID{sample_cols}"
        gene_cpm += "gene1\t100\t110\t95\t200\t220\t210\n"
        gene_cpm += "gene2\t50\t55\t48\t100\t110\t105\n"
        gene_cpm += "gene3\t75\t80\t72\t150\t160\t155\n"
        _write(cohort / "gene_cpm.tsv", gene_cpm)

        tx_cpm = f"TXNAME{sample_cols}"
        tx_cpm += "tx1\t100\t110\t95\t200\t220\t210\n"
        tx_cpm += "tx2\t50\t55\t48\t100\t110\t105\n"
        tx_cpm += "tx3\t75\t80\t72\t150\t160\t155\n"
        _write(cohort / "transcript_cpm.tsv", tx_cpm)

    samples = tmp_path / "samples"
    samples.mkdir()
    mod_summaries = tmp_path / "mod_summaries"
    mod_summaries.mkdir()
    sqanti = tmp_path / "sqanti"
    sqanti.mkdir()

    de_dir = None
    if de_qc is not None:
        de_dir = tmp_path / "de_analysis"
        de_dir.mkdir()
        _write(de_dir / "de_qc_stats.json", json.dumps(de_qc))
        for contrast_name in de_qc.get("contrasts", {}):
            contrast_dir = de_dir / contrast_name
            contrast_dir.mkdir()
            _write(
                contrast_dir / "results_dge.tsv",
                (
                    "GENEID\tnewGeneClass\tgene_name\tbaseMean\tlog2FoldChange\tlfcSE"
                    "\tstat\tpvlaue\tpadj\n"
                    "gene1\tannotation\tgene2\t100.0\t1.0\t-0.02\t-8.5\t0.001\t0.05\n"
                )
            )
            _write(
                contrast_dir / "results_dtu_transcript.tsv",
                "featureID\tgroupID\tlog2FoldChange\tpvalue\tpadj\texonBaseMean\n"
                "tx1\tgene1\t1.0\t0.01\t0.05\t20.0\n"
            )

    out_report = tmp_path / "wf-transcriptomes-report.html"
    argv = [
        str(out_report),
        "--metadata",
        str(metadata),
        "--cohort_dir",
        str(cohort),
        "--ref_summary",
        str(reference / "annotation_reference_summary.json"),
        "--samples_dir",
        str(samples),
        "--mod_summary_dir",
        str(mod_summaries),
        "--sqanti_dir",
        str(sqanti),
        "--versions",
        str(versions),
        "--params",
        str(params),
    ]
    if de_dir is not None:
        argv.extend(["--de_dir", str(de_dir)])
    return report.argparser().parse_args(argv), out_report


def test_report_main_accepts_optional_file_sentinels(monkeypatch, tmp_path):
    """The report entry point should tolerate null-object sentinel files."""
    tables = []
    monkeypatch.setattr(report.labs, "LabsReport", _FakeReport)
    monkeypatch.setattr(report, "Tabs", _FakeTabs)
    monkeypatch.setattr(report, "p", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "pre", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "_create_warning_banner", lambda *args, **kwargs: None)
    monkeypatch.setattr(report.fastcat, "SeqSummary", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report.DataTable,
        "from_pandas",
        staticmethod(lambda table, *args, **kwargs: tables.append(table.copy())),
    )

    metadata = _write(
        tmp_path / "metadata.json",
        json.dumps([{"alias": "sampleA", "has_stats": False}]),
    )
    params = _write(tmp_path / "params.json", "{}")
    versions = tmp_path / "versions"
    versions.mkdir()
    _write(versions / "versions.txt", "tool,1.0\n")

    cohort = tmp_path / "cohort"
    cohort.mkdir()
    reference = cohort / "reference"
    reference.mkdir()
    _write(
        reference / "annotation_reference_summary.json",
        json.dumps({
            "seqname_overlap": ["chr1"],
            "only_in_annotation": ["chrMissing"],
            "only_in_reference": ["chrExtra"],
            "annotation": {
                "kept_records": 10,
                "excluded_unstranded_records": 2,
                "sanitised_attribute_records": 1,
                "unstranded_examples": ["chr1\tsim\ttranscript\t1\t4\t.\t.\t."],
            },
            "reference_build_hints": ["GRCh38"],
            "annotation_build_hints": [],
            "reference_provider_hints": [],
            "annotation_provider_hints": [],
            "warnings": ["Warning: Some seqnames are present in the annotation."],
        }),
    )

    samples = tmp_path / "samples"
    samples.mkdir()
    (samples / "OPTIONAL_FILE").touch()

    mod_summaries = tmp_path / "mod_summaries"
    mod_summaries.mkdir()
    (mod_summaries / "OPTIONAL_FILE").touch()

    sqanti = tmp_path / "sqanti"
    sqanti.mkdir()
    (sqanti / "OPTIONAL_FILE").touch()

    out_report = tmp_path / "wf-transcriptomes-report.html"
    args = report.argparser().parse_args(
        [
            str(out_report),
            "--metadata",
            str(metadata),
            "--cohort_dir",
            str(cohort),
            "--ref_summary",
            str(cohort / "reference" / "annotation_reference_summary.json"),
            "--samples_dir",
            str(samples),
            "--mod_summary_dir",
            str(mod_summaries),
            "--sqanti_dir",
            str(sqanti),
            "--versions",
            str(versions),
            "--params",
            str(params),
        ]
    )

    report.main(args)

    assert out_report.exists()
    assert any("Overlapping seqnames" in table.to_string() for table in tables)
    assert any(
        "Annotation attributes sanitised" in table.to_string() for table in tables
    )
    assert any("GRCh38" in table.to_string() for table in tables)


def test_report_main_renders_modified_base_summary_tables(monkeypatch, tmp_path):
    """Per-sample modified base summaries should render as report tables."""
    tables = []
    raw_calls = []

    monkeypatch.setattr(report.labs, "LabsReport", _FakeReport)
    monkeypatch.setattr(report, "Tabs", _FakeTabs)
    monkeypatch.setattr(report, "p", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "pre", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "raw", lambda value: raw_calls.append(value) or value)
    monkeypatch.setattr(report, "_create_warning_banner", lambda *args, **kwargs: None)
    monkeypatch.setattr(report.fastcat, "SeqSummary", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report.DataTable,
        "from_pandas",
        staticmethod(lambda table, *args, **kwargs: tables.append(table.copy())),
    )

    args, out_report = _build_report_args(tmp_path)
    mod_summaries = tmp_path / "mod_summaries"
    _write(
        mod_summaries / "not_the_sample_name.mods.summary.tsv",
        (
            "sample\tfull_mod_code\tmod_code\tmod_label\tvalid_coverage\t"
            "modified_calls\tcanonical_calls\tother_calls\tdelete_calls\t"
            "fail_calls\tdiff_calls\tnocall_calls\tmodification_percent\n"
            "sampleA\tA:a\ta\tm6A\t12\t6\t6\t0\t0\t3\t0\t0\t50.00\n"
        ),
    )
    _write(
        mod_summaries / "sampleB.mods.summary.tsv",
        (
            "sample\tfull_mod_code\tmod_code\tmod_label\tvalid_coverage\t"
            "modified_calls\tcanonical_calls\tother_calls\tdelete_calls\t"
            "fail_calls\tdiff_calls\tnocall_calls\tmodification_percent\n"
            "sampleB\tC:m\tm\tm5C\t20\t5\t15\t0\t0\t0\t0\t0\t25.00\n"
        ),
    )

    report.main(args)

    assert out_report.exists()
    matrix_html = next(
        value for value in raw_calls
        if "<table class='mod-summary-matrix'>" in value
    )
    assert "sampleA" in matrix_html
    assert "sampleB" in matrix_html
    assert "m6A" in matrix_html
    assert "m5C" in matrix_html
    assert "50.00%" in matrix_html
    assert "25.00%" in matrix_html
    assert "6 modified" in matrix_html
    assert "12 valid" in matrix_html
    assert "5 modified" in matrix_html
    assert "20 valid" in matrix_html


def test_report_main_handles_degenerate_bambu_qc_and_read_summary(
    monkeypatch,
    tmp_path,
):
    """Tiny/empty stats should not crash the report rendering path."""
    tables = []
    banners = []

    monkeypatch.setattr(report.labs, "LabsReport", _FakeReport)
    monkeypatch.setattr(report, "Tabs", _FakeTabs)
    monkeypatch.setattr(report, "p", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "pre", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report.fastcat,
        "SeqSummary",
        lambda *args, **kwargs: (_ for _ in ()).throw(KeyError(1)),
    )
    monkeypatch.setattr(
        report,
        "_create_warning_banner",
        lambda message, level="warning": banners.append((level, message)),
    )
    monkeypatch.setattr(
        report.DataTable,
        "from_pandas",
        staticmethod(lambda table, *args, **kwargs: tables.append(table.copy())),
    )

    metadata = _write(
        tmp_path / "metadata.json",
        json.dumps([{"alias": "sampleA", "has_stats": True}]),
    )
    params = _write(tmp_path / "params.json", "{}")
    versions = tmp_path / "versions"
    versions.mkdir()
    _write(versions / "versions.txt", "tool,1.0\n")

    cohort = tmp_path / "cohort"
    cohort.mkdir()
    reference = cohort / "reference"
    reference.mkdir()
    _write(
        reference / "annotation_reference_summary.json",
        json.dumps(
            {
                "seqname_overlap": ["chr1"],
                "only_in_annotation": [],
                "only_in_reference": [],
                "annotation": {
                    "kept_records": 10,
                    "excluded_unstranded_records": 0,
                    "sanitised_attribute_records": 0,
                },
                "warnings": [],
            }
        ),
    )
    samples = tmp_path / "samples"
    samples.mkdir()
    sample_a = samples / "sampleA"
    sample_a.mkdir()
    _write(
        sample_a / "bambu_qc_stats.json",
        json.dumps(
            {
                "samples": 1,
                "library_sizes": {"sampleA": 0},
                "min_library_size": 0,
                "max_library_size": 0,
                "median_library_size": 0,
                "library_size_ratio": "NA",
                "total_transcripts_before_filter": 0,
                "total_transcripts_after_filter": 0,
                "transcripts_filtered": 0,
                "median_transcripts_detected": 0,
                "total_genes_after_filter": 0,
                "transcriptome_mode": "discover",
                "ndr_used": "automatic",
            }
        ),
    )

    sqanti = tmp_path / "sqanti"
    sqanti.mkdir()
    (sqanti / "OPTIONAL_FILE").touch()

    alignment_stats = tmp_path / "alignment_stats"
    alignment_stats.mkdir()

    out_report = tmp_path / "wf-transcriptomes-report.html"
    args = report.argparser().parse_args(
        [
            str(out_report),
            "--metadata",
            str(metadata),
            "--stats",
            str(alignment_stats),
            "--cohort_dir",
            str(cohort),
            "--ref_summary",
            str(cohort / "reference" / "annotation_reference_summary.json"),
            "--samples_dir",
            str(samples),
            "--sqanti_dir",
            str(sqanti),
            "--versions",
            str(versions),
            "--params",
            str(params),
        ]
    )

    report.main(args)

    assert out_report.exists()
    assert any(
        level == "info" and "Read summary plots could not be rendered" in message
        for level, message in banners
    )
    assert any(
        "Library size ratio (max/min)" in table.to_string()
        and "N/A" in table.to_string()
        for table in tables
    )


def test_report_main_uses_cohort_bambu_qc_for_multi_sample_inputs(
    monkeypatch,
    tmp_path,
):
    """Multi-sample runs should render bambu QC from cohort-level outputs."""
    tables = []

    monkeypatch.setattr(report.labs, "LabsReport", _FakeReport)
    monkeypatch.setattr(report, "Tabs", _FakeTabs)
    monkeypatch.setattr(report, "p", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "pre", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "_create_warning_banner", lambda *args, **kwargs: None)
    monkeypatch.setattr(report.fastcat, "SeqSummary", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report.DataTable,
        "from_pandas",
        staticmethod(lambda table, *args, **kwargs: tables.append(table.copy())),
    )

    metadata = _write(
        tmp_path / "metadata.json",
        json.dumps(
            [
                {"alias": "sampleA", "has_stats": False},
                {"alias": "sampleB", "has_stats": False},
            ]
        ),
    )
    params = _write(tmp_path / "params.json", "{}")
    versions = tmp_path / "versions"
    versions.mkdir()
    _write(versions / "versions.txt", "tool,1.0\n")

    cohort = tmp_path / "cohort"
    cohort.mkdir()
    reference = cohort / "reference"
    reference.mkdir()
    _write(
        reference / "annotation_reference_summary.json",
        json.dumps(
            {
                "seqname_overlap": ["chr1"],
                "only_in_annotation": [],
                "only_in_reference": [],
                "annotation": {
                    "kept_records": 10,
                    "excluded_unstranded_records": 0,
                    "sanitised_attribute_records": 0,
                },
                "warnings": [],
            }
        ),
    )
    _write(
        cohort / "bambu_qc_stats.json",
        json.dumps(
            {
                "samples": 2,
                "library_sizes": {"sampleA": 1200, "sampleB": 900},
                "min_library_size": 900,
                "max_library_size": 1200,
                "median_library_size": 1050,
                "library_size_ratio": 1.3333,
                "total_transcripts_before_filter": 100,
                "total_transcripts_after_filter": 80,
                "transcripts_filtered": 20,
                "median_transcripts_detected": 70,
                "total_genes_after_filter": 60,
                "transcriptome_mode": "discover",
                "ndr_used": 0.1,
            }
        ),
    )

    samples = tmp_path / "samples"
    samples.mkdir()
    sample_a = samples / "sampleA"
    sample_a.mkdir()
    _write(
        sample_a / "bambu_qc_stats.json",
        json.dumps({"samples": 1, "library_sizes": {"sampleA": 5}}),
    )
    _write(
        sample_a / "transcript_metadata.tsv",
        "TXNAME\tGENEID\n"
        "tx1\tgene1\n"
        "tx2\tgene1\n"
        "tx3\tgene2\n",
    )
    sample_b = samples / "sampleB"
    sample_b.mkdir()
    _write(
        sample_b / "bambu_qc_stats.json",
        json.dumps({"samples": 1, "library_sizes": {"sampleB": 7}}),
    )
    _write(
        sample_b / "transcript_metadata.tsv",
        "TXNAME\tGENEID\n"
        "txA\tgeneA\n"
        "txB\tgeneB\n",
    )

    sqanti = tmp_path / "sqanti"
    sqanti.mkdir()
    (sqanti / "OPTIONAL_FILE").touch()

    out_report = tmp_path / "wf-transcriptomes-report.html"
    args = report.argparser().parse_args(
        [
            str(out_report),
            "--metadata",
            str(metadata),
            "--cohort_dir",
            str(cohort),
            "--ref_summary",
            str(cohort / "reference" / "annotation_reference_summary.json"),
            "--samples_dir",
            str(samples),
            "--sqanti_dir",
            str(sqanti),
            "--versions",
            str(versions),
            "--params",
            str(params),
        ]
    )

    report.main(args)

    assert out_report.exists()
    assert any(
        "Samples analyzed" in table.to_string() and "2" in table.to_string()
        for table in tables
    )
    assert any(
        "Library size ratio (max/min)" in table.to_string()
        and "1.33x" in table.to_string()
        for table in tables
    )
    assert any(
        "Sample" in table.columns
        and "Library Size" in table.columns
        and {"sampleA", "sampleB"}.issubset(set(table["Sample"].tolist()))
        and "1,200" in table.to_string()
        and "900" in table.to_string()
        for table in tables
    )
    per_sample_metric_tables = [
        table
        for table in tables
        if list(table.columns) == ["Metric", "Value"]
        and set(table["Metric"].tolist()) == {"Transcripts", "Genes"}
    ]
    assert any(
        set(zip(table["Metric"], table["Value"])) == {
            ("Transcripts", 3),
            ("Genes", 2),
        }
        for table in per_sample_metric_tables
    )
    assert any(
        set(zip(table["Metric"], table["Value"])) == {
            ("Transcripts", 2),
            ("Genes", 2),
        }
        for table in per_sample_metric_tables
    )


def test_report_main_renders_statistical_methods_and_warnings(
    monkeypatch,
    tmp_path,
):
    """DE/DTU QC report renders fallback methods and warning banners."""
    tables = []
    headings = []
    banners = []

    monkeypatch.setattr(report.labs, "LabsReport", _FakeReport)
    monkeypatch.setattr(report, "Tabs", _FakeTabs)
    monkeypatch.setattr(report, "p", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "pre", lambda *args, **kwargs: None)
    monkeypatch.setattr(report.fastcat, "SeqSummary", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report,
        "h4",
        lambda label: (headings.append(label), _NullContext())[1],
    )
    monkeypatch.setattr(
        report,
        "_create_warning_banner",
        lambda message, level="warning": banners.append((level, message)),
    )
    monkeypatch.setattr(
        report.DataTable,
        "from_pandas",
        staticmethod(lambda table, *args, **kwargs: tables.append(table.copy())),
    )

    de_qc = {
        "total_samples": 6,
        "condition_column": "condition",
        "reference_level": "control",
        "covariates": ["batch"],
        "num_contrasts": 2,
        "sample_size_warnings": "none",
        "samples_per_group": {"control": 3, "treated": 3},
        "contrasts": {
            "condition_treated_vs_control": {
                "n_target": 3,
                "n_reference": 3,
                "dge_significant_fdr05": 10,
                "dge_upregulated": 6,
                "dge_downregulated": 4,
                "dtu_status": "SUCCESS",
                "dtu_significant_genes": 2,
                "deseq2_dispersion_fallback": {
                    "applied": True,
                    "method_used": "gene-wise",
                    "reason": "recoverable",
                    "diagnostic_file": (
                        "DESeq2_dispersion_fallback_"
                        "condition_treated_vs_control.txt"
                    ),
                },
                "dexseq_dispersion_method": "local",
                "dexseq_covariates_dropped": ["batch"],
            },
            "condition_treated2_vs_control": {
                "n_target": 3,
                "n_reference": 3,
                "dge_significant_fdr05": 4,
                "dge_upregulated": 3,
                "dge_downregulated": 1,
                "dtu_status": "FAILED",
            },
        },
    }

    args, out_report = _build_report_args(tmp_path, de_qc=de_qc)
    report.main(args)

    assert out_report.exists()
    assert "Statistical Methods & Warnings" in headings
    assert any("DESeq2 dispersion" in table.columns for table in tables)
    assert any(
        "gene-wise (fallback)" in table.to_string()
        for table in tables
        if "DESeq2 dispersion" in table.columns
    )
    assert any(
        "batch" in table.to_string()
        for table in tables
        if "DEXSeq covariates dropped" in table.columns
    )
    assert any("gene-wise dispersion fallback" in msg.lower() for _, msg in banners)
    assert any("covariates dropped" in msg.lower() for _, msg in banners)
    assert any(
        level == "danger" and "DTU Analysis Failed" in msg
        for level, msg in banners
    )


def test_report_main_tolerates_missing_statistical_fields(monkeypatch, tmp_path):
    """Older DE QC JSON without new fallback fields should still render."""
    tables = []
    headings = []

    monkeypatch.setattr(report.labs, "LabsReport", _FakeReport)
    monkeypatch.setattr(report, "Tabs", _FakeTabs)
    monkeypatch.setattr(report, "p", lambda *args, **kwargs: None)
    monkeypatch.setattr(report, "pre", lambda *args, **kwargs: None)
    monkeypatch.setattr(report.fastcat, "SeqSummary", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report,
        "h4",
        lambda label: (headings.append(label), _NullContext())[1],
    )
    monkeypatch.setattr(report, "_create_warning_banner", lambda *args, **kwargs: None)
    monkeypatch.setattr(
        report.DataTable,
        "from_pandas",
        staticmethod(lambda table, *args, **kwargs: tables.append(table.copy())),
    )

    legacy_de_qc = {
        "total_samples": 4,
        "condition_column": "condition",
        "reference_level": "control",
        "covariates": "none",
        "num_contrasts": 1,
        "sample_size_warnings": "none",
        "samples_per_group": {"control": 2, "treated": 2},
        "contrasts": {
            "condition_treated_vs_control": {
                "n_target": 2,
                "n_reference": 2,
                "dge_significant_fdr05": 1,
                "dge_upregulated": 1,
                "dge_downregulated": 0,
                "dtu_status": "SUCCESS",
            }
        },
    }

    args, out_report = _build_report_args(tmp_path, de_qc=legacy_de_qc)
    report.main(args)

    assert out_report.exists()
    assert "Statistical Methods & Warnings" in headings
    method_tables = [
        table
        for table in tables
        if "DESeq2 dispersion" in table.columns
    ]
    assert method_tables
    assert "parametric" in method_tables[0].to_string()


def test_round_de_table_formats_supported_numeric_columns():
    """DE/DTU preview tables should round known numeric columns for display."""
    table = pd.DataFrame(
        {
            "GENEID": ["gene1"],
            "baseMean": [123.4567],
            "log2FoldChange": [0.00001234],
            "lfcSE": [0.98765],
            "stat": [-45.6789],
            "pvalue": [0.00001234],
            "padj": [0.123456],
            "other": [7.89123],
        }
    )

    rounded = report._round_de_table(table)

    assert rounded.loc[0, "baseMean"] == "123.457"
    assert rounded.loc[0, "log2FoldChange"] == "1.234e-05"
    assert rounded.loc[0, "lfcSE"] == "0.988"
    assert rounded.loc[0, "stat"] == "-45.679"
    assert rounded.loc[0, "pvalue"] == "1.234e-05"
    assert rounded.loc[0, "padj"] == "0.123"
    assert rounded.loc[0, "other"] == 7.89123


def test_round_de_table_uses_scientific_notation_when_fixed_decimal_would_zero():
    """Tiny non-zero values should not display as 0.000 in DE/DTU tables."""
    table = pd.DataFrame(
        {
            "GENEID": ["gene1"],
            "pvalue": [0.00012],
            "padj": [0.00049],
        }
    )

    rounded = report._round_de_table(table)

    assert rounded.loc[0, "pvalue"] == "1.200e-04"
    assert rounded.loc[0, "padj"] == "4.900e-04"


def test_sorted_transcript_abundance_table_uses_sample_alias_columns(tmp_path):
    """Transcript abundance sorting should use only real sample alias columns."""
    tx_counts_file = _write(
        tmp_path / "transcript_counts.tsv",
        (
            "TXNAME\tsampleA\tsampleB\tgene_length\tannotation_score\n"
            "tx_low\t1\t1\t10000\t999\n"
            "tx_high\t5\t5\t10\t1\n"
            "tx_mid\t2\t2\t5000\t500\n"
        ),
    )

    sorted_table = report._sorted_transcript_abundance_table(
        tx_counts_file, ["sampleA", "sampleB"]
    )

    assert sorted_table["TXNAME"].tolist() == ["tx_high", "tx_mid", "tx_low"]
