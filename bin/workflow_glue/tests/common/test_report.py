"""Tests for the workflow report entry point."""

import json
from pathlib import Path

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

    samples = tmp_path / "samples"
    samples.mkdir()
    sqanti = tmp_path / "sqanti"
    sqanti.mkdir()
    alignment_stats = tmp_path / "alignment_stats"
    alignment_stats.mkdir()
    (samples / "OPTIONAL_FILE").touch()
    (sqanti / "OPTIONAL_FILE").touch()
    (alignment_stats / "OPTIONAL_FILE").touch()

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
                "GENEID\tlog2FoldChange\tpadj\n"
                "gene1\t1.0\t0.01\n",
            )
            _write(
                contrast_dir / "results_dtu_transcript.tsv",
                "featureID\tgroupID\tpadj\n"
                "tx1\tgene1\t0.05\n",
            )

    out_report = tmp_path / "wf-transcriptomes-report.html"
    argv = [
        str(out_report),
        "--metadata",
        str(metadata),
        "--alignment_stats_dir",
        str(alignment_stats),
        "--cohort_dir",
        str(cohort),
        "--samples_dir",
        str(samples),
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

    sqanti = tmp_path / "sqanti"
    sqanti.mkdir()
    (sqanti / "OPTIONAL_FILE").touch()

    alignment_stats = tmp_path / "alignment_stats"
    alignment_stats.mkdir()
    (alignment_stats / "OPTIONAL_FILE").touch()

    out_report = tmp_path / "wf-transcriptomes-report.html"
    args = report.argparser().parse_args(
        [
            str(out_report),
            "--metadata",
            str(metadata),
            "--alignment_stats_dir",
            str(alignment_stats),
            "--cohort_dir",
            str(cohort),
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
    assert any("Overlapping seqnames" in table.to_string() for table in tables)
    assert any(
        "Annotation attributes sanitised" in table.to_string() for table in tables
    )
    assert any("GRCh38" in table.to_string() for table in tables)


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
    _write(
        cohort / "bambu_qc_stats.json",
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

    samples = tmp_path / "samples"
    samples.mkdir()
    (samples / "OPTIONAL_FILE").touch()

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
            "--alignment_stats_dir",
            str(alignment_stats),
            "--stats",
            str(alignment_stats),
            "--cohort_dir",
            str(cohort),
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
