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


def test_pychopper_tables_uses_sample_directory_names(tmp_path):
    """Pychopper summaries should be keyed by sample alias."""
    pychopper_dir = tmp_path / "pychopper"
    sample_dir = pychopper_dir / "sampleA_pychopper_output"
    sample_dir.mkdir(parents=True)
    _write(
        sample_dir / "pychopper_summary.tsv",
        "Classification\tValue\nFull length\t10\nUnclassified\t2\n",
    )

    tables = report._pychopper_tables(pychopper_dir)

    assert list(tables) == ["sampleA"]
    assert list(tables["sampleA"]["Classification"]) == [
        "Full length",
        "Unclassified",
    ]
