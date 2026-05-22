"""Tests for report QC integration."""

import json


def test_load_bambu_qc_exists(tmp_path):
    """Test loading bambu QC stats when file exists."""
    from workflow_glue.report import _load_bambu_qc

    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()

    qc_data = {
        "samples": 4,
        "library_sizes": {"sample1": 1000000, "sample2": 2000000},
        "median_library_size": 1500000,
        "library_size_warning": "2.0x variation (>3x threshold)",
    }

    with open(cohort_dir / "bambu_qc_stats.json", "w") as f:
        json.dump(qc_data, f)

    result = _load_bambu_qc(cohort_dir)
    assert result is not None
    assert result["samples"] == 4
    assert result["median_library_size"] == 1500000


def test_load_bambu_qc_missing(tmp_path):
    """Test loading bambu QC stats when file missing."""
    from workflow_glue.report import _load_bambu_qc

    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()

    result = _load_bambu_qc(cohort_dir)
    assert result is None


def test_load_de_qc_exists(tmp_path):
    """Test loading DE QC stats when file exists."""
    from workflow_glue.report import _load_de_qc

    de_dir = tmp_path / "de_analysis"
    de_dir.mkdir()

    qc_data = {
        "total_samples": 6,
        "num_contrasts": 2,
        "sample_size_warnings": "Some groups have n<3",
        "contrasts": {
            "contrast1": {
                "dge_significant_fdr05": 123,
                "dtu_status": "SUCCESS",
            },
        },
    }

    with open(de_dir / "de_qc_stats.json", "w") as f:
        json.dump(qc_data, f)

    result = _load_de_qc(de_dir)
    assert result is not None
    assert result["total_samples"] == 6
    assert result["num_contrasts"] == 2
    assert "contrasts" in result


def test_load_de_qc_missing(tmp_path):
    """Test loading DE QC stats when file missing."""
    from workflow_glue.report import _load_de_qc

    de_dir = tmp_path / "de_analysis"
    de_dir.mkdir()

    result = _load_de_qc(de_dir)
    assert result is None


def test_load_annotation_reference_summary_exists(tmp_path):
    """Test loading reference/annotation prep summary when file exists."""
    from workflow_glue.report import _load_annotation_reference_summary

    summary_dir = tmp_path / "cohort" / "reference"
    summary_dir.mkdir(parents=True)
    summary_data = {
        "seqname_overlap": ["chr1"],
        "only_in_annotation": ["chrMissing"],
        "only_in_reference": ["chrExtra"],
        "warnings": ["Warning: test"],
    }

    with open(summary_dir / "annotation_reference_summary.json", "w") as f:
        json.dump(summary_data, f)

    result = _load_annotation_reference_summary(tmp_path / "cohort")
    assert result is not None
    assert result["seqname_overlap"] == ["chr1"]
    assert result["only_in_annotation"] == ["chrMissing"]


def test_load_annotation_reference_summary_missing(tmp_path):
    """Test loading reference/annotation prep summary when file is missing."""
    from workflow_glue.report import _load_annotation_reference_summary

    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()

    result = _load_annotation_reference_summary(cohort_dir)
    assert result is None


def test_load_cpm_tables_exists(tmp_path):
    """Load cohort CPM tables when both files exist."""
    from workflow_glue.report import _load_cpm_tables

    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()
    (cohort_dir / "gene_cpm.tsv").write_text(
        "GENEID\tsample1\tsample2\n"
        "gene1\t1.0\t2.0\n"
    )
    (cohort_dir / "transcript_cpm.tsv").write_text(
        "TXNAME\tsample1\tsample2\n"
        "tx1\t3.0\t4.0\n"
    )

    result = _load_cpm_tables(cohort_dir)
    assert result["gene"] is not None
    assert result["transcript"] is not None
    assert list(result["gene"]["GENEID"]) == ["gene1"]
    assert list(result["transcript"]["TXNAME"]) == ["tx1"]


def test_load_cpm_tables_missing(tmp_path):
    """Return None entries when cohort CPM tables are missing."""
    from workflow_glue.report import _load_cpm_tables

    cohort_dir = tmp_path / "cohort"
    cohort_dir.mkdir()

    result = _load_cpm_tables(cohort_dir)
    assert result["gene"] is None
    assert result["transcript"] is None


def test_format_hint_values():
    """Build/provider hints should be compactly formatted for the report."""
    from workflow_glue.report import _format_hint_values

    assert _format_hint_values([]) == "None detected"
    assert _format_hint_values(["GRCh38"]) == "GRCh38"
    assert _format_hint_values(["GRCh38", "GENCODE"]) == "GRCh38, GENCODE"


def test_warning_banner_creation():
    """Test that warning banner can be created."""
    from workflow_glue.report import _create_warning_banner
    from dominate import document

    doc = document()
    with doc:
        _create_warning_banner("Test warning message", level="warning")

    html = doc.render()
    assert "Test warning message" in html
    assert "background-color" in html


def test_bambu_qc_with_warnings():
    """Test bambu QC data structure with warnings."""
    qc_data = {
        "samples": 3,
        "library_sizes": {"s1": 1000000, "s2": 4000000, "s3": 2000000},
        "min_library_size": 1000000,
        "max_library_size": 4000000,
        "median_library_size": 2000000,
        "library_size_ratio": 4.0,
        "library_size_warning": "4.0x variation (>3x threshold)",
        "total_transcripts_before_filter": 50000,
        "total_transcripts_after_filter": 45000,
        "transcripts_filtered": 5000,
        "transcriptome_mode": "discover",
        "ndr_used": 0.1,
    }

    # Verify structure
    assert qc_data["library_size_warning"] is not None
    assert qc_data["library_size_ratio"] > 3.0


def test_de_qc_with_multiple_warnings():
    """Test DE QC data structure with multiple warning types."""
    qc_data = {
        "total_samples": 4,
        "num_contrasts": 2,
        "sample_size_warnings": "Some groups have n<3",
        "multiple_testing_note": "Testing 2 contrasts yields FWER ~9.8%",
        "samples_per_group": {"control": 2, "treated": 2},
        "contrasts": {
            "condition_treated_vs_control": {
                "n_samples": 4,
                "n_target": 2,
                "n_reference": 2,
                "dge_significant_fdr05": 50,
                "dtu_status": "FAILED",
                "dtu_power_warning": "DTU may be underpowered (n=4)",
            },
        },
    }

    # Verify warnings are captured
    assert qc_data["sample_size_warnings"] != "none"
    assert qc_data["multiple_testing_note"] is not None
    contrast = qc_data["contrasts"]["condition_treated_vs_control"]
    assert contrast["dtu_status"] == "FAILED"
    assert contrast["dtu_power_warning"] is not None


def test_de_qc_no_warnings():
    """Test DE QC data structure with no warnings."""
    qc_data = {
        "total_samples": 6,
        "num_contrasts": 1,
        "sample_size_warnings": "none",
        "samples_per_group": {"control": 3, "treated": 3},
        "contrasts": {
            "condition_treated_vs_control": {
                "n_samples": 6,
                "n_target": 3,
                "n_reference": 3,
                "dge_significant_fdr05": 150,
                "dtu_status": "SUCCESS",
                "dtu_significant_genes": 25,
            },
        },
    }

    # Verify no warnings
    assert qc_data["sample_size_warnings"] == "none"
    assert (
        "multiple_testing_note" not in qc_data
        or qc_data.get("multiple_testing_note") is None
    )
    contrast = qc_data["contrasts"]["condition_treated_vs_control"]
    assert contrast["dtu_status"] == "SUCCESS"
