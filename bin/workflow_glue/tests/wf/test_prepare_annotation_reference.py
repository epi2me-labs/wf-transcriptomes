"""Tests for annotation/reference preparation."""

import gzip
from pathlib import Path

import pytest
from workflow_glue import get_components
from workflow_glue import prepare_annotation_reference


def _write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def _write_gzip(path, text):
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)
    return path


SIMPLE_REFERENCE = ">chr1\nAAAA\n"

SIMPLE_GTF = 'chr1\tsim\ttranscript\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'

GFF3_INPUT = """##gff-version 3
chr1\tsim\tgene\t1\t10\t.\t+\t.\tID=gene1
chr1\tsim\tmRNA\t1\t10\t.\t+\t.\tID=transcript1;Parent=gene1
chr1\tsim\texon\t1\t5\t.\t+\t.\tParent=transcript1
chr1\tsim\texon\t6\t10\t.\t+\t.\tParent=transcript1
"""

GFF3_INPUT_WITH_UNKNOWN_STRAND = """##gff-version 3
chr1\tsim\tgene\t1\t10\t.\t?\t.\tID=gene1
chr1\tsim\tmRNA\t1\t10\t.\t+\t.\tID=transcript1;Parent=gene1
chr1\tsim\texon\t1\t5\t.\t+\t.\tParent=transcript1
chr1\tsim\texon\t6\t10\t.\t+\t.\tParent=transcript1
"""

NCBI_REFERENCE = """>chr1 GRCh38 RefSeq primary assembly
AAAA
>chrExtra
TTTT
"""

NCBI_GTF = (
    "#!genome-build GRCh38\n"
    "#!genome-version GRCh38.p14 RefSeq\n"
    'chr1\tRefSeq\ttranscript\t1\t4\t.\t-\t.\tgene_id "transcript_id"; '
    'transcript_id "rna-XM_000001.1";\n'
    'chr1\tRefSeq\texon\t1\t4\t.\t-\t.\ttranscript_id "rna-XM_000001.1"; '
    'exon_number "1";\n'
)

MIXED_STRANDED_GTF = (
    "# header\n"
    'chr1\tsim\ttranscript\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
    'chrMissing\tsim\ttranscript\t1\t4\t.\t.\t.\tgene_id "g2"; '
    'transcript_id "t2";\n'
)

MOUSE_REFERENCE = ">chr1 dna:GRCm39 primary assembly\nAAAA\n"

MOUSE_GTF = (
    "#!genome-build GRCm39\n"
    'chr1\tGENCODE\ttranscript\t1\t4\t.\t+\t.\tgene_id "g1"; '
    'transcript_id "t1";\n'
)


@pytest.mark.parametrize(
    "ref_data,ref_name,annot_data,annot_name",
    [
        (SIMPLE_REFERENCE, "reference.fa", SIMPLE_GTF,
         "annotation.gtf"),
        (SIMPLE_REFERENCE, "reference.fa", SIMPLE_GTF,
         "annotation.gtf.gz"),
        (SIMPLE_REFERENCE, "reference.fa", GFF3_INPUT,
         "annotation.gff3.gz"),
        (NCBI_REFERENCE, "reference.fna", NCBI_GTF,
         "annotation.gtf.gz"),
        (SIMPLE_REFERENCE + ">chrExtra\nTTTT\n", "reference.fa",
         MIXED_STRANDED_GTF, "annotation.gtf"),
        (MOUSE_REFERENCE, "GRCm39.genome.fa", MOUSE_GTF,
         "gencode.vM33.annotation.gtf"),
    ],
    ids=[
        "simple_gtf", "gzipped_gtf", "gff3_input",
        "ncbi_format", "mixed_stranded", "mouse_gencode"],
)
def test_prepare_all_formats_produce_valid_outputs(
    ref_data, ref_name, annot_data, annot_name, tmp_path
):
    """All supported input formats produce annotation.gtf and reference.fasta."""
    # write reference
    ref_path = tmp_path / ref_name
    _write(ref_path, ref_data)

    # write annotation
    annot_path = tmp_path / annot_name
    if annot_name.endswith(".gz"):
        _write_gzip(annot_path, annot_data)
    else:
        _write(annot_path, annot_data)

    out_dir = tmp_path / "prepared"
    summary = prepare_annotation_reference.prepare_annotation_reference(
        annot_path, ref_path, out_dir
    )

    # output files must exist
    assert (out_dir / "annotation.gtf").exists()
    assert (ref_path).exists()
    assert (out_dir / "annotation_reference_summary.json").exists()

    # paths returned in summary must match
    assert summary["annotation"]["prepared"] == str(out_dir / "annotation.gtf")
    assert summary["reference"]["input"] == str(ref_path)

    # outputs must be valid for downstream tools (Bambu, SQANTI)
    reference = Path(summary["reference"]["input"])
    annotation = Path(summary["annotation"]["prepared"])
    assert reference.read_text(encoding="utf-8").startswith(">")
    assert summary["analysis_seqnames"]["has_overlap"] is True

    # all annotation records must be stranded with gene_id and transcript_id
    records = [
        line.rstrip("\n").split("\t")
        for line in annotation.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert records
    for fields in records:
        assert len(fields) == 9
        assert fields[6] in {"+", "-"}
        assert 'gene_id "' in fields[8]
        assert 'transcript_id "' in fields[8]


def test_ncbi_format_sanitises_gene_ids(tmp_path):
    """NCBI gene_id='transcript_id' should be sanitised to actual transcript_id."""
    reference = _write(tmp_path / "reference.fna", NCBI_REFERENCE)
    annotation = _write_gzip(tmp_path / "annotation.gtf.gz", NCBI_GTF)

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    prepared_text = Path(summary["annotation"]["prepared"]).read_text(encoding="utf-8")
    assert 'gene_id "transcript_id"' not in prepared_text
    assert prepared_text.count('gene_id "rna-XM_000001.1"') == 2
    assert summary["annotation"]["sanitised_attribute_records"] == 2
    assert any(
        "Sanitised 2 annotation records" in warning
        for warning in summary["warnings"]
    )


def test_gff3_conversion_via_gffread(tmp_path):
    """GFF3 inputs invoke gffread and convert to GTF format."""
    reference = _write(tmp_path / "reference.fa", SIMPLE_REFERENCE)
    annotation = _write(tmp_path / "annotation.gff3", GFF3_INPUT)

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    assert summary["annotation"]["was_gff"] is True

    prepared_text = Path(summary["annotation"]["prepared"]).read_text(encoding="utf-8")
    assert 'gene_id "gene1"' in prepared_text
    assert 'transcript_id "transcript1"' in prepared_text
    assert "\texon\t" in prepared_text


def test_normalise_unknown_gff_strands_rewrites_question_mark(tmp_path):
    """Unknown GFF strand '?' should be rewritten to '.'."""
    input_path = _write(tmp_path / "annotation.gff3", GFF3_INPUT_WITH_UNKNOWN_STRAND)
    output_path = tmp_path / "annotation_sanitised.gff"
    init_annotation = _write(tmp_path / "annotation_init.gtf", SIMPLE_GTF)
    ann = prepare_annotation_reference.Annotation(
        init_annotation,
        tmp_path,
    )

    used_path, normalised = ann._normalise_unknown_gff_strands(input_path, output_path)

    assert used_path == output_path
    assert normalised == 1
    output_text = output_path.read_text(encoding="utf-8")
    assert "\t?\t" not in output_text
    assert "\t.\t" in output_text


def test_gff3_unknown_strand_with_real_gffread(tmp_path):
    """Integration: real gffread accepts input after '?' strand normalisation."""
    reference = _write(tmp_path / "reference.fa", SIMPLE_REFERENCE)
    annotation = _write(tmp_path / "annotation.gff3", GFF3_INPUT_WITH_UNKNOWN_STRAND)

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    prepared_text = Path(summary["annotation"]["prepared"]).read_text(encoding="utf-8")
    assert summary["annotation"]["normalised_unknown_strand_records"] == 1
    assert "\t?\t" not in prepared_text
    assert summary["annotation"]["kept_records"] > 0


def test_unstranded_records_filtered_and_saved_separately(tmp_path):
    """Unstranded entries removed from main annotation, saved to separate file."""
    reference = _write(tmp_path / "reference.fa", SIMPLE_REFERENCE)
    annotation = _write(
        tmp_path / "annotation.gtf",
        SIMPLE_GTF
        + 'chr1\tsim\ttranscript\t5\t8\t.\t.\t.\t'
        + 'gene_id "g2"; transcript_id "t2";\n',
    )

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    prepared = Path(summary["annotation"]["prepared"]).read_text(encoding="utf-8")
    unstranded = Path(
        summary["annotation"]["unstranded_path"]).read_text(encoding="utf-8")

    assert 'transcript_id "t1"' in prepared
    assert 'transcript_id "t2"' not in prepared
    assert 'transcript_id "t2"' in unstranded
    assert summary["annotation"]["excluded_unstranded_records"] == 1
    assert any("Excluded 1 unstranded" in w for w in summary["warnings"])


def test_seqname_warnings_reflect_filtered_annotation(tmp_path):
    """Seqname overlap checks use filtered annotation, not unfiltered."""
    reference = _write(
        tmp_path / "reference.fa", SIMPLE_REFERENCE + ">chrExtra\nTTTT\n")
    annotation = _write(tmp_path / "annotation.gtf", MIXED_STRANDED_GTF)

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    # chrMissing was filtered as unstranded, so NOT in seqname warnings
    assert "chrMissing" not in summary["seqnames"]["only_in_annotation"]
    assert "chrExtra" in summary["seqnames"]["only_in_reference"]
    assert not any(
        "chrMissing" in w and "seqnames" in w.lower()
        for w in summary["warnings"]
    )


def test_invalid_gff_for_conversion_raises_error(tmp_path):
    """Malformed GFF input should fail during gffread conversion."""
    reference = _write(tmp_path / "reference.fa", SIMPLE_REFERENCE)
    annotation = _write(tmp_path / "annotation.gff3", "chr1\tsim\tgene\t1\t10\n")

    with pytest.raises(
        ValueError,
        match="(Failed to convert annotation to GTF|Prepared annotation is empty)",
    ):
        prepare_annotation_reference.prepare_annotation_reference(
            annotation,
            reference,
            tmp_path / "prepared",
        )


def test_no_seqname_overlap_raises_error(tmp_path):
    """Completely mismatched reference and annotation should abort."""
    reference = _write(tmp_path / "reference.fa", ">chrOther\nAAAA\n")
    annotation = _write(tmp_path / "annotation.gtf", SIMPLE_GTF)

    with pytest.raises(ValueError, match="No overlapping seqnames"):
        prepare_annotation_reference.prepare_annotation_reference(
            annotation,
            reference,
            tmp_path / "prepared",
        )


def test_malformed_gtf_raises_error(tmp_path):
    """GTF with wrong number of columns should abort."""
    reference = _write(tmp_path / "reference.fa", SIMPLE_REFERENCE)
    annotation = _write(
        tmp_path / "annotation.gtf",
        'chr1\tsim\ttranscript\t1\t4\t.\t+\t.\n',
    )

    with pytest.raises(ValueError, match="expected 9 GTF columns"):
        prepare_annotation_reference.prepare_annotation_reference(
            annotation,
            reference,
            tmp_path / "prepared",
        )


def test_all_records_unstranded_raises_error(tmp_path):
    """If filtering removes all records, abort before downstream tools."""
    reference = _write(tmp_path / "reference.fa", SIMPLE_REFERENCE)
    annotation = _write(
        tmp_path / "annotation.gtf",
        'chr1\tsim\ttranscript\t1\t4\t.\t.\t.\tgene_id "g1"; transcript_id "t1";\n',
    )

    with pytest.raises(ValueError, match="reference annotation"):
        prepare_annotation_reference.prepare_annotation_reference(
            annotation,
            reference,
            tmp_path / "prepared",
        )


def test_conflicting_build_hints_warn(tmp_path):
    """Mismatched genome build hints produce warning but do not fail."""
    reference = _write(
        tmp_path / "reference.fa",
        ">chr1 GRCh38 primary assembly\nAAAA\n",
    )
    annotation = _write(
        tmp_path / "annotation.gtf",
        "#!genome-build GRCh37\n" + SIMPLE_GTF,
    )

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    assert set(summary["reference_build_hints"]) == {"GRCh38"}
    assert set(summary["annotation_build_hints"]) == {"GRCh37"}
    assert any("different genome build hints" in w for w in summary["warnings"])


def test_conflicting_provider_hints_warn(tmp_path):
    """Mismatched provider hints produce warning but do not fail."""
    reference = _write(tmp_path / "ensembl.GRCh38.fa", ">chr1\nAAAA\n")
    annotation = _write(tmp_path / "refseq.GRCh38.annotation.gtf", SIMPLE_GTF)

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    assert any("different provider hints" in w for w in summary["warnings"])


def test_filename_hints_captured(tmp_path):
    """Build and provider hints extracted from filenames."""
    reference = _write(tmp_path / "reference.hg38.fa", ">chr1\nAAAA\n")
    annotation = _write(tmp_path / "gencode.GRCh37.annotation.gtf", SIMPLE_GTF)

    summary = prepare_annotation_reference.prepare_annotation_reference(
        annotation,
        reference,
        tmp_path / "prepared",
    )

    assert "GRCh38" in summary["reference_build_hints"]
    assert "GRCh37" in summary["annotation_build_hints"]
    assert "GENCODE" in summary["annotation_provider_hints"]


def test_cli_component_discoverable():
    """workflow-glue discovers prepare_annotation_reference command."""
    components = get_components(allowed_components=["prepare_annotation_reference"])
    assert "prepare_annotation_reference" in components
