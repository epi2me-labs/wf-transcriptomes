"""Prepare reference and annotation files for transcriptome analysis."""

import gzip
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from .util import get_named_logger, wf_parser  # noqa: ABS101


UNSTRANDED_WARNING = (
    "Warning: Unstranded entries found and excluded from downstream "
    "transcriptome analysis."
)

BUILD_PATTERNS = {
    "GRCh38": re.compile(r"(?<![A-Za-z0-9])(?:GRCh38|hg38)(?![A-Za-z0-9])", re.I),
    "GRCh37": re.compile(r"(?<![A-Za-z0-9])(?:GRCh37|hg19)(?![A-Za-z0-9])", re.I),
    "GRCm39": re.compile(r"(?<![A-Za-z0-9])GRCm39(?![A-Za-z0-9])", re.I),
    "GRCm38": re.compile(r"(?<![A-Za-z0-9])GRCm38(?![A-Za-z0-9])", re.I),
}

PROVIDER_PATTERNS = {
    "GENCODE": re.compile(r"(?<![A-Za-z0-9])gencode(?![A-Za-z0-9])", re.I),
    "Ensembl": re.compile(r"(?<![A-Za-z0-9])ensembl(?![A-Za-z0-9])", re.I),
    "RefSeq": re.compile(r"(?<![A-Za-z0-9])(?:refseq|ncbi)(?![A-Za-z0-9])", re.I),
}


def _is_gzip(path):
    """Return whether path looks gzip-compressed."""
    return Path(path).suffix.lower() == ".gz"


def _is_gff(path):
    """Return whether path looks like a GFF/GFF3 annotation."""
    p = Path(path)
    # check all suffixes for GFF extensions
    suffixes = {s.lower() for s in p.suffixes}
    return bool(suffixes & {".gff", ".gff3"})


def _find_patterns(text, patterns):
    """Find matching pattern names in text."""
    return {name for name, pattern in patterns.items() if pattern.search(text)}


class PreparedReference:
    """A reference genome prepared for analysis."""

    def __init__(self, input_path, work_dir):
        """Initialize and prepare reference genome."""
        self.input_path = Path(input_path)
        self.work_dir = Path(work_dir)
        self.output_path = self.work_dir / "reference.fasta"
        self._seqnames = None
        self._prepare()

    def _prepare(self):
        """Prepare reference (decompress if gzipped)."""
        if _is_gzip(self.input_path):
            with gzip.open(self.input_path, "rb") as src, open(
                self.output_path, "wb"
            ) as dst:
                shutil.copyfileobj(src, dst)
        else:
            shutil.copyfile(self.input_path, self.output_path)

        if not self.output_path.exists() or self.output_path.stat().st_size < 1:
            raise ValueError(f"Prepared file is empty: {self.output_path}")

    @property
    def seqnames(self):
        """Extract and cache sorted FASTA sequence names."""
        if self._seqnames is None:
            ids = set()
            with open(self.output_path, encoding="utf-8") as f:
                for line in f:
                    if line.startswith(">"):
                        ids.add(line[1:].split()[0])
            self._seqnames = sorted(ids)
        return self._seqnames

    def detect_hints(self):
        """Detect build and provider hints from filename and headers."""
        builds = set()
        providers = set()

        # check filename
        builds.update(_find_patterns(self.input_path.name, BUILD_PATTERNS))
        providers.update(_find_patterns(self.input_path.name, PROVIDER_PATTERNS))

        # check headers
        with open(self.output_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith(">"):
                    builds.update(_find_patterns(line, BUILD_PATTERNS))
                    providers.update(_find_patterns(line, PROVIDER_PATTERNS))

        return {
            "builds": sorted(builds),
            "providers": sorted(providers),
        }


class AttributeSanitizer:
    """Fixes common GTF attribute problems."""

    ATTRIBUTE_RE = re.compile(r'(^|;\s*){key}\s+"([^"]*)";?')

    @classmethod
    def _pattern_for(cls, key):
        """Return regex pattern for a GTF attribute key."""
        return re.compile(cls.ATTRIBUTE_RE.pattern.format(key=re.escape(key)))

    @classmethod
    def _extract(cls, attr_field, key):
        """Extract a quoted GTF attribute value."""
        match = cls._pattern_for(key).search(attr_field)
        if not match:
            return None
        value = match.group(2)
        # normalize whitespace and remove quotes/semicolons
        value = re.sub(r"\s+", " ", value.replace('"', "").replace(";", "")).strip()
        return value or None

    @classmethod
    def _replace(cls, attr_field, key, value):
        """Replace a quoted GTF attribute value."""
        pattern = cls._pattern_for(key)
        return pattern.sub(
            lambda m: f'{m.group(1)}{key} "{value}"; ',  # noqa: E702
            attr_field,
            count=1,
        )

    @classmethod
    def fix(cls, attr_field):
        """Fix gene_id if missing or set to literal 'transcript_id'."""
        transcript_id = cls._extract(attr_field, "transcript_id")
        if transcript_id is None:
            return attr_field, False

        gene_id = cls._extract(attr_field, "gene_id")

        # add gene_id if missing
        if gene_id is None:
            return (
                f'gene_id "{transcript_id}"; {attr_field.strip()}',  # noqa: E702
                True,
            )

        # fix gene_id="transcript_id" literal
        if gene_id == "transcript_id":
            return cls._replace(attr_field, "gene_id", transcript_id), True

        return attr_field, False


class Annotation:
    """An annotation file with preparation and filtering capabilities."""

    MAX_UNSTRANDED_EXAMPLES = 20

    def __init__(self, input_path, work_dir, gffread="gffread"):
        """Initialize and convert annotation to GTF."""
        self.input_path = Path(input_path)
        self.work_dir = Path(work_dir)
        self.is_gff = _is_gff(input_path)
        self.gffread = gffread

        # outputs
        self.unfiltered_path = None
        self.prepared_path = None
        self._convert_to_gtf()

    def _open_gtf(self, path):
        """Open a GTF file, handling gzip transparently."""
        if _is_gzip(path):
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, encoding="utf-8")

    def _read_headers(self, path):
        """Yield comment/header lines from a GTF file."""
        with self._open_gtf(path) as f:
            for line in f:
                if line.startswith("#"):
                    yield line.rstrip("\n")

    def _read_records(self, path):
        """Yield (line_text, fields_list) for all lines in a GTF file."""
        with self._open_gtf(path) as f:
            for line in f:
                stripped = line.rstrip("\n")
                if not stripped or stripped.startswith("#"):
                    yield line, None
                else:
                    fields = stripped.split("\t")
                    if len(fields) != 9:
                        raise ValueError(
                            f"Malformed annotation record has "
                            f"{len(fields)} columns; expected 9 GTF "  # noqa: E702
                            f"columns: {stripped}"
                        )
                    yield line, fields

    def _extract_seqnames(self, path):
        """Extract sorted unique seqnames from a GTF file."""
        ids = set()
        with self._open_gtf(path) as f:
            for line in f:
                stripped = line.rstrip("\n")
                if stripped and not stripped.startswith("#"):
                    fields = stripped.split("\t")
                    if len(fields) >= 1:
                        ids.add(fields[0])
        return sorted(ids)

    def _convert_to_gtf(self):
        """Convert GFF/GFF3 to GTF using gffread, or use input if already GTF."""
        if self.is_gff:
            dest = self.work_dir / "annotation_unfiltered.gtf"

            if _is_gzip(self.input_path):
                # decompress as gff doesn't support gzip input
                name_without_gz = self.input_path.name[:-3]
                suffix = Path(name_without_gz).suffix or ".gff"
                intermediate = self.work_dir / f"annotation_input{suffix}"
                with gzip.open(self.input_path, "rb") as src, open(
                    intermediate, "wb"
                ) as dst:
                    shutil.copyfileobj(src, dst)
            else:
                intermediate = self.input_path

            # use gffread to convert to GTF
            result = subprocess.run(
                [self.gffread, "-T", str(intermediate), "-o", str(dest)],
                check=False,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                details = "\n".join(
                    part
                    for part in (result.stdout, result.stderr)
                    if part
                )
                raise ValueError(
                    "Failed to convert annotation to GTF with "
                    "gffread:\n" + details
                )

            if not dest.exists() or dest.stat().st_size < 1:
                raise ValueError(f"Prepared annotation is empty: {dest}")

            self.unfiltered_path = dest
        else:
            # already GTF
            if not self.input_path.exists() or self.input_path.stat().st_size < 1:
                raise ValueError(f"Prepared annotation is empty: {self.input_path}")
            self.unfiltered_path = self.input_path

    @property
    def seqnames(self):
        """Get seqnames from prepared (stranded) annotation."""
        if self.prepared_path is None:
            raise ValueError("Must call filter_stranded() before accessing seqnames")
        return self._extract_seqnames(self.prepared_path)

    def filter_stranded(self):
        """Filter to stranded records, return statistics."""
        dest = self.work_dir / "annotation.gtf"
        unstranded_dest = self.work_dir / "unstranded_annotation.gtf"

        stats = {
            "total_records": 0,
            "kept_records": 0,
            "excluded_unstranded_records": 0,
            "sanitised_attribute_records": 0,
            "unstranded_examples": [],
        }

        unstranded_file = None
        try:
            with open(dest, "w", encoding="utf-8") as out:
                for line, fields in self._read_records(self.unfiltered_path):
                    # pass through headers and blank lines
                    if fields is None:
                        out.write(line)
                        continue

                    stats["total_records"] += 1

                    if fields[6] in {"+", "-"}:
                        # fix attributes if needed of stranded records
                        fixed_attrs, was_fixed = AttributeSanitizer.fix(fields[8])
                        if was_fixed:
                            stats["sanitised_attribute_records"] += 1
                            fields[8] = fixed_attrs

                        out.write("\t".join(fields) + "\n")
                        stats["kept_records"] += 1
                    else:
                        # collect unstranded records for posterity
                        stats["excluded_unstranded_records"] += 1
                        if (
                            len(stats["unstranded_examples"])
                            < self.MAX_UNSTRANDED_EXAMPLES
                        ):
                            stats["unstranded_examples"].append(
                                line.rstrip("\n")
                            )

                        if unstranded_file is None:
                            unstranded_file = open(
                                unstranded_dest, "w", encoding="utf-8"
                            )
                        unstranded_file.write(line)
        finally:
            if unstranded_file:
                unstranded_file.close()

        # clean up empty unstranded file
        if stats["excluded_unstranded_records"] == 0 and unstranded_dest.exists():
            unstranded_dest.unlink()

        self.prepared_path = dest
        stats["unstranded_path"] = (
            str(unstranded_dest) if stats["excluded_unstranded_records"] > 0 else None
        )
        return stats

    def detect_hints(self):
        """Detect build and provider hints from filename and headers."""
        builds = set()
        providers = set()

        # check filename
        builds.update(_find_patterns(self.input_path.name, BUILD_PATTERNS))
        providers.update(_find_patterns(self.input_path.name, PROVIDER_PATTERNS))

        # check headers
        for line in self._read_headers(self.input_path):
            builds.update(_find_patterns(line, BUILD_PATTERNS))
            providers.update(_find_patterns(line, PROVIDER_PATTERNS))

        return {
            "builds": sorted(builds),
            "providers": sorted(providers),
        }


def validate_seqname_overlap(annotation_seqnames, reference_seqnames):
    """Return overlap statistics for annotation/reference seqnames."""
    ann_set = set(annotation_seqnames)
    ref_set = set(reference_seqnames)
    matches = sorted(ann_set & ref_set)
    only_in_annotation = sorted(ann_set - ref_set)
    only_in_reference = sorted(ref_set - ann_set)

    return {
        "has_overlap": bool(matches),
        "matches": matches,
        "only_in_annotation": only_in_annotation,
        "only_in_reference": only_in_reference,
    }


def build_warnings(filter_stats, seqnames, ref_hints, ann_hints):
    """Build all warning messages."""
    warnings = []

    if filter_stats["excluded_unstranded_records"]:
        warnings.extend([
            UNSTRANDED_WARNING,
            (
                f"Excluded {filter_stats['excluded_unstranded_records']} "
                f"unstranded annotation records; "  # noqa: E702
                f"kept {filter_stats['kept_records']} stranded records."
            ),
            "A sample of unstranded entries:",
            *filter_stats["unstranded_examples"],
        ])

    if filter_stats["sanitised_attribute_records"]:
        warnings.append(
            f"Warning: Sanitised {filter_stats['sanitised_attribute_records']} "
            "annotation records with missing or malformed gene_id attributes."
        )

    if seqnames["only_in_annotation"]:
        warnings.extend([
            "Warning: Some seqnames are present in the annotation but not the genome:",
            *seqnames["only_in_annotation"][:5],
        ])
    if seqnames["only_in_reference"]:
        warnings.extend([
            "Warning: Some seqnames are present in the genome but not the annotation:",
            *seqnames["only_in_reference"][:5],
        ])

    ref_builds = ref_hints["builds"]
    ann_builds = ann_hints["builds"]
    if (
        ref_builds
        and ann_builds
        and not (set(ref_builds) & set(ann_builds))
    ):
        warnings.append(
            "Warning: Reference and annotation contain different genome "
            f"build hints: reference={', '.join(ref_builds)}; "  # noqa: E702
            f"annotation={', '.join(ann_builds)}."
        )

    ref_providers = ref_hints["providers"]
    ann_providers = ann_hints["providers"]
    if (
        ref_providers
        and ann_providers
        and not (set(ref_providers) & set(ann_providers))
    ):
        warnings.append(
            "Warning: Reference and annotation contain different "
            f"provider hints: "  # noqa: E702
            f"reference={', '.join(ref_providers)}; "  # noqa: E702
            f"annotation={', '.join(ann_providers)}."
        )

    return warnings


def prepare_annotation_reference(annotation, reference, out_dir, gffread="gffread"):
    """Prepare annotation.gtf and reference.fasta in out_dir."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)

    # prepare reference genome
    ref = PreparedReference(reference, out_dir)

    # prepare annotation (convert to GTF if needed, then filter to stranded)
    ann = Annotation(annotation, out_dir, gffread=gffread)
    filter_stats = ann.filter_stranded()

    # check seqname overlap between prepared annotation and reference
    seqnames = validate_seqname_overlap(ann.seqnames, ref.seqnames)

    if not seqnames["has_overlap"]:
        raise ValueError(
            "ERROR: No overlapping seqnames were found between the "
            "reference annotation and the reference genome.\n"
            "Annotation ID examples:\n"
            + "\n".join(seqnames["only_in_annotation"][:5])
            + "\n"
            + "Reference ID examples:\n"
            + "\n".join(seqnames["only_in_reference"][:5])
        )

    # heuristically detect build and provider hints from filenames and headers
    ref_hints = ref.detect_hints()
    ann_hints = ann.detect_hints()
    warnings = build_warnings(filter_stats, seqnames, ref_hints, ann_hints)

    # build a summary
    summary = {
        "annotation": {
            "input": str(ann.input_path),
            "prepared": str(ann.prepared_path),
            "unfiltered": str(ann.unfiltered_path),
            "was_gff": ann.is_gff,
            **filter_stats,
        },
        "reference": {
            "input": str(ref.input_path),
            "prepared": str(ref.output_path),
        },
        "seqnames": seqnames,
        "analysis_seqnames": seqnames,
        "seqname_overlap": seqnames["matches"],
        "only_in_annotation": seqnames["only_in_annotation"],
        "only_in_reference": seqnames["only_in_reference"],
        "reference_build_hints": ref_hints["builds"],
        "annotation_build_hints": ann_hints["builds"],
        "reference_provider_hints": ref_hints["providers"],
        "annotation_provider_hints": ann_hints["providers"],
        "warnings": warnings,
    }

    summary_path = out_dir / "annotation_reference_summary.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")

    return summary


def main(args):
    """Run the annotation/reference preparation command."""
    logger = get_named_logger("prepRef")
    summary = prepare_annotation_reference(
        args.annotation,
        args.reference,
        args.out_dir,
        gffread=args.gffread,
    )

    for warning in summary["warnings"]:
        sys.stdout.write(f"{warning}\n")

    logger.info(
        "Prepared annotation %s and reference %s.",
        summary["annotation"]["prepared"],
        summary["reference"]["prepared"],
    )


def argparser():
    """Argument parser for the preparation entry point."""
    parser = wf_parser("prepare_annotation_reference")
    parser.add_argument(
        "--annotation", required=True, help="Reference annotation GTF/GFF/GFF3."
    )
    parser.add_argument("--reference", required=True, help="Reference genome FASTA.")
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory for prepared annotation/reference files.",
    )
    parser.add_argument("--gffread", default="gffread", help="gffread executable.")
    return parser
