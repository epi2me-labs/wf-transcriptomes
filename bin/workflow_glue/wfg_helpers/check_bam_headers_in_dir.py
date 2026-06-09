"""Check (u)BAM files for `@SQ` lines whether they are the same in all headers."""
import os
from pathlib import Path
import sys

import pysam

from ..util import get_named_logger, wf_parser  # noqa: ABS101


READ_SCAN_LIMIT = 10000


def alignment_has_cigar_op(alignment, cigar_op):
    """Return whether an alignment has a splice/reference-skip CIGAR op."""
    return any(op == cigar_op for op, _ in (alignment.cigartuples or []))


def alignment_has_tags(alignment, tags):
    """Return whether an alignment carries all tags in a list of SAM tags."""
    alignment_tags = {tag.upper() for tag, _ in alignment.get_tags()}
    return all(tag in alignment_tags for tag in tags)


def check_n_reads(
    alignment_file,
    check_modbase_tags=False,
    check_splice_cigars=False,
    read_scan_limit=READ_SCAN_LIMIT,
):
    """Extract bounded read-level facts from an open XAM file."""
    read_facts = {
        "has_reads": False,
    }
    if check_splice_cigars:
        read_facts["has_splice_cigars"] = False
    if check_modbase_tags:
        read_facts["has_modbase_tags"] = False

    for i, alignment in enumerate(alignment_file.fetch(until_eof=True)):
        if i >= read_scan_limit:
            break
        read_facts["has_reads"] = True
        if check_splice_cigars and not read_facts["has_splice_cigars"]:
            read_facts["has_splice_cigars"] = alignment_has_cigar_op(
                alignment, cigar_op=3)
        if check_modbase_tags and not read_facts["has_modbase_tags"]:
            read_facts["has_modbase_tags"] = alignment_has_tags(
                alignment, tags=["MM", "ML"])

        if all(read_facts.values()):
            break

    return read_facts


def check_header(alignment_file, check_ref=False):
    """Extract header information from a BAM/CRAM file."""
    # Extract SQ lines, comparing only SN/LN/M5 elements
    # (see CW-4842 - ignore different SQ.UR values)
    sq_lines = [{
        "SN": sq["SN"],
        "LN": sq["LN"],
        "M5": sq.get("M5"),
    } for sq in alignment_file.header.get("SQ", [])]

    hd_lines = alignment_file.header.get("HD")

    xam_reflen = None
    if check_ref:
        xam_reflen = set(zip(alignment_file.references, alignment_file.lengths))

    return sq_lines, hd_lines, xam_reflen


def compare_ref_lengths(xam_reflen, ref_reflen, logger):
    """Compare reference lengths from FASTA and BAM/CRAM SQ lines."""
    diff = ref_reflen.symmetric_difference(xam_reflen)
    if len(diff) > 0:
        rows = ["sequence_name sequence_length in_ref in_xam"]
        rows += [
            f"{reflen[0]} {reflen[1]} "
            f"{'1' if reflen in ref_reflen else '0'} "
            f"{'1' if reflen in xam_reflen else '0'}"
            for reflen in diff
        ]
        logger.info("\n".join(rows) + "\n")

        logger.info(
            " [WARN] In the input XAM there is at least one (name, length)"
            " pair that does not map 1:1 between the input alignment"
            " and reference. XAM will be realigned.\n"
        )
        return False
    else:
        logger.info(
            "[OKAY] Input alignment and reference sequences map 1:1\n"
        )
        return True


def main(args):
    """Run the entry point."""
    logger = get_named_logger("checkBamHdr")

    # Check the reference if given.
    ref_reflen = None
    if args.ref:
        try:
            ref = pysam.FastaFile(
                    args.ref,
                    filepath_index=args.ref_idx)
        except (ValueError, IOError):
            logger.info(
                f"[FAIL] {args.ref} reference input file could not be"
                " read. Is it in the right format?"
            )
            sys.exit(os.EX_NOINPUT)

        # get ref length for comparison with BAM/CRAM SQ lines
        ref_reflen = set(zip(ref.references, ref.lengths))

    if not args.input_path.is_dir():
        raise ValueError(f"Input path '{args.input_path}' must be a directory.")

    target_files = sorted(list(args.input_path.glob("*")))
    if not target_files:
        raise ValueError(f"No files found in input directory '{args.input_path}'.")
    # Loop over target files and check if there are `@SQ` lines in all headers or not.
    # Set `is_unaligned` accordingly.
    # Detect mixed headers by comparing @SQ lines between files.
    # Set mixed_sq_headers to True if they differ across files
    # (either some files lack @SQ lines, or files have different @SQ lines).
    # Also check if there is the SO line, to validate whether the file is (un)sorted.
    first_sq_lines = None
    mixed_sq_headers = False
    sorted_xam = False
    xam_reflen = None
    requires_realign = False
    any_sq_lines = False  # Track if any file has @SQ lines
    has_reads = False
    # a bit hacky to check these specific XAM properties here but this is a very good
    # place to do it. in future i'd like to abstract this out to be more like how sample
    # sheets are checked, with delegation to the workflow scope for determining checks
    has_splice_cigars = False
    has_modbase_tags = False

    for xam_file in target_files:
        try:
            xam_fh = pysam.AlignmentFile(xam_file, check_sq=False)
        except (ValueError, IOError):
            # File couldn't be opened
            logger.error(f"Failed to open {xam_file}")
            continue

        with xam_fh:
            sq_lines, hd_lines, xam_reflen = check_header(
                xam_fh,
                check_ref=bool(args.ref),
            )
            if not (has_reads and has_splice_cigars and has_modbase_tags):
                read_facts = check_n_reads(
                    xam_fh,
                    check_splice_cigars=not has_splice_cigars,
                    check_modbase_tags=not has_modbase_tags,
                )
                has_reads = has_reads or read_facts["has_reads"]
                has_splice_cigars = (
                    has_splice_cigars or read_facts["has_splice_cigars"])
                has_modbase_tags = (
                    has_modbase_tags or read_facts["has_modbase_tags"])

        if sq_lines:
            any_sq_lines = True

        # Check if it is sorted.
        # When there is more than one BAM, merging/sorting
        # will happen regardless of this flag.
        if hd_lines is not None and hd_lines.get('SO') == 'coordinate':
            sorted_xam = True

        if first_sq_lines is None:
            # this is the first file
            first_sq_lines = sq_lines
        else:
            # this is a subsequent file; check with the first `@SQ` lines
            if sq_lines != first_sq_lines:
                mixed_sq_headers = True
                break

        if ref_reflen and not mixed_sq_headers and not requires_realign:
            if sq_lines:
                # if a reference was given and there are sq lines (indicating aligned)
                # compare ref and xam header lengths to check alignment is to this ref
                requires_realign = not compare_ref_lengths(
                    xam_reflen, ref_reflen, logger)
            else:
                requires_realign = True

    # Assumption: presence of @SQ lines indicates aligned data.
    # Set `is_unaligned` to `True` if there were no mixed headers and none of the
    # files had `@SQ` lines, or if found to be incorrectly aligned based on
    # SQ/ref length comparison.
    is_unaligned = (not mixed_sq_headers and not any_sq_lines) or requires_realign

    # write `is_unaligned`, `is_sorted` and `mixed_sq_headers`
    # out so that they can be set as env.
    # variables and handled downstream.
    sys.stdout.write(
        f"IS_UNALIGNED={int(is_unaligned)};" +
        f"MIXED_SQ_HEADERS={int(mixed_sq_headers)};" +
        f"IS_SORTED={int(sorted_xam)};" +
        f"HAS_READS={int(has_reads)};" +
        f"HAS_SPLICE_CIGARS={int(has_splice_cigars)};" +
        f"HAS_MODBASE_TAGS={int(has_modbase_tags)};"
    )
    logger.info(f"Checked (u)BAM heads in '{args.input_path}'.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("check_bam_headers_in_dir")
    parser.add_argument("input_path", type=Path, help="Path to target directory")
    parser.add_argument(
        "--ref", type=Path, help="Optional reference FASTA file", required=False)
    parser.add_argument(
        "--ref_idx", type=Path,
        help="Optional reference index FASTA file required for pysam", required=False)
    return parser
