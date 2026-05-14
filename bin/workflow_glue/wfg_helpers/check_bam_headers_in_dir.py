"""Check (u)BAM files for `@SQ` lines whether they are the same in all headers."""
import os
from pathlib import Path
import sys

import pysam

from ..util import get_named_logger, wf_parser  # noqa: ABS101


def extract_header_info(xam_file, has_reads, check_ref=False):
    """Extract header information from a BAM/CRAM file."""
    try:
        f = pysam.AlignmentFile(xam_file, check_sq=False)
    except (ValueError, IOError):
        return None, None, None, False

    with f:
        # Extract SQ lines, comparing only SN/LN/M5 elements
        # (see CW-4842 - ignore different SQ.UR values)
        sq_lines = [{
            "SN": sq["SN"],
            "LN": sq["LN"],
            "M5": sq.get("M5"),
        } for sq in f.header.get("SQ", [])]

        hd_lines = f.header.get("HD")

        xam_reflen = None
        if check_ref:
            xam_reflen = set(zip(f.references, f.lengths))
        # Whilst we have file open check for atleast 1 reads
        # unless reads already found
        if not has_reads:
            try:
                next(f.fetch(until_eof=True))
                has_reads = True
            except StopIteration:
                has_reads = False

    return sq_lines, hd_lines, xam_reflen, has_reads


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

    for xam_file in target_files:
        sq_lines, hd_lines, xam_reflen, has_reads = extract_header_info(
            xam_file, has_reads, check_ref=bool(args.ref))

        if sq_lines is None:
            # File couldn't be opened
            logger.error(f"Failed to open {xam_file}")
            continue

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
        f"HAS_READS={int(has_reads)};"
    )
    logger.info(f"Checked (u)BAM headers in '{args.input_path}'.")


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
