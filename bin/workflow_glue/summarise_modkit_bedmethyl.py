#!/usr/bin/env python
"""Summarise per-sample modkit bedMethyl outputs."""

from collections import defaultdict
import csv
import gzip

from .util import wf_parser  # noqa: ABS101


SUMMARY_FIELDS = (
    "sample",
    "full_mod_code",
    "mod_code",
    "mod_label",
    "valid_coverage",
    "modified_calls",
    "canonical_calls",
    "other_calls",
    "delete_calls",
    "fail_calls",
    "diff_calls",
    "nocall_calls",
    "modification_percent",
)


def load_mod_code_labels(label_tsv):
    """Load friendly labels keyed by full mod code."""
    labels = {}
    with open(label_tsv, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            mod_code = row.get("mod_code", "").strip()
            label = row.get("label", "").strip()
            if mod_code and label:
                labels[mod_code] = label
    return labels


def open_text_maybe_gzip(path):
    """Open plain text or gzipped text files transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, encoding="utf-8")


def build_full_code_lookup(mod_codes):
    """Map stripped mod codes back to their full primary_base:mod_code form."""
    lookup = {}
    for full_code in mod_codes.split(","):
        full_code = full_code.strip()
        if not full_code:
            continue
        if ":" in full_code:
            _, mod_code = full_code.split(":", 1)
        else:
            mod_code = full_code
        lookup[mod_code] = full_code
    return lookup


def summarise_bedmethyl(bedmethyl, sample, mod_codes, labels):
    """Aggregate a bedMethyl file into one row per modification code."""
    full_code_lookup = build_full_code_lookup(mod_codes)
    counts = defaultdict(
        lambda: {
            "valid_coverage": 0,
            "modified_calls": 0,
            "canonical_calls": 0,
            "other_calls": 0,
            "delete_calls": 0,
            "fail_calls": 0,
            "diff_calls": 0,
            "nocall_calls": 0,
        }
    )

    # TODO would be nice to tidy this up as an ezcharts component
    with open_text_maybe_gzip(bedmethyl) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            mod_code = fields[3]
            stats = counts[mod_code]
            stats["valid_coverage"] += int(fields[9])
            stats["modified_calls"] += int(fields[11])
            stats["canonical_calls"] += int(fields[12])
            stats["other_calls"] += int(fields[13])
            stats["delete_calls"] += int(fields[14])
            stats["fail_calls"] += int(fields[15])
            stats["diff_calls"] += int(fields[16])
            stats["nocall_calls"] += int(fields[17])

    rows = []
    for mod_code, stats in sorted(counts.items()):
        full_mod_code = full_code_lookup.get(mod_code, mod_code)
        valid_coverage = stats["valid_coverage"]
        modification_percent = (
            (stats["modified_calls"] / valid_coverage) * 100
            if valid_coverage else 0.0
        )
        rows.append(
            {
                "sample": sample,
                "full_mod_code": full_mod_code,
                "mod_code": mod_code,
                "mod_label": labels.get(full_mod_code, mod_code),
                **stats,
                "modification_percent": f"{modification_percent:.2f}",
            }
        )
    return rows


def write_summary(rows, output):
    """Write summary rows to a TSV file."""
    with open(output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main(args):
    """Run the entry point."""
    labels = load_mod_code_labels(args.mod_code_labels)
    rows = summarise_bedmethyl(
        args.bedmethyl,
        args.sample,
        args.mod_codes,
        labels,
    )
    write_summary(rows, args.output)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("summarise_modkit_bedmethyl")
    parser.add_argument("bedmethyl", help="Input modkit bedMethyl file")
    parser.add_argument("sample", help="Sample alias")
    parser.add_argument("mod_codes", help="Comma-separated full modification codes")
    parser.add_argument("mod_code_labels", help="TSV mapping full mod codes to labels")
    parser.add_argument("output", help="Output TSV")
    return parser
