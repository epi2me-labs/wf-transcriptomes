#!/usr/bin/env python
"""Create workflow report."""

import argparse

from aplanat.components import fastcat
from aplanat.report import HTMLReport


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    args = parser.parse_args()

    report = HTMLReport(
        "Workflow Template Sequencing report",
        ("Results generated through the wf-template nextflow "
            "workflow by Oxford Nanopore Technologies"))

    report.add_section(
        section=fastcat.full_report(args.summaries))

    report.markdown('''
### About

**Oxford Nanopore Technologies products are not intended for use for health
assessment or to diagnose, treat, mitigate, cure or prevent any disease or
condition.**

This report was produced using the
[epi2me-labs/wf-template](https://github.com/epi2me-labs/wf-template).  The
workflow can be run using `nextflow epi2me-labs/wf-template --help`

---
''')

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
