#!/usr/bin/env python
"""Create workflow report."""

import argparse

from aplanat.components import fastcat
from aplanat.report import WFReport
import conda_versions


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    args = parser.parse_args()

    report = WFReport(
        "Workflow Template Sequencing report", "wf-template",
        revision=args.revision, commit=args.commit)

    report.add_section(
        section=fastcat.full_report(args.summaries))

    section = report.add_section()
    section.markdown('''
### Software versions
The table below highlights versions of key software used within the analysis.
''')
    req = [
        'python', 'aplanat', 'pysam', 'fastcat']
    versions = conda_versions.scrape_data(
        as_dataframe=True, include=req)
    section.table(versions[['Name', 'Version', 'Build']], index=False)

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
