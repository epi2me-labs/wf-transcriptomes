#!/usr/bin/env python
"""Create workflow report."""

import argparse

from aplanat import annot, hist, report
from bokeh.layouts import gridplot
import numpy as np
import pandas as pd


def read_files(summaries):
    """Combine a list of files into a single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep="\t"))
    return pd.concat(dfs)


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    args = parser.parse_args()

    report_doc = report.HTMLReport(
        "Workflow Template Sequencing report",
        ("Results generated through the wf-template nextflow "
            "workflow by Oxford Nanopore Technologies"))

    report_doc.markdown('''
### Read Quality control
This section displays basic QC metrics indicating read data quality.
''')

    np_blue = '#0084A9'

    # read length summary
    seq_summary = read_files(args.summaries)
    total_bases = seq_summary['sequence_length_template'].sum()
    mean_length = total_bases / len(seq_summary)
    median_length = np.median(seq_summary['sequence_length_template'])
    datas = [seq_summary['sequence_length_template']]
    length_hist = hist.histogram(
        datas, colors=[np_blue], bins=100,
        title="Read length distribution.",
        x_axis_label='Read Length / bases',
        y_axis_label='Number of reads',
        xlim=(0, 2000))
    length_hist = annot.subtitle(
        length_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_length, median_length))

    datas = [seq_summary['mean_qscore_template']]
    mean_q, median_q = np.mean(datas[0]), np.median(datas[0])
    q_hist = hist.histogram(
        datas, colors=[np_blue], bins=100,
        title="Read quality score",
        x_axis_label="Quality score",
        y_axis_label="Number of reads",
        xlim=(4, 25))
    q_hist = annot.subtitle(
        q_hist,
        "Mean: {:.0f}. Median: {:.0f}".format(
            mean_q, median_q))

    report_doc.plot(gridplot([[length_hist, q_hist]]))

    # Footer section
    report_doc.markdown('''
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
    report_doc.write(args.report)


if __name__ == "__main__":
    main()
