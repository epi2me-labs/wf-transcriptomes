#!/usr/bin/env python
"""Plot a seqkit alignment stats file."""

import argparse
import warnings

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import six

matplotlib.use('Agg')

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns

warnings.resetwarnings()
_ = sns

# Parse command line arguments:
parser = argparse.ArgumentParser(
    description="""Plot a seqkit alignment stats file.""")
parser.add_argument(
    '-r', metavar='report_pdf', type=str, help="Report PDF (stats.pdf).",
    default="stats.pdf")
parser.add_argument(
    'input', metavar='input_tsv', type=str, help="Input TSV.")


class Report:
    """Class for plotting utilities on the top of matplotlib.

    Plots are saved in the specified file through the PDF backend.
    """

    def __init__(self, pdf):
        """
        Init Report with a matplotlib PdfPages instance.

        :param self: object.
        :param pdf: Output pdf.
        :returns: The report object.
        :rtype: Report
        """
        self.pdf = pdf
        self.plt = plt
        self.pages = PdfPages(pdf)

    def _set_properties_and_close(self, fig, title, xlab, ylab):
        """Set title, axis labels and close the figure.

        :param self: object.
        :param fig: The current figure.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :returns: None
        :rtype: object
        """
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        plt.title(title)
        self.pages.savefig(fig)
        plt.close(fig)

    def plot_boxplots(self, data_map, title="", xlab="", ylab="",
                      xticks_rotation=0, xticks_fontsize=5):
        """Plot multiple pairs of data arrays.

        :param self: object.
        :param data_map:
            A dictionary with labels as keys and lists as data values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param xticks_rotation: Rotation value for x tick labels.
        :param xticks_fontsize: Fontsize for x tick labels.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()
        plt.boxplot(list(data_map.values()))
        plt.xticks(np.arange(len(data_map)) + 1, data_map.keys(),
                   rotation=xticks_rotation, fontsize=xticks_fontsize)
        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_bars_simple(self, data_map, title="", xlab="", ylab="", alpha=0.6,
                         xticks_rotation=0, auto_limit=False):
        """Plot simple bar chart from input dictionary.

        :param self: object.
        :param data_map: A dictionary with labels as keys and data as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param alpha: Alpha value.
        :param xticks_rotation: Rotation value for x tick labels.
        :param auto_limit: Set y axis limits automatically.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        labels = list(data_map.keys())
        data = list(data_map.values())
        positions = np.arange(len(labels))
        plt.bar(positions, data, align='center', alpha=alpha)
        plt.xticks(positions, labels, rotation=xticks_rotation)

        if auto_limit:
            low, high = min(data), max(data)
            plt.ylim([(low - 0.5 * (high - low)), (high + 0.5 * (high - low))])

        self._set_properties_and_close(fig, title, xlab, ylab)

    def plot_histograms(self, data_map, title="", xlab="", ylab="", bins=50,
                        alpha=0.7, legend_loc='best', legend=True,
                        vlines=None):
        """Plot histograms of multiple data arrays.

        :param self: object.
        :param data_map:
            A dictionary with labels as keys and data arrays as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param bins: Number of bins.
        :param alpha: Transparency value for histograms.
        :param legend_loc: Location of legend.
        :param legend: Plot legend if True.
        :param vlines:
            Dictionary with labels and positions of vertical lines to draw.
        :returns: None
        :rtype: object
        """
        fig = plt.figure()

        for label, data in six.iteritems(data_map):
            if len(data) > 0:
                plt.hist(data, bins=bins, label=label, alpha=alpha)
        if vlines is not None:
            for label, pos in six.iteritems(vlines):
                plt.axvline(x=pos, label=label)
        if legend:
            plt.legend(loc=legend_loc)

        self._set_properties_and_close(fig, title, xlab, ylab)

    def close(self):
        """Close PDF backend.

        Do not forget to call this at the end of your
        script or your output will be damaged!

        :param self: object
        :returns: None
        :rtype: object
        """
        self.pages.close()


if __name__ == '__main__':
    args = parser.parse_args()

    plotter = Report(args.r)

    # Plot overview panel:
    stats = pd.read_csv(args.input, sep="\t")
    perc = stats[["PrimAlnPerc", "MultimapPerc"]].copy()
    num = stats[["PrimAln", "SecAln", "SupAln", "Unmapped", "TotalReads",
                "TotalRecords"]].copy()

    perc.plot(kind='bar')
    plt.title("Percent primary and multimapping reads")
    plt.tight_layout()
    plotter.pages.savefig()

    num.plot(kind='bar')
    plt.title("Number of alignment records")
    plt.tight_layout()
    plotter.pages.savefig()

    plotter.close()
