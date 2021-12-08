#!/usr/bin/env python

"""Plot a gffcompare stats file."""

import argparse
from collections import OrderedDict
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
    description="""Plot a gffcompare stats file.""")
parser.add_argument(
    '-r', metavar='report_pdf', type=str,
    help="Report PDF (plot_gffcmp_stats.pdf).",
    default="plot_gffcmp_stats.pdf")
parser.add_argument(
    '-t', metavar='tracking_tsv', type=str,
    help="Tracking file produced by gffcompare.", default=None)
parser.add_argument(
    'input', metavar='input_txt', type=str,
    help="Input gffcompare stats file.")


class Report:
    """Matplotlib plotting utilities."""

    def __init__(self, pdf):
        """Init class with PdfPahges instance.

        Plots are saved in the specified file through the PDF backend.

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
        :param data_map: A dictionary with labels as keys and lists as data
            values.
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
        :param data_map: A dictionary with labels as keys and data arrays
            as values.
        :param title: Figure title.
        :param xlab: X axis label.
        :param ylab: Y axis label.
        :param bins: Number of bins.
        :param alpha: Transparency value for histograms.
        :param legend_loc: Location of legend.
        :param legend: Plot legend if True.
        :param vlines: Dictionary with labels and positions of vertical lines
            to draw.
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


def _parse_stat_line(sl):
    """Parse a stats line."""
    res = {}
    tmp = sl.split(':')[1]
    tmp = tmp.split('|')
    res['sensitivity'] = float(tmp[0].strip())
    res['precision'] = float(tmp[1].strip())
    return res


def _parse_matching_line(line):
    """Parse a metching line."""
    tmp = line.split(':')[1].strip()
    return int(tmp)


def _parse_mn_line(line):
    """Parse a miss or novel line."""
    res = {}
    tmp = line.split(':')[1].strip()
    tmp = tmp.split('/')
    res['value'] = int(tmp[0])
    tmp = tmp[1].split('(')
    res['value_total'] = int(tmp[0].strip())
    res['percent'] = float(tmp[1].split('%)')[0])
    return res


def _parse_total_line(line):
    """Parse a total line."""
    res = {}
    tmp = line.split(':')[1].strip()
    tmp = tmp.split('in')
    res['transcripts'] = int(tmp[0].strip())
    tmp = tmp[1].split('loci')
    res['loci'] = int(tmp[0].strip())
    tmp = int(tmp[1].split('(')[1].split(' ')[0])
    res['me_transcripts'] = tmp
    return res


def parse_gffcmp_stats(txt):
    """Parse a gffcompare stats file.

    :param txt: Path to the gffcompare stats file.
    :returns: Return as tuple of dataframes containing:
        perfromance statistics, match statistics, miss statistics,
        novel statistics, total statistics.
    :rtype: tuple
    """
    sensitivity = []
    precision = []
    level = []

    matching = OrderedDict()

    missed_level = []
    missed = []
    missed_total = []
    missed_percent = []

    novel_level = []
    novel = []
    novel_total = []
    novel_percent = []

    total_target = []
    total_loci = []
    total_transcripts = []
    total_multiexonic = []

    fh = open(txt, 'r')
    for line in fh:
        line = line.strip()
        if len(line) == 0:
            continue
        # Parse totals:
        if line.startswith('#     Query mRNAs'):
            total_target.append('Query')
            r = _parse_total_line(line)
            total_loci.append(r['loci'])
            total_transcripts.append(r['transcripts'])
            total_multiexonic.append(r['me_transcripts'])

        if line.startswith('# Reference mRNAs '):
            total_target.append('Reference')
            r = _parse_total_line(line)
            total_loci.append(r['loci'])
            total_transcripts.append(r['transcripts'])
            total_multiexonic.append(r['me_transcripts'])

        # Parse basic statistics:
        if line.startswith('Base level'):
            st = _parse_stat_line(line)
            level.append('Base')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Exon level'):
            st = _parse_stat_line(line)
            level.append('Exon')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Intron level'):
            st = _parse_stat_line(line)
            level.append('Intron')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Intron chain level'):
            st = _parse_stat_line(line)
            level.append('Intron chain')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Transcript level'):
            st = _parse_stat_line(line)
            level.append('Transcript')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])
        if line.startswith('Locus level'):
            st = _parse_stat_line(line)
            level.append('Locus')
            sensitivity.append(st['sensitivity'])
            precision.append(st['precision'])

        # Parse match statistics:
        if line.startswith('Matching intron chains'):
            m = _parse_matching_line(line)
            matching['Intron chains'] = [m]
        if line.startswith('Matching transcripts'):
            m = _parse_matching_line(line)
            matching['Transcripts'] = [m]
        if line.startswith('Matching loci'):
            m = _parse_matching_line(line)
            matching['Loci'] = [m]

        # Parse missing statistics:
        if line.startswith('Missed exons'):
            missed_level.append('Exons')
            r = _parse_mn_line(line)
            missed.append(r['value'])
            missed_total.append(r['value_total'])
            missed_percent.append(r['percent'])
        if line.startswith('Missed introns'):
            missed_level.append('Introns')
            r = _parse_mn_line(line)
            missed.append(r['value'])
            missed_total.append(r['value_total'])
            missed_percent.append(r['percent'])
        if line.startswith('Missed loci'):
            missed_level.append('Loci')
            r = _parse_mn_line(line)
            missed.append(r['value'])
            missed_total.append(r['value_total'])
            missed_percent.append(r['percent'])

        # Parse novel statistics:
        if line.startswith('Novel exons'):
            novel_level.append('Exons')
            r = _parse_mn_line(line)
            novel.append(r['value'])
            novel_total.append(r['value_total'])
            novel_percent.append(r['percent'])
        if line.startswith('Novel introns'):
            novel_level.append('Introns')
            r = _parse_mn_line(line)
            novel.append(r['value'])
            novel_total.append(r['value_total'])
            novel_percent.append(r['percent'])
        if line.startswith('Novel loci'):
            novel_level.append('Loci')
            r = _parse_mn_line(line)
            novel.append(r['value'])
            novel_total.append(r['value_total'])
            novel_percent.append(r['percent'])

    fh.close()

    df_stats = pd.DataFrame(OrderedDict(
        [('Sensitivity', sensitivity), ('Precision', precision)]), index=level)
    df_match = pd.DataFrame(matching, index=['Matching'])

    df_miss = pd.DataFrame(
        OrderedDict(
                [('Total', missed_total),
                 ('Missed', missed),
                 ('Percent missed', missed_percent)]), index=missed_level)

    df_novel = pd.DataFrame(
        OrderedDict(
                [('Total', novel_total),
                 ('Novel', novel),
                 ('Percent novel', novel_percent)]), index=novel_level)

    df_total = pd.DataFrame(OrderedDict(
        [('Loci', total_loci), ('Transcripts', total_transcripts),
         ('Multiexonic', total_multiexonic)]), index=total_target)

    return df_stats, df_match, df_miss, df_novel, df_total


if __name__ == '__main__':
    args = parser.parse_args()

    stats, match, miss, novel, total = parse_gffcmp_stats(args.input)
    tracking = pd.read_csv(args.t, sep="\t", header=None, usecols=[0, 3],
                           names=['Count', 'Overlaps'])

    tracking = tracking.groupby("Overlaps").count().reset_index()
    tracking = tracking.sort_values("Overlaps")

    plotter = Report(args.r)

    # Plot overview panel:

    plt.figure(1)
    plt.subplot(2, 2, 1)
    total.plot(ax=plt.gca(), kind='barh', sharex=False, title='Totals')
    plt.tight_layout()

    plt.subplot(2, 2, 2)
    stats.plot(ax=plt.gca(), kind='barh', legend=True, sharex=False,
               title='Performance').legend(loc='best')
    plt.tight_layout()

    plt.subplot(2, 2, 3)

    miss.copy().drop(
        'Percent missed', axis=1).plot(
            ax=plt.gca(), kind='barh',
            legend=True, sharex=False,
            title='Missed')

    plt.tight_layout()

    plt.subplot(2, 2, 4)

    novel.copy().drop(
        'Percent novel', axis=1).plot(
            ax=plt.gca(), kind='barh', legend=True, sharex=False,
            title='Novel')

    plt.tight_layout()
    plotter.pages.savefig()

    # Plot individual panels:

    total.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    stats.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    match.plot(kind='barh', subplots=True, legend=False)
    plt.tight_layout()
    plotter.pages.savefig()

    miss.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    novel.plot(kind='barh', subplots=True, legend=False, sharex=False)
    plt.tight_layout()
    plotter.pages.savefig()

    def fix_names(s):
        """Map trancript classification codes."""
        names = {
            '=': 'ExactMatch:=',
            'c': 'Contained:c',
            'k': 'ReverseContained:k',
            'm': 'RetainedIntron:m',
            'n': 'PartRetainedIntron:n',
            'j': 'PartialMatch:j',
            'e': 'TransFragMatch:e',
            's': 'OppositeMatch:s',
            'o': 'OtherSameStrand:o',
            'x': 'ExonicOpposite:o',
            'y': 'RefInIntrons:y',
            'p': 'PolymeraseRunon:p',
            'r': 'Repeat:r',
            'u': 'Intergenic:u',
            'i': 'FullyIntronic:i',
            }
        return names[s]

    # Plot overlaps panel:
    tracking.Overlaps = tracking.Overlaps.apply(fix_names)
    tracking = tracking.set_index("Overlaps")
    tracking.plot(kind='bar', title="Overlaps detected by gffcompare",
                  colormap='Paired')
    plt.tight_layout()
    plotter.pages.savefig()
    tracking["Percent"] = tracking.Count * 100 / tracking.Count.sum()
    tracking[["Percent"]].plot(
        kind='bar', title="Overlaps detected by gffcompare", colormap='Paired')
    plt.tight_layout()
    plotter.pages.savefig()

    plotter.close()
