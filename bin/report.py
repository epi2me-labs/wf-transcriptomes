#!/usr/bin/env python
"""Create workflow report."""

import argparse
from collections import OrderedDict

from aplanat import bars, lines
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource
from bokeh.palettes import Category10_10
from bokeh.plotting import figure
from bokeh.transform import dodge
import numpy as np
import pandas as pd


def simple_hbar(df, y, right, title="", color=Colors.cerulean,
                fig_kwargs={}, plot_kwargs={}):
    """Create a simple barplot.

    :param groups: the grouping variable (the x-axis values).
    :param values: the data for bars are drawn (the y-axis values).
    :param kwargs: kwargs for bokeh figure.

    Move to planat when it's working?
    """
    defaults = {
        'output_backend': 'webgl',
        'plot_height': 300, 'plot_width': 600}
    defaults.update(fig_kwargs)

    p = figure(y_range=df[y], height=250, title=title,
               toolbar_location=None, tools="")

    plot_kwargs.update({'height': 0.2})
    p.hbar(y=df[y], right=df[right], **plot_kwargs)

    return p


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


def grouped_bar(df, title=""):
    """Create grouped bar plot from pandas dataframe.

    :param pandas.DataFrame
        Index:
            str: the x group labels - groups cluserted using these
        Columns:
           numeric: sub-groups of data - each sub group has same colour

    :returns bokaoh.plotting.figure instance
    """
    min_ = 0
    max_ = df.to_numpy().max()
    max_ = max_ + (max_ * 0.3)  # Add some padding at top of plot for legends
    yrange = int(min_), int(max_)

    df['x_groups'] = df.index
    df = df.reset_index(drop=True)
    source = ColumnDataSource(data=df)

    p = figure(x_range=df['x_groups'], y_range=yrange, height=250, title=title,
               toolbar_location=None, tools="")
    i = 0
    # Use the dodge method to plot groups of bars
    # https://docs.bokeh.org/en/latest/docs/user_guide/categorical.html
    dodge_range = (-0.25, 0.25)
    current_dodge = dodge_range[0]
    dodge_increment = abs(dodge_range[0] - dodge_range[1]) \
        / (len(df.columns) - 1)

    for col in df.columns:

        num_colors = df.shape[1] - 1
        colors = list(zip(*[[Category10_10[x]] * (len(df.columns) - 1)
                      for x in range(num_colors)]))
        colors = [item for sublist in colors for item in sublist]

        if col == 'x_groups':
            continue
        color = colors[i]
        i += 1

        width = df.size / 60
        p.vbar(x=dodge('x_groups', current_dodge, range=p.x_range), top=col,
               width=width, source=source, color=color, legend_label=col)
        current_dodge += dodge_increment

    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.legend.location = "top_left"
    p.legend.orientation = "horizontal"
    return p


def workflow_plots(report, df_aln_stats_file,
                   gff_cmp_stats_file, gff_cmp_tracking_file):
    """Create various sections and plots in a WfReport.

    :param report: aplanat WFReport
    :param df_aln_stats_file: alignment stats. Output of `seqkit bam -s`
    :param gff_cmp_stats_file: gffcompare stats file
    :param gff_cmp_tracking_file: gffcompare tracking file
    :return: None
    """
    df_aln_stats = pd.read_csv(df_aln_stats_file, sep='\t')
    df_aln_stats = df_aln_stats.select_dtypes([np.number]).dropna(axis=1)

    section = report.add_section()
    section.markdown('''
    ### Read mapping summary

    Output of [seqkit](https://bioinf.shenwei.me/seqkit/) bam -s''')

    section.table(df_aln_stats)

    # Percentage primary and secondary mapping
    df_perc = df_aln_stats[['PrimAlnPerc', 'MultimapPerc']]
    bar_perc = bars.simple_bar(
        df_perc.columns.values, df_perc.iloc[0].values,
        title='% primary and multimapping reads', colors=Colors.cerulean)

    # Counts of read mapping class
    df_counts = df_aln_stats.drop(columns=['PrimAlnPerc', 'MultimapPerc'])
    bar_counts = bars.simple_bar(
        df_counts.columns.values, df_counts.iloc[0].values,
        title='Number of alignment records', colors=Colors.cerulean
    )

    grid = gridplot([bar_perc, bar_counts], ncols=2,
                    plot_width=400, plot_height=400)
    section.plot(grid)

    # If gffcompare has not been run, finish report here
    if not gff_cmp_stats_file or not gff_cmp_tracking_file:
        return

    stats, _, miss, novel, total = \
        parse_gffcmp_stats(gff_cmp_stats_file)

    # Plot overview panel:
    section = report.add_section()
    section.markdown('''
    ### Annotation summary

    The following plots summarize some of the output from
    [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

    * **Totals**:
    Comparison of the number of stringtie-generated
    transcripts, multiexonic transcripts and
    loci (I'm not exactly sure what defines this class at the moment) between
    reference

    * **Performance**:
    How accurate are the query transcript annotations with respect to the
    reference at various levels.

    * **Missed**:
    Features present in the reference, but absent in the query

    * **Novel**:
    Features present in the query transcripts, but absent in the reference
    ''')

    bar_totals = grouped_bar(total, title="Totals")
    bar_performance = grouped_bar(stats, title="Performance")
    bar_missed = grouped_bar(miss, title="Missed")
    bar_novel = grouped_bar(novel, title="Novel")

    grid = gridplot([bar_totals, bar_performance, bar_missed, bar_novel],
                    ncols=2, plot_width=400, plot_height=400)
    section.plot(grid)

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
    section = report.add_section()
    section.markdown('''
    ## Query transfrag class assignments

    The classes that are assinged by
    [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml),
    which describe the relationship between query transfrag and the most
    similar reference transcript.

    [This diagram](https://ccb.jhu.edu/software/stringtie/
    gffcompare_codes.png) illustrates the different classes.
        ''')

    tracking = pd.read_csv(gff_cmp_tracking_file, sep="\t", header=None,
                           usecols=[0, 3], names=['Count', 'Overlaps'])

    tracking = tracking.groupby("Overlaps").count().reset_index()
    tracking = tracking.sort_values("Overlaps")
    tracking.Overlaps = tracking.Overlaps.apply(fix_names)
    tracking["Percent"] = tracking.Count * 100 / tracking.Count.sum()

    tracking_bar = simple_hbar(
        tracking, 'Overlaps', 'Count', title="totals")
    tracking_bar_perc = simple_hbar(
        tracking, 'Overlaps', 'Percent', title='percent'
    )

    grid = gridplot([tracking_bar, tracking_bar_perc], ncols=2,
                    plot_width=400, plot_height=400)
    section.plot(grid)


def pychopper_plots(report, df):
    """Make plots from pychopper.cdna_classifier.py.

    :param report: aplanat WFReport
    :param df: result DataFrame
    """
    section = report.add_section()
    section.markdown('''
    ### pychopper summary statisitcs

    The following plots summarize the output of [cdna_classifier.py]
    (https://github.com/nanoporetech/pychopper)

    * **Classification of output reads**:
        * Primers_found: Reads with primers found in correct orientation at
        both ends.
        * Rescue: Reads 'rescued' from fused reads
        * Unusable: Read with missing or incorrect primer orientation
    * **Strand of oriented reads**:
        * Strand of read relative to the mRNA
    * **Strand of rescued read**:
        * Strand of read that were rescued from fused reads
    * **Number of primer alignment hits in unclassified reads**:
        * Note: Need to look into what this means
    * **Number of primer alignment hits in rescued reads**:
        * Note: Need to look into what this means
    * **Number of usable segments per rescued read**:
        * Number of usable segments (primer-flanked, correctly oriented
        regions) per fused read.
    * **Usable bases as a function of cutoff**:
        * The cutoff value supplied to the primer alignment tool.
        Note: What are usabel bases in this conext
    * ** Log10 length distribution of trimmed away sequences**:
        * todo
    ''')

    def g(df, index, title):
        df_ = df[df.index == index]
        groups = df_.Name.values
        bar_ = bars.simple_bar(
            groups, df_['Value'].values,
            title=title, colors=Colors.cerulean)
        return bar_

    df1 = df.set_index('Name', drop=True)
    df1 = df1.T[['Primers_found', 'Rescue', 'Unusable']]
    bar_class = bars.simple_bar(
            df1.columns.values, df1.iloc[0].values,
            title='Classification of output reads', colors=Colors.cerulean)

    plots = [
        bar_class,
        g(df, 'Strand', 'Strand of oriented read'),
        g(df, 'RescueStrand', 'Strand of rescued reads'),
        g(df, 'UnclassHitNr', 'Number of hits in unclassified reads'),
        g(df, 'RescueHitNr', 'Number of hits in rescued reads'),
        g(df, 'RescueSegmentNr', 'Number of usable segments per rescued read')
        ]

    q = round(df.loc['Parameter', 'Value'], 4)
    df_at = df[df.index == 'AutotuneSample'].astype('float')

    # Add vertical line at x=q
    ymin, ymax = df_at['Value'].min(), df_at['Value'].max()
    plots.append(lines.line([df_at['Name'].values.tolist(), [q, q]],
                            [df_at['Value'].values.tolist(), [ymin, ymax]],
                            title=("Usable bases as function of cutoff(q).Best"
                                   " q={}").format(q), colors=['blue', 'red']
                            ))
    df_unusable = df[df.index == 'Unusable'].astype('float')

    plots.append(lines.line([np.log10(1 + df_unusable['Name'])],
                            [df_unusable['Value']],
                            title=("Log10 length distribution of trimmed away"
                                   " sequences.")
                            ))

    grid = gridplot(plots, ncols=2,
                    plot_width=400, plot_height=400)
    section.plot(grid)


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument("summaries", nargs='+', help="Read summary file.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")

    parser.add_argument(
        "--alignment_stats", required=True,
        help="TSV summary file of alignment statistics")

    parser.add_argument(
        "--gffcompare_tracking", required=False, default=None,
        help="TSV summary file of alignment statistics")
    parser.add_argument(
        "--gffcompare_stats", required=False, default=None,
        help="TSV summary file of alignment statistics")

    parser.add_argument(
        "--pychop_report", required=True,
        help="TSV summary file of pychopper statistics")
    args = parser.parse_args()

    report = WFReport(
        "Workflow for assembling transcript isoforms", "wf-isoforms",
        revision=args.revision, commit=args.commit)

    # Add reads summary section
    report.add_section(
        section=fastcat.full_report(args.summaries))

    # workflow-specific plotting
    workflow_plots(report, args.alignment_stats,  args.gffcompare_stats,
                   args.gffcompare_tracking)

    if args.pychop_report:
        df_chop_stats = pd.read_csv(args.pychop_report, sep='\t', index_col=0)
        pychopper_plots(report, df_chop_stats)

    # Arguments and software versions
    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
