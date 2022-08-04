#!/usr/bin/env python
"""Create workflow report."""

import argparse
from collections import Counter, defaultdict, OrderedDict
import math
from pathlib import Path

from aplanat import bars, hist, lines
from aplanat.components import simple as scomponents
from aplanat.components.fastcat import read_length_plot, read_quality_plot
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, Legend, Panel, Tabs
from bokeh.models.widgets import DataTable, TableColumn
from bokeh.palettes import Category10_10
from bokeh.plotting import figure
from bokeh.transform import dodge
import gffutils
import numpy as np
import pandas as pd
import sigfig


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


def grouped_bar(df, title="", tilted_xlabs=False):
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

    p = figure(
        x_range=df['x_groups'], y_range=yrange, height=250, title=title,
        toolbar_location=None, tools="")
    i = 0
    # Use the dodge method to plot groups of bars
    # https://docs.bokeh.org/en/latest/docs/user_guide/categorical.html
    dodge_range = (-0.25, 0.25)
    current_dodge = dodge_range[0]
    dodge_increment = abs(dodge_range[0] - dodge_range[1]) / \
        (len(df.columns) - 1)

    legend_it = []

    if tilted_xlabs:
        p.xaxis.major_label_orientation = math.pi / 4

    for col in df.columns:

        num_colors = df.shape[1] - 1
        colors = list(zip(
            *[[Category10_10[x]] * (len(df.columns) - 1)
                for x in range(num_colors)]))
        colors = [item for sublist in colors for item in sublist]

        if col == 'x_groups':
            continue
        color = colors[i]
        i += 1

        width = df.size / 60
        v = p.vbar(
            x=dodge('x_groups', current_dodge, range=p.x_range), top=col,
            width=width, source=source, color=color)
        current_dodge += dodge_increment
        legend_it.append([col, [v]])

    legend = Legend(items=legend_it)
    p.add_layout(legend, 'right')
    p.x_range.range_padding = 0.1
    p.xgrid.grid_line_color = None
    p.legend.location = "top_left"
    p.legend.orientation = "vertical"
    return p


def gff_compare_plots(report, gffcompare_outdirs: Path, sample_ids):
    """Create various sections and plots in a WfReport.

    :param report: aplanat WFReport
    :param gffcompare_outdirs: List of output directories from run_gffcompare
    :return: None

    TODO: split this into separate functions
    """
    # If any of the gffcompare dirs are empty, skip this section
    if not all([any(Path(x).iterdir()) for x in gffcompare_outdirs]):
        return

    # Plot overview panel:
    section = report.add_section()
    gffcompare_md = ('''
    ### Annotation summary

    The following plots summarize some of the output from
    [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)

    * **Totals**:
    Comparison of the number of stringtie-generated
    transcripts, multiexonic transcripts and loci between reference and query.

    * **Performance**:
    How accurate are the query transcript annotations with respect to the
    reference at various levels.

    * **Missed**:
    Features present in the reference, but absent in the query

    * **Novel**:
    Features present in the query transcripts, but absent in the reference
    ''')

    tabs = []
    gff_fails = False
    for id_, dir_ in zip(sample_ids, gffcompare_outdirs):
        stats, _, miss, novel, total = \
            parse_gffcmp_stats(dir_ / 'str_merged.stats')

        if not any([x.empty for x in [stats, miss, novel, total]]):
            bar_totals = grouped_bar(total, title="Totals")
            bar_performance = grouped_bar(
                stats, title="Performance", tilted_xlabs=True)
            bar_missed = grouped_bar(miss, title="Missed")
            bar_novel = grouped_bar(novel, title="Novel")
            tabs.append(Panel(
                child=gridplot(
                    [bar_totals, bar_performance, bar_missed, bar_novel],
                    ncols=2, width=350, height=260), title=id_))
        else:
            gff_fails = True

    if gff_fails:
        gffcompare_md += ('''
                    __Warning__: Some gffcompare summary cannot be shown.
                    This could be due to  incompatible reference fasta and gff
                    files.''')
    cover_panel = Tabs(tabs=tabs)
    section.markdown(gffcompare_md)
    section.plot(cover_panel)

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

    # Plot overlaps panel:
    section = report.add_section()
    section.markdown('''
    ### Query transfrag classes

    The classes that are assigned by
    [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml),
    which describe the relationship between query transfrag and the most
    similar reference transcript.

    [This diagram](https://ccb.jhu.edu/software/stringtie/
    gffcompare_codes.png) illustrates the different classes.
    ''')

    tracking_dfs = []

    print(gffcompare_outdirs)
    track_files = [x / 'str_merged.tracking' for x in gffcompare_outdirs]

    df_tracking = load_sample_data(
        track_files, sample_ids,
        read_func=lambda x: pd.read_csv(
            x, sep="\t", header=None,
            usecols=[0, 3],
            names=['Count', 'Overlaps']))

    tabs = []
    for id_, df_track in df_tracking.groupby('sample_id'):
        tracking = df_track.groupby("Overlaps").count().reset_index()
        tracking.Overlaps = tracking.Overlaps.map(names)
        tracking['Percent'] = tracking.Count * 100 / tracking.Count.sum()
        tracking = tracking.sort_values("Count", ascending=False)
        track_bar = bars.simple_hbar(
            list(reversed(tracking['Overlaps'].values.tolist())),
            list(reversed(tracking['Percent'].values.tolist())),
            colors=Colors.cerulean, title=id_)

        tracking_dfs.append(tracking)

        tracking['Description'] = pd.Series(tracking.Overlaps.apply(
            lambda x: x.split(':')[0]))

        tracking['Code'] = pd.Series(tracking.Overlaps.apply(
            lambda x: x.split(':')[1]))

        tracking.drop(columns=['sample_id', 'Overlaps'], inplace=True)
        tracking = tracking[['Code', 'Description', 'Count', 'Percent']]
        tracking = tracking.round({
            'Percent': 2
        })

        cols = [TableColumn(
            field=Ci, title=Ci, width=100) for Ci in tracking.columns]

        track_table = DataTable(
            columns=cols, source=ColumnDataSource(tracking),
            index_position=None, width=500)

        tabs.append(Panel(
            child=gridplot([track_bar, track_table], ncols=2), title=id_)
        )

    cover_panel = Tabs(tabs=tabs)
    section.plot(cover_panel)

    def plot_isoforms_per_tpm_bin(
            df_code, class_code, sample_id, geomspace=False):
        """Make plots of number of isoforms per TPM coverage bin."""
        max_ = int(sigfig.round(df_code.TPM.max(), 2))
        if geomspace:
            bins = [math.ceil(x) for x in np.geomspace(10, max_, num=15)]
        else:
            bins = np.linspace(10, max_, 15)
        bins = np.unique(bins)  # Low max_ can end up with duplicated bins

        groups = pd.cut(df_code.TPM, bins).value_counts()
        df_code.to_csv('dfcode.csv')
        df_temp = pd.DataFrame.from_dict(dict(
            x=[x.mid for x in groups.index], y=groups.values
        ))
        df_temp.sort_values(by='x', inplace=True)

        x = [str(math.ceil(x)) for x in df_temp.x]
        y = df_temp.y

        reads_per_iso_geom = \
            bars.simple_bar(x, y,
                            title="{} - Num isoforms/TPM bin - "
                                  "gffcompare class code - '{}'".format(
                                      sample_id, class_code),
                            colors=Colors.cerulean,
                            x_axis_label='TPM',
                            y_axis_label='Number of isoforms')

        reads_per_iso_geom.xaxis.major_label_orientation = math.pi / 2.8
        return reads_per_iso_geom

    log_plots = defaultdict(list)

    try:
        tmap_files = [next(x.glob('*.tmap')) for x in gffcompare_outdirs]
    except StopIteration:
        print("Cannot find .tmap files in {}".format(gffcompare_outdirs))
        return

    df_tmap = load_sample_data(tmap_files, sample_ids)

    for id_, df in df_tmap.groupby('sample_id'):

        log_plots[id_].append(plot_isoforms_per_tpm_bin(
            df, 'all', id_, geomspace=True))

        for class_code, df_code in df.groupby('class_code'):
            log_plots[id_].append(plot_isoforms_per_tpm_bin(
                df_code, class_code, id_, geomspace=True))

    tabs = []
    for id_, sample_plots in log_plots.items():
        tabs.append(Panel(
            child=gridplot(sample_plots, ncols=2), title=id_))

    section.markdown('''
    ### Read coverage by gffcompare transfrag class''')
    cover_panel = Tabs(tabs=tabs)
    section.plot(cover_panel)

    return df_tmap


def pychopper_plots(report, pychop_report):
    """Make plots from pychopper output.

    :param report: aplanat WFReport
    :param pychop_report: path to pychopper stats file
    """
    section = report.add_section()
    section.markdown('''
    ### Pychopper summary statisitcs

    The following plots summarize [pychopper] output
    (https://github.com/nanoporetech/pychopper)

    * **Pr.found**: Reads with primers found in correct orientation at
     both ends.
    * **Resc**: Reads 'rescued' from fused reads
    * **Unusable**: Read with missing or incorrect primer orientation
    * **+/-**: Orientation of reads relative to the mRNA
    ''')

    plots = []
    df = pd.read_csv(pychop_report, sep='\t', index_col=0)

    for id_, df in df.groupby('sample_id'):
        df1 = df.set_index('Name', drop=True)
        df1 = df1.T[['Primers_found', 'Rescue', 'Unusable']]
        df1.rename(columns={
            'Primers_found': 'Pr.found',
            'Rescue': 'Resc',
            'Unusable': 'Un'}, inplace=True)
        df2 = df[df.index == 'Strand']
        bar_chop = bars.simple_bar(
            df1.columns.values.tolist() + df2.Name.values.tolist(),
            df1.iloc[0].values.tolist() + df2.Value.values.tolist(),
            title='{} - Pychopper stats'.format(id_),
            colors=Colors.cerulean)

        plots.extend([bar_chop])
    grid = gridplot(
        plots, ncols=4, width=300, height=300)
    section.plot(grid)


def cluster_quality(cluster_qc_dir, report, sample_ids):
    """Make cluster quality section."""
    section = report.add_section()
    section.markdown('''
    ### De novo clustering quality

    This section shows plots relating to the clustering quality performed
    by isONclust2. The full length reads are mapped to a reference genome
    to create a ground truth of reads mapped to clusters. This is then compared
    to the de novo-generated clusters, and the following statistics are
    generated.

    * [Homogeneity](https://scikit-learn.org/stable/modules/generated/
    sklearn.metrics.homogeneity_score.html): Penalises over-clustering.

    * [Completeness](https://scikit-learn.org/stable/modules/generated/
    sklearn.metrics.completeness_score.html): Penalises under-clustering.

    * [V-measure](https://clusteringjl.readthedocs.io/en/latest/vmeasure.html):
    The harmonic mean of the homogeneity and completeness

    * [Adjusted Rand Index](https://scikit-learn.org/stable/modules/generated/
    sklearn.metrics.adjusted_rand_score.html): Intuitively,  measures the
    percentage of read pairs correctly clustered, normalized so that a perfect
    clustering = 1 and a random cluster assignment achieves = 0

    * NonSingleton: Clusters with multiple reads
    * Singleton: Clusters consisting of a single read (These do not contribute
    to the final transcript calling - I need to check this!)

    ''')

    tabs = []
    for id_, cluster_dir in zip(sample_ids, cluster_qc_dir):
        plots = []
        for fn in ['v_ari_com_hom.csv', 'sing_nonsing.csv']:
            # Skip the next two plots for now
            # 'class_sizes1.csv', 'class_sizes2.csv']:
            df = pd.read_csv(Path(cluster_dir) / fn)
            bar = bars.simple_bar(
                df.Statistic.values.tolist(), df.Value.values.tolist(),
                colors=Colors.cerulean
            )
            bar.xaxis.major_label_orientation = math.pi / 2.8
            plots.append(bar)
        tabs.append(Panel(
            child=gridplot(plots, ncols=4,
                    width=300, height=300), title=id_))

    cover_panel = Tabs(tabs=tabs)
    section.plot(cover_panel)


def transcript_table(report, df_tmaps, covr_threshold):
    """Create searchable table of transcripts."""
    section = report.add_section()

    # Should we put data from each sample into it's own table or have it
    # all in single table and sample_id column? Currently it's the latter

    # drop some columns for the big table and do some filtering
    section.markdown('''
    ### Isoforms table

    Low coverage transcripts are removed to speed up the table viewing. <br>
    Coverage threshold can be set with the parameter
    `transcript_table_cov_thresh`.
    ''')

    df = df_tmaps.drop(
        columns=[
            'FPKM', 'qry_gene_id', 'major_iso_id', 'ref_match_len', 'TPM'])

    if len(df) == 0:
        print("No transcripts found")
        section.markdown("No transcripts found")
        return

    df.sort_values('cov', ascending=True, inplace=True)
    counts = list(range(len(df)))

    # Filter on coverage threshold
    df = df[df['cov'] >= covr_threshold]
    if len(df) < 200:  # Min size of table should be 200. Don't filter
        df = df.sort_values('cov', ascending=False).iloc[:, 0:200]
        covr_threshold = 0

    # Keep Isoforms with coverage > threshold
    vline_x = np.argmax(df['cov'] > covr_threshold)
    vline_y = [0, df['cov'].max()]

    cov_plt = lines.line(
        [counts, [vline_x, vline_x]],  # x-values
        [df['cov'].values.tolist(), vline_y],  # y-values
        title=(
            "Read Coverage. Threshold = {}x coverage".format(
                covr_threshold)
        ), x_axis_label='Num Isoforms',
        y_axis_label='Coverage',
        colors=['blue', 'red'])

    section.plot(cov_plt)

    # Make a column of number of isoforms in parent gene
    gb = df.groupby(['ref_gene_id', 'sample_id']).count()
    # gb = gb.set_index(['ref_gene_id', 'sample_id'])
    gb.rename(columns={'ref_id': 'num_isoforms'}, inplace=True)

    df['parent gene iso num'] = df.apply(
        lambda x: gb.loc[(x.ref_gene_id, x.sample_id), 'num_isoforms'], axis=1)
    # Uncalssified transcritps should not be lumped togetehr
    df.loc[df.class_code == 'u', 'parent gene iso num'] = None

    df.sort_values('parent gene iso num', inplace=True, ascending=True)

    section.table(df, index=False)


def transcriptome_summary(report, gffs, sample_ids, denovo=False):
    """
    Plot transcriptome summaries.

    Some of this data is available via gffcompare output, but the de novo
    pipeline skips that, so we do it al here.

    We do not report exon number for the denovo assembly yet. This is because
    in this case, the gff annotation is generated by aligning to the CDS not
    the genome.

    :param report: aplanat WFReport
    :param gffs: list of paths to gff transcriptome annotations
    :param sample_ids: list of sample ids
    :param denovo: whether annotation was generated by de novo pipeline or not
    """
    # test.db gets written to the git repo.
    section = report.add_section()
    section.markdown('''
    ### Transcriptome summary
    ''')

    tabs = []
    for id_, gff in zip(sample_ids, gffs):

        plots = []

        db = gffutils.create_db(
            gff, dbfn=':memory:', force=True, keep_order=True,
            merge_strategy='merge', sort_attribute_values=True
        )

        num_transcripts = db.count_features_of_type('transcript')
        num_genes = db.count_features_of_type('gene')

        transcript_lens = []
        exons_per_transcript = Counter()
        isoforms_per_gene = []

        for g in db.features_of_type('gene'):

            n_isos = len(list(db.children(g, featuretype='transcript')))
            isoforms_per_gene.append(n_isos)

            for t in db.children(
                    g, featuretype='transcript', order_by='start'):
                tr_len = 0
                exons = list(enumerate(db.children(t, featuretype='exon')))
                if len(exons) == 0:
                    continue
                for nx, ex in exons:
                    tr_len += abs(ex.end - ex.start)

                exons_per_transcript[nx] += 1

                transcript_lens.append(tr_len)

        bar_isos = hist.histogram(
            [isoforms_per_gene], colors=[Colors.cerulean],
            title="isoforms per gene")
        bar_isos.xaxis.axis_label = "Num. isoforms"
        bar_isos.yaxis.axis_label = "Num. genes"

        bar_isos.xaxis.major_label_orientation = math.pi / 2.8
        plots.append(bar_isos)

        box = bars.boxplot_series(
            [id_] * len(transcript_lens), transcript_lens,
            width=70, ylim=(min(transcript_lens), max(transcript_lens)),
            title='transcript lengths')
        plots.append(box)

        if not denovo:
            x, y = zip(*sorted(exons_per_transcript.items()))

            fig = figure(title="Exons per transcript")
            fig.vbar(
                x, top=list(y), color=Colors.cerulean)
            fig.xaxis.axis_label = 'Num. exons'
            fig.yaxis.axis_label = 'Num. genes'

            fig.xaxis.major_label_orientation = math.pi / 2.8
            plots.append(fig)

        df_sum = pd.DataFrame.from_dict(
            {'Total genes': [num_genes],
             'Total transcripts': [num_transcripts],
             'Max trans. len': max(transcript_lens),
             'Min trans. len': min(transcript_lens)}).T
        df_sum.reset_index(drop=False, inplace=True)

        df_sum.columns = [' ', 'count']
        cols = [TableColumn(
            field=Ci, title=Ci, width=80) for Ci in df_sum.columns]
        data_table = DataTable(
            columns=cols, source=ColumnDataSource(df_sum),
            index_position=None, width=180)
        plots.append(data_table)

        tabs.append(Panel(
            child=gridplot(plots, ncols=4,
                           width=300, height=300), title=id_))

    cover_panel = Tabs(tabs=tabs)
    section.plot(cover_panel)


def load_sample_data(files, sample_ids, read_func=None):
    """Load CSVs into dataframe, and assign sample_id column."""
    df_ = pd.DataFrame()
    if not files:
        return None
    for id_, x in zip(sample_ids, files):
        if read_func:
            d = read_func(x)
        else:
            d = pd.read_csv(x, sep='\t+')
        d['sample_id'] = id_
        df_ = pd.concat([df_, d])
    return df_


def seq_stats_tabs(report, sample_ids, stats):
    """Make tabs of sequence summaries by sample."""
    tabs = []
    for id_, summ in sorted(zip(sample_ids, stats)):
        df_sum = pd.read_csv(summ, index_col=False, sep='\t')
        rlp = read_length_plot(df_sum)
        rqp = read_quality_plot(df_sum)
        grid = gridplot(
            [rlp, rqp], ncols=2, sizing_mode="stretch_width")

        tabs.append(Panel(child=grid, title=id_))
    section = report.add_section()
    section.markdown("""
    ### Sequence summaries""")
    section.plot(Tabs(tabs=tabs))


def jaffal_table(report, result_csv):
    """Make a table of fusion transcripts identified by JAFFAL."""
    cols = [
        'sample_id', 'fusion genes', 'chrom1', 'chrom2', 'spanning reads',
        'classification', 'known']

    df = pd.read_csv(result_csv)
    sid_col = df.pop('sample_id')
    df.insert(0, 'sample_id', sid_col)

    df = df[cols]
    df['chroms'] = df.chrom1.astype(str) + ':' + df.chrom2.astype(str)
    df.rename(columns={
        'spanning reads': 'nreads',
        'fusion genes': 'genes'}, inplace=True)
    df.drop(columns=['chrom1', 'chrom2'], inplace=True)
    section = report.add_section()
    section.markdown("""
    ### JAFFAL fusion transcript summary

    This table summarizes putative fusion transcripts identified
    by [JAFFAL](https://github.com/Oshlack/JAFFA/).

    * genes: the gene symbols of the fusion partners
    * nreads: The number of reads supporting the fusion
    * classification: JAFFAL's classification
    * known: whether this fusion is in the given set of known gene fusions
    * chroms: the respective, original chromosome location of the two partner
    genes
    """)
    section.table(df)


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", help="Report output file")
    parser.add_argument("--summaries", nargs='+', help="Read summary file.")
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
        "--alignment_stats", required=False, default=None, nargs='*',
        help="TSV summary file of alignment statistics")
    parser.add_argument(
        "--gff_annotation", required=True, nargs='+',
        help="transcriptome annotation gff file")
    parser.add_argument(
        "--gffcompare_dir", required=False, default=None, nargs='*',
        help="gffcompare outout dir")
    parser.add_argument(
        "--pychop_report", required=False, default=None,
        help="TSV summary file of pychopper statistics")
    parser.add_argument(
        "--sample_ids", required=True, nargs='+',
        help="List of sample ids")
    parser.add_argument(
        "--transcript_table_cov_thresh", required=False, type=int, default=50,
        help="Isoforms without this support will be excluded from the table")
    parser.add_argument(
        "--cluster_qc_dirs", required=False, type=str, default=None, nargs='*',
        help="Directory with various cluster quality csvs")
    parser.add_argument(
        "--jaffal_csv", required=False, type=str, default=None,
        help="Path to JAFFAL results csv")
    parser.add_argument('--denovo', dest='denovo', action='store_true')

    args = parser.parse_args()

    sample_ids = args.sample_ids

    report = WFReport(
        "Transcript isoform report", "wf-transcriptomes",
        revision=args.revision, commit=args.commit)

    # QC
    seq_stats_tabs(report, args.sample_ids, args.summaries)

    if args.alignment_stats is not None:
        df_aln_stats = load_sample_data(args.alignment_stats, sample_ids)
        section = report.add_section()
        section.markdown('''
          ### Read mapping summary

          Summary of minimap2 mapping from
          [seqkit](https://bioinf.shenwei.me/seqkit/)
          `seqkit bam -s`''')

        section.table(df_aln_stats)

    if args.pychop_report is not None:
        pychopper_plots(report, args.pychop_report)

    # Results
    transcriptome_summary(
        report, args.gff_annotation, sample_ids, denovo=args.denovo)

    df_tmaps = gff_compare_plots(
        report,
        [Path(x) for x in args.gffcompare_dir],
        sample_ids)

    report.write(args.report)

    if df_tmaps is not None:
        transcript_table(report, df_tmaps, args.transcript_table_cov_thresh)

    if args.cluster_qc_dirs is not None:
        cluster_quality(args.cluster_qc_dirs, report, sample_ids)

    if args.jaffal_csv:
        jaffal_table(report, args.jaffal_csv)

    # Arguments and software versions
    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    report.write(args.report)


if __name__ == "__main__":
    main()
