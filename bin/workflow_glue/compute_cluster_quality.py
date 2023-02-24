#!/usr/bin/env python
"""Generate cluster quality data."""

# Adapted form script by Kristoffer Sahlin for
# isONclust: https://github.com/ksahlin/isONclust

from collections import defaultdict
import math
from pathlib import Path
import sys

import matplotlib
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import pysam
from sklearn.metrics.cluster import (
    adjusted_rand_score, completeness_score,
    homogeneity_score, v_measure_score)

from .util import wf_parser  # noqa: ABS101

matplotlib.use('Agg')


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("compute_cluster_quality")
    parser.add_argument(
        '--clusters',
        type=str,
        help='Inferred clusters (tsv file)')
    parser.add_argument(
        '--classes',
        type=str,
        help='A sorted and indexed bam file.')
    parser.add_argument(
        '--ctsv',
        default=None,
        type=str,
        help='Write true classes in this TSV file.')
    parser.add_argument(
        '--simulated',
        action="store_true",
        help='Simulated data, we can simply read correct classes '
             'from the ref field.')
    parser.add_argument(
        '--ont',
        action="store_true",
        help='ONT data, parsing accessions differently.')
    parser.add_argument(
        '--modified_ont',
        action="store_true",
        help='ONT data preprocessed accessions, parsing '
             'accessions differently.')
    parser.add_argument('--outfile', type=str, help='Output file with results')
    parser.add_argument(
        '--report',
        type=str,
        help='Output PDF file with report')
    parser.add_argument('--sizes', type=str, help='Cluster sizes')
    parser.add_argument(
        '--raw_data_out',
        type=str,
        help='dir to save raw data for plotting')

    return parser


def parse_inferred_clusters_tsv(tsv_file, args):
    """parse_inferred_clusters_tsv."""
    infile = open(tsv_file, "r")
    infile.readline()
    clusters = {}
    for line in infile:
        cluster_id, _, read_acc = line.strip().split("\t")
        if args.simulated:
            read_acc = "_".join([item for item in read_acc.split("_")[:-1]])
        elif args.ont:
            read_acc = read_acc.split(" ")[0]
        elif args.modified_ont:
            read_acc = read_acc
        else:
            read_acc = read_acc.split("_strand")[0]

        clusters[read_acc] = int(cluster_id)
    return clusters


def parse_true_clusters(ref_file):
    """parse_true_clusters."""
    classes = defaultdict(dict)
    ref_id_to_chrom = {}
    alignment_counter = defaultdict(int)

    prev_chrom = -1
    curr_class_id = -1
    prev_class_stop = -1
    prev_read_id = ""
    unique_reads = set()
    unclassified = 0
    for read in ref_file.fetch(until_eof=True):
        unique_reads.add(read.query_name)
        if read.is_unmapped:
            unclassified += 1
            continue
        # deal with supplementary alignments!!
        if read.is_secondary or read.is_supplementary:
            continue
        assert prev_read_id != read.query_name

        chrom = read.reference_name
        if chrom != prev_chrom:
            curr_class_id += 1
            classes[read.query_name] = curr_class_id
            prev_chrom = chrom
            prev_class_stop = read.reference_end

        else:
            read_ref_start = read.reference_start
            if read_ref_start > prev_class_stop:
                curr_class_id += 1
                classes[read.query_name] = curr_class_id
                prev_class_stop = read.reference_end
            else:
                classes[read.query_name] = curr_class_id

        prev_class_stop = max(read.reference_end, prev_class_stop)
        prev_read_id = read.query_name

        # classes[read.query_name] = int(read.reference_id)  # chrom
        ref_id_to_chrom[int(read.reference_id)] = chrom
        alignment_counter[int(read.reference_id)] += 1
        # if chrom not in class_ranges:
        #     class_ranges[chrom] = {}

        # for start, stop in class_ranges[chrom]:
        #     if start <= read_ref_start and  read_ref_end <= stop:
        #         # entirly within
        #     elif

        # classes[read.query_name] =  #read.reference_name.split("|")[0]
        # #"|".join(read.reference_name.split("|")[:2])
    return classes, len(unique_reads), unclassified


def parse_true_clusters_simulated(ref_file):
    """parse_true_clusters_simulated."""
    classes = defaultdict(dict)
    for read in ref_file.fetch(until_eof=True):
        # by gene id
        classes[read.query_name] = read.reference_name.split("|")[0]
        # by transcript id
        # classes[read.query_name] = read.reference_name.split("|")[1]

    return classes


def compute_v_measure(clusters, classes):
    """compute_v_measure."""
    class_list, cluster_list = [], []
    # not_found_id = 1000000
    clustered_but_unaligned = 0
    for read in clusters:
        if read in classes:
            class_list.append(classes[read])
            cluster_list.append(clusters[read])
        else:
            clustered_but_unaligned += 1

    # added the unprocessed reads to the measure
    not_clustered = set(classes.keys()) - set(clusters.keys())
    highest_cluster_id = max(clusters.values())
    highest_cluster_id += 1
    for read in not_clustered:
        class_list.append(classes[read])
        cluster_list.append(highest_cluster_id)
        highest_cluster_id += 1

    v_score = v_measure_score(class_list, cluster_list)
    compl_score = completeness_score(class_list, cluster_list)
    homog_score = homogeneity_score(class_list, cluster_list)
    ari = adjusted_rand_score(class_list, cluster_list)

    sys.stdout("Not included in clustering but aligned:", len(not_clustered))
    sys.stdout(
        "v:",
        v_score,
        "Completeness:",
        compl_score,
        "Homogeneity:",
        homog_score)
    sys.stdout(
        "Nr reads clustered but unaligned "
        "(i.e., no class and excluded from v-measure): ",
        clustered_but_unaligned)
    return v_score, compl_score, homog_score, clustered_but_unaligned, ari


def compute_v_measure_non_singleton_classes(clusters, classes):
    """V measure for non-singleton classes."""
    max_cluster_id = max(clusters.values())
    new_id = max_cluster_id + 1
    classes_dict = {}
    for read_acc, cl_id in classes.items():
        if cl_id not in classes_dict:
            classes_dict[cl_id] = [read_acc]
        else:
            classes_dict[cl_id].append(read_acc)

    nontrivial_classes_reads = []
    for cl_id in classes_dict:
        if len(classes_dict[cl_id]) < 5:
            continue
        else:
            for read in classes_dict[cl_id]:
                nontrivial_classes_reads.append(read)

    class_list, cluster_list = [], []
    for read in nontrivial_classes_reads:
        if read in clusters:
            class_list.append(classes[read])
            cluster_list.append(clusters[read])
        else:
            class_list.append(classes[read])
            cluster_list.append(new_id)
            new_id += 1

    v_score = v_measure_score(class_list, cluster_list)
    compl_score = completeness_score(class_list, cluster_list)
    homog_score = homogeneity_score(class_list, cluster_list)
    nr_filtered_classes = len(
        [1 for cl_id in classes_dict if len(classes_dict[cl_id]) >= 5])
    sys.stdout(
        "NONTRIvIAL CLASSES: v:",
        v_score,
        "Completeness:",
        compl_score,
        "Homogeneity:",
        homog_score)
    sys.stdout("NUMBER OF CLASSES (FILTERED):", len(
        [1 for cl_id in classes_dict if len(classes_dict[cl_id]) >= 5]))
    return v_score, compl_score, homog_score, nr_filtered_classes


def compute_v_measure_non_singletons(clusters, classes):
    """V measure for non-singletons."""
    cluster_dict = {}
    for read_acc, cl_id in clusters.items():
        if cl_id not in cluster_dict:
            cluster_dict[cl_id] = [read_acc]
        else:
            cluster_dict[cl_id].append(read_acc)

    nontrivial_clustered_reads = []
    for cl_id in cluster_dict:
        if len(cluster_dict[cl_id]) <= 1:
            continue
        else:
            for read in cluster_dict[cl_id]:
                nontrivial_clustered_reads.append(read)

    class_list, cluster_list = [], []
    # not_found_id = 1000000
    clustered_but_unaligned = 0
    for read in nontrivial_clustered_reads:
        if read in classes:
            class_list.append(classes[read])
            cluster_list.append(clusters[read])
        else:
            # sys.stdout("Read was clustered but unaligned:", read)
            clustered_but_unaligned += 1

    v_score = v_measure_score(class_list, cluster_list)
    compl_score = completeness_score(class_list, cluster_list)
    homog_score = homogeneity_score(class_list, cluster_list)
    sys.stdout(
        "NONTRIvIAL CLUSTERS: v:",
        v_score,
        "Completeness:",
        compl_score,
        "Homogeneity:",
        homog_score)
    sys.stdout(
        "NONTRIvIAL CLUSTERS: Nr reads clustered but unaligned "
        "(i.e., no class and excluded from v-veasure): ",
        clustered_but_unaligned)
    return v_score, compl_score, homog_score, clustered_but_unaligned


def percentile(n, percent, key=lambda x: x):
    """
    Find the percentile of a list of values.

    @parameter n - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value
    from each element of N.

    @return - the percentile of the values
    """
    if not n:
        return None
    k = (len(n) - 1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(n[int(k)])
    d0 = key(n[int(f)]) * (c - k)
    d1 = key(n[int(c)]) * (k - f)
    return d0 + d1

# end of http://code.activestate.com/recipes/511478/ }}}


def get_cluster_information(clusters, classes):
    """Get cluster info."""
    # class distribution
    class_dict = {}
    for read_acc, class_id in classes.items():
        if class_id not in class_dict:
            class_dict[class_id] = [read_acc]
        else:
            class_dict[class_id].append(read_acc)

    total_nr_classes = len(class_dict)
    class_distribution = sorted([len(cl) for cl in class_dict.values()])
    singleton_classes = set(
        [acc_list[0] for cl_id, acc_list in class_dict.items()
         if len(acc_list) == 1])
    min_class_size = min(class_distribution)
    max_class_size = max(class_distribution)
    mean_class_size = sum(class_distribution) / float(len(class_distribution))
    median_class_size = class_distribution[int(len(class_distribution) / 2)] \
        if len(class_distribution) % 2 == 1 else (
        class_distribution[int(len(class_distribution) / 2)] +
        class_distribution[int(len(class_distribution) / 2) - 1]) / 2.0

    upper_75_class_size = percentile(class_distribution, 0.75)
    median_class_size = percentile(class_distribution, 0.5)

    tot_size = sum(class_distribution)
    e_class_size = sum(
        [c_s**2 for c_s in class_distribution]) / float(tot_size)
    tot_iterated_size = 0
    for c_s in class_distribution[::-1]:
        tot_iterated_size += c_s
        if tot_iterated_size >= tot_size / 2.0:
            n50_class_size = c_s
            break

    # cluster distribution
    cluster_dict = {}
    for read_acc, cl_id in clusters.items():
        if cl_id not in cluster_dict:
            cluster_dict[cl_id] = [read_acc]
        else:
            cluster_dict[cl_id].append(read_acc)

    cluster_distribution = sorted([len(cl) for cl in cluster_dict.values()])

    # in case unclustered reads are missing from output (as for isoseq3)
    omitted_from_output_singletons = set(classes.keys()) - set(clusters.keys())
    cluster_distribution = [
        1 for i in range(
            len(omitted_from_output_singletons))] + cluster_distribution
    total_nr_clusters = len(cluster_distribution)

    singleton_clusters = set(
        [acc_list[0] for cl_id, acc_list in cluster_dict.items()
         if len(acc_list) == 1])
    min_cluster_size = min(cluster_distribution)
    max_cluster_size = max(cluster_distribution)
    mean_cluster_size = sum(cluster_distribution) / \
        float(len(cluster_distribution))

    upper_75_cluster_size = percentile(cluster_distribution, 0.75)
    median_cluster_size = percentile(cluster_distribution, 0.5)

    tot_size = sum(cluster_distribution)
    e_cluster_size = sum(
        [c_s**2 for c_s in cluster_distribution]) / float(tot_size)
    tot_iterated_size = 0
    for c_s in cluster_distribution[::-1]:
        tot_iterated_size += c_s
        if tot_iterated_size >= tot_size / 2.0:
            n50_cluster_size = c_s
            break

    unaligned_but_nontrivially_clustered = set(
        clusters.keys()) - singleton_clusters - set(classes.keys())

    # not_considered = set([read for read in classes if read not in clusters ])

    not_clustered_classes = defaultdict(int)
    clustered_classes = defaultdict(int)
    reads_not_clustered = defaultdict(list)
    for read in clusters:
        if read in classes:
            class_id = classes[read]
        else:
            class_id = "unaligned"

        if read in singleton_clusters:
            not_clustered_classes[class_id] += 1
            reads_not_clustered[class_id].append(read)
        else:
            clustered_classes[class_id] += 1

    sys.stdout("UNCLUSTERED:", "Tot classes:", len(not_clustered_classes))
    sys.stdout("CLUSTERED:", "Tot classes:", len(clustered_classes))
    sys.stdout("MIXED:", "Tot classes containing both:", len(
        set(clustered_classes.keys()) & set(not_clustered_classes.keys())))
    sys.stdout("Total number of classes (unique gene ID):", total_nr_classes)
    return (
        total_nr_classes - len(singleton_classes),
        len(singleton_classes),
        min_class_size,
        max_class_size,
        mean_class_size,
        median_class_size,
        total_nr_clusters,
        len(singleton_clusters) + len(omitted_from_output_singletons),
        min_cluster_size,
        max_cluster_size,
        mean_cluster_size,
        median_cluster_size,
        len(unaligned_but_nontrivially_clustered),
        upper_75_class_size,
        upper_75_cluster_size,
        e_class_size,
        n50_class_size,
        e_cluster_size,
        n50_cluster_size)


def main(args):
    """Entry point."""
    clusters = parse_inferred_clusters_tsv(args.clusters, args)
    if not clusters:
        outfile = open(args.outfile, "w")
        outfile.write("No clusters created\n")
        outfile.close()
        return

    if args.simulated:
        ref_file = pysam.AlignmentFile(args.classes, "r", check_sq=False)
        classes = parse_true_clusters_simulated(ref_file)
        # by simulation we know classes of all reads, they are therefore the
        # same number.
        tot_nr_reads = len(classes)
    else:
        ref_file = pysam.AlignmentFile(args.classes, "rb", check_sq=False)
        classes, tot_nr_reads, unclassified = parse_true_clusters(ref_file)

    v_score, compl_score, homog_score, clustered_but_unaligned, ari = \
        compute_v_measure(clusters, classes)

    (
        nr_non_singleton_classes, singleton_classes, min_class_size,
        max_class_size, mean_class_size, median_class_size, total_nr_clusters,
        singleton_clusters, min_cluster_size, max_cluster_size,
        mean_cluster_size, median_cluster_size,
        unaligned_but_nontrivially_clustered,
        upper_75_class_size, upper_75_cluster_size, e_class_size,
        n50_class_size, e_cluster_size,
        n50_cluster_size
    ) = get_cluster_information(clusters, classes)

    outfile = open(args.outfile, "w")

    outfile.write("CLASSES\n")

    # reads, unaligned, classes, singleton, min,max, mean,median

    outfile.write(
        "{0},{1},{2},{3},{4},{5},{6},{7}\n".format(
            "tot_nr_reads",
            "unclassified",
            "nr_non_singleton_classes",
            "singleton_classes",
            "upper_75_class_size",
            "median_class_size",
            "e_class_size",
            "n50_class_size"))
    outfile.write(
        "{0},{1},{2},{3},{4},{5},{6},{7}\n".format(
            tot_nr_reads,
            unclassified,
            nr_non_singleton_classes,
            singleton_classes,
            upper_75_class_size,
            median_class_size,
            e_class_size,
            n50_class_size))

    # reads_nontrivially_clustered_(%), Singletons_(%),
    # reads_nontrivially_clustered_but_unaligned, v, c,h ,v_nt, c_nt,h_nt,
    # non_singleton_clusters, min, max, median, mean

    reads_nontrivially_clustered_percent = round(
        100 * (float(tot_nr_reads - singleton_clusters) / tot_nr_reads), 1)
    # round(1.0 - reads_nontrivially_clustered_percent, 2)
    reads_nontrivially_clustered_but_unaligned = \
        unaligned_but_nontrivially_clustered
    v, c, h = round(v_score, 3), round(compl_score, 3), round(homog_score, 3)

    non_singleton_clusters = total_nr_clusters - singleton_clusters

    sys.stdout(
        "NONTRIVIAL CLUSTERS: ", (total_nr_clusters - singleton_clusters))

    outfile.write("CLUSTERS\n")
    outfile.write(
        "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n".format(
            "v",
            "c",
            "h",
            "ARI",
            "reads_nontrivially_clustered_percent",
            "reads_nontrivially_clustered_but_unaligned",
            "non_singleton_clusters",
            "singleton_clusters",
            "upper_75_cluster_size",
            "median",
            "e_cluster_size",
            "n50_cluster_size"))
    outfile.write(
        "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11}\n".format(
            v,
            c,
            h,
            ari,
            reads_nontrivially_clustered_percent,
            reads_nontrivially_clustered_but_unaligned,
            non_singleton_clusters,
            singleton_clusters,
            upper_75_cluster_size,
            median_cluster_size,
            e_cluster_size,
            n50_cluster_size))
    outfile.close()

    if args.ctsv is not None:
        cfh = open(args.ctsv, "w")
        cfh.write("Read\tCluster\n")
        for r, c in classes.items():
            cfh.write("{}\t{}\n".format(r, c))
        cfh.flush()
        cfh.close()

    dfc = pd.DataFrame(
        {
            'Statistic': [
                'v-measure',
                'ARI',
                'Completeness',
                'Homogeneity'],
            'value': [
                v,
                ari,
                c,
                h]}).set_index('Statistic')
    dfn = pd.DataFrame(
        {'Statistic': ['NonSingleton',
                       'Singletons'],
         'value': [non_singleton_clusters,
                   singleton_clusters]}).set_index('Statistic')
    dfs = pd.DataFrame(
        {
            'Statistic': [
                'Upper75ClsSize',
                'Upper75ClassSize',
                'MedianClsSize',
                'MedianClassSize'],
            'value': [
                upper_75_cluster_size,
                upper_75_class_size,
                median_cluster_size,
                median_class_size]}).set_index('Statistic')
    dfs2 = pd.DataFrame(
        {'Statistic': ['N50ClsSize', 'N50ClassSize'],
            'value': [n50_cluster_size, n50_class_size]
         }).set_index('Statistic')

    rdo = Path(args.raw_data_out)
    dfc.to_csv(rdo / 'v_ari_com_hom.csv')
    dfn.to_csv(rdo / 'sing_nonsing.csv')
    dfs.to_csv(rdo / 'class_sizes1.csv')
    dfs2.to_csv(rdo / 'class_sizes2.csv')

    pages = PdfPages(args.report)

    yd = 7 * 2.5

    ax = dfc.plot(kind='bar', fontsize=7, rot=0)
    for p in ax.patches:
        ax.annotate(
            "{:.2f}".format(
                p.get_height()),
            (p.get_x() +
             p.get_width() /
             2.0,
             p.get_height()),
            ha='center')
    pages.savefig()
    plt.clf()

    ax = dfn.plot(kind='bar', fontsize=7, rot=0)
    for p in ax.patches:
        ax.annotate(
            "{:.2f}".format(
                p.get_height()),
            (p.get_x() +
             p.get_width() /
             2.0,
             p.get_height()),
            ha='center')
    pages.savefig()
    plt.clf()

    ax = dfs.plot(kind='bar', fontsize=7, rot=0)
    for p in ax.patches:
        ax.annotate(
            "{:.2f}".format(
                p.get_height()),
            (p.get_x() +
             p.get_width() /
             2.0,
             p.get_height() +
             yd),
            ha='center')
    pages.savefig()
    plt.clf()
    ax = dfs2.plot(kind='bar', fontsize=7, rot=0)
    for p in ax.patches:
        ax.annotate(
            "{:.2f}".format(
                p.get_height()),
            (p.get_x() +
             p.get_width() /
             2.0,
             p.get_height() +
             yd),
            ha='center')
    pages.savefig()
    plt.clf()

    pages.close()
