import sys
import time
import os
from os import path
import pandas as pd
from glob import glob
import re
from threading import Lock
from itertools import zip_longest
from collections import OrderedDict, namedtuple

if not workflow.overwrite_configfiles:
    configfile: "config.yml"

WORKDIR = path.abspath(path.join(config["workdir_top"], config["pipeline"]))
workdir: WORKDIR
SNAKEDIR = path.dirname(workflow.snakefile)

include: "snakelib/utils.snake"

in_fastq = config["reads_fastq"]
if not path.isabs(in_fastq):
	in_fastq = path.join(SNAKEDIR, in_fastq)


class Node:
    def __init__(self, Id, File, Left, Right, Parent, Level):
        self.Id = Id
        self.File = File
        self.Left = Left
        self.Right = Right
        self.Parent = Parent
        self.Level = Level
        self.Done = False
        self.RightSide = False
    def __repr__(self):
        return "Node:{} Level: {} File: {} Done: {} Left: {} Right: {} Parent: {}".format(self.Id, self.Level, self.File, self.Done, self.Left.Id if self.Left is not None else None, self.Right.Id if self.Right is not None else None, self.Parent.Id if self.Parent is not None else None)

def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)

def build_job_tree(batch_dir):
    batches = glob("{}/isONbatch_*.cer".format(batch_dir))
    batch_ids = [int(re.search('/isONbatch_(.*)\.cer$', x).group(1)) for x in batches]
    global LEVELS
    global JOB_TREE
    LEVELS = OrderedDict()
    LEVELS[0] = []
    for Id, bf in sorted(zip(batch_ids, batches), key=lambda x: x[0]):
        n = Node(Id, "clusters/isONcluster_{}.cer".format(Id), None, None, None, 0)
        n.Done = True
        JOB_TREE[Id] = n
        LEVELS[0].append(n)
    level = 0
    max_id = LEVELS[0][-1].Id
    while len(LEVELS[level]) != 1:
        next_level = level + 1
        LEVELS[next_level] = []
        for l, r in grouper(2, LEVELS[level]):
            if r is None:
                LEVELS[level].pop()
                l.Level += 1
                LEVELS[next_level].append(l)
                continue
            max_id += 1
            new_batch = "clusters/isONcluster_{}.cer".format(max_id)
            new_node = Node(max_id, new_batch, l, r, None, next_level)
            l.Parent = new_node
            r.Parent =  new_node
            r.RightSide = True
            LEVELS[next_level].append(new_node)
            JOB_TREE[max_id] = new_node
        level = next_level
    global ROOT
    ROOT = JOB_TREE[len(JOB_TREE)-1].Id
    JOB_TREE[ROOT].RightSide = True
    print("Merge clustering job tree nodes:",file=sys.stderr)
    for n in JOB_TREE.values():
        print("\t{}".format(n),file=sys.stderr)


def generate_rules(levels, snk):
    init_template = """
rule cluster_job_%d:
    input:
        left = "sorted/batches/isONbatch_%d.cer",
    output: "clusters/isONcluster_%d.cer"
    shell: "isONclust2 cluster -x %s -v -Q -l %s -o %s %s; sync"

    """
    template = """
rule cluster_job_%d:
    input:
        left = "clusters/isONcluster_%d.cer",
        right = "clusters/isONcluster_%d.cer",
    output: "clusters/isONcluster_%d.cer"
    shell: "isONclust2 cluster -x %s -v -Q -l %s -r %s -o %s %s; sync"

    """

    link_template="""
rule link_root:
    input: "clusters/isONcluster_%d.cer",
    output: "clusters/isONcluster_ROOT.cer"
    shell: "ln -s `realpath {input}` {output}"

    """

    fh = open(snk, "w")
    for nr, l in levels.items():
        for n in l:
            purge = "-z" if n.RightSide else ""
            if nr == 0 or n.Left is None or n.Right is None:
                jr = init_template % (n.Id, n.Id, n.Id, config["cls_mode"], "{input.left}", "{output}", purge)
                fh.write(jr)
            else:
                jr = template % (n.Id, n.Left.Id, n.Right.Id, n.Id, config["cls_mode"], "{input.left}", "{input.right}", "{output}", purge)
                fh.write(jr)
    global ROOT
    fh.write(link_template % ROOT)
    fh.write("\nROOT = %d" % ROOT)
    fh.flush()
    fh.close()

def count_fastq_bases(fname, size=128000000):
    fh = open(fname, "r")
    count = 0
    while True:
        b = fh.read(size)
        if not b:
            break
        count += b.count("A")
        count += b.count("T")
        count += b.count("G")
        count += b.count("C")
        count += b.count("U")
    fh.close()
    return count

def preprocess_reads(fq):
    pc_opts = config["pychopper_opts"]
    concat = config["concatenate"]
    thr = config["cores"]

    out_fq = "processed_reads/full_length_reads.fq"
    if os.path.isdir("processed_reads"):
        return out_fq

    shell("mkdir -p processed_reads")
    if concat:
        print("Concatenating reads under directory: " + fq)
        shell("find %s  -regextype posix-extended -regex '.*\.(fastq|fq)$' -exec cat {{}} \\; > processed_reads/input_reads.fq" % fq)
    else:
        shell("ln -s `realpath %s` processed_reads/input_reads.fq" % fq)

    if config["run_pychopper"]:
        print("Running pychopper of fastq file: processed_reads/input_reads.fq")
        shell("(cd processed_reads; cdna_classifier.py -t %d %s input_reads.fq full_length_reads.fq)" % (thr, pc_opts))
    else:
        shell("ln -s `realpath processed_reads/input_reads.fq` processed_reads/full_length_reads.fq")

    return out_fq


ROOT = None
DYNAMIC_RULES="job_rules.snk"
if ((not os.path.isfile(os.path.join(WORKDIR,"sorted","sorted_reads.fastq"))) or (not os.path.isfile(os.path.join(SNAKEDIR, DYNAMIC_RULES)))):
    print("Preprocessing read in fastq file:", in_fastq)
    proc_fastq = preprocess_reads(in_fastq)
    print("Counting records in input fastq:", proc_fastq)
    nr_bases = count_fastq_bases(proc_fastq)
    print("Bases in input: {} megabases".format(int(nr_bases/10**6)))
    if config['batch_size'] < 0:
        config['batch_size'] = int(nr_bases/1000/config["cores"])

    print("Batch size is: {}".format(config['batch_size']))

    init_cls_options = """ --batch-size {} --kmer-size {} --window-size {} --min-shared {} --min-qual {}\
                         --mapped-threshold {} --aligned-threshold {} --min-fraction {} --min-prob-no-hits {} -M {} -P {} -g {} -c {} -F {} """
    init_cls_options = init_cls_options.format(config["batch_size"], config["kmer_size"], config["window_size"], config["min_shared"], config["min_qual"], \
                        config["mapped_threshold"], config["aligned_threshold"], config["min_fraction"], config["min_prob_no_hits"], config["batch_max_seq"], config["consensus_period"],
    config["consensus_minimum"], config["consensus_maximum"], config["min_left_cls"])

    shell("""
            rm -fr clusters sorted
            mkdir -p sorted; isONclust2 sort {} -v -o sorted {};
            mkdir -p clusters;
        """.format(init_cls_options, proc_fastq))
    JOB_TREE = OrderedDict()
    LEVELS = None

    build_job_tree("sorted/batches")
    generate_rules(LEVELS, "{}/job_rules.snk".format(SNAKEDIR))

include: DYNAMIC_RULES



rule all:
    input: rules.link_root.output
    output: directory("final_clusters")
    shell:
        """ isONclust2 dump -v -i sorted/sorted_reads_idx.cer -o final_clusters {input}; sync """