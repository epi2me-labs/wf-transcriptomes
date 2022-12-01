#!/usr/bin/env python
"""Dynamically generate isONclust2 processes."""
from collections import OrderedDict
from glob import glob
from itertools import zip_longest
from pathlib import Path
import re
import subprocess as sub


class Node:
    """Node."""

    def __init__(self, node_id, file_, left, right, parent, level):
        """Set node attaributes."""
        self.Id = node_id
        self.File = file_
        self.Left = left
        self.Right = right
        self.Parent = parent
        self.Level = level
        self.Done = False
        self.RightSide = False

    def __repr__(self):
        """Get string repr of a node."""
        return "Node:{} Level: {} File: {} Done: " \
            "{} Left: {} Right: {} Parent: {}".format(
                self.Id, self.Level,
                self.File, self.Done, self.Left.Id if
                self.Left is not None else None,
                self.Right.Id if self.Right is not None else None,
                self.Parent.Id if self.Parent is not None else None)


def grouper(n, iterable, fillvalue=None):
    """
    Group adjacent nodes.

    grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx.
    """
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def build_job_tree():
    """Build a job tree of nodes."""
    job_tree = OrderedDict()
    batches = glob("batches/isONbatch_*.cer")
    batch_ids = [
        int(re.search(
            'batches/isONbatch_(.*)\\.cer$', x).group(1))
        for x in batches]
    levels = OrderedDict()
    levels[0] = []
    for id_, bf in sorted(zip(batch_ids, batches), key=lambda x: x[0]):
        n = Node(
            id_,
            "clusters/isONcluster_{}.cer".format(id_),
            None,
            None,
            None,
            0)
        n.Done = True
        job_tree[id_] = n
        levels[0].append(n)
    level = 0
    max_id = levels[0][-1].Id
    while len(levels[level]) != 1:  # Final level will be link
        next_level = level + 1
        levels[next_level] = []
        for l_, r in grouper(2, levels[level]):
            if r is None:  # End of a level
                levels[level].pop()  # remove last node?
                l_.Level += 1  # ncrement level
                levels[next_level].append(l_)  # Add the left to the next level
                continue
            max_id += 1
            new_batch = "clusters/isONcluster_{}.cer".format(max_id)
            new_node = Node(max_id, new_batch, l_, r, None, next_level)
            l_.Parent = new_node
            r.Parent = new_node
            r.RightSide = True
            levels[next_level].append(new_node)
            job_tree[max_id] = new_node
        level = next_level
    root = job_tree[len(job_tree) - 1].Id
    job_tree[root].RightSide = True

    return job_tree, levels


def main():
    """Entry point."""
    Path('clusters').mkdir(exist_ok=True)
    job_tree, levels = build_job_tree()

    init_template = (
        'isONclust2 cluster -x {} -v -Q -l batches/isONbatch_{}.cer '
        '-o clusters/isONcluster_{}.cer {}; '
        'sync;\n')
    template = (
        'isONclust2 cluster -x {} -v -Q -l clusters/isONcluster_{}.cer '
        '-r clusters/isONcluster_{}.cer -o clusters/isONcluster_{}.cer '
        '{}; sync\n')

    for nr, l in levels.items():
        jobs_out = 'jobs_level_{}.sh'.format(nr)
        with open(jobs_out, 'w') as fh:
            for n in l:
                purge = "-z" if n.RightSide else ""
                if nr == 0 or n.Left is None or n.Right is None:
                    jr = init_template.format('sahlin', n.Id, n.Id, purge)
                    fh.write(jr)
                else:
                    jr = template.format(
                        'sahlin', n.Left.Id, n.Right.Id, n.Id, purge)
                    fh.write(jr)
        # Run a level in parallel
        cmd = "parallel < {}".format(jobs_out)
        sub.call(cmd, shell=True)
    sub.call((
        "ln -s `realpath clusters/isONcluster_{}.cer` "
        "isONcluster_ROOT.cer".format(n.Id)), shell=True)


if __name__ == '__main__':
    # The cwd should be the process dir that contains 'batches/'
    main()
