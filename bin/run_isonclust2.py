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

    def __init__(self, Id, File, Left, Right, Parent, Level):
        """Set node attaributes."""
        self.Id = Id
        self.File = File
        self.Left = Left
        self.Right = Right
        self.Parent = Parent
        self.Level = Level
        self.Done = False
        self.RightSide = False

    def __repr__(self):
        """Get string repr of a node."""
        return "Node:{} Level: {} File: {} Done: {} Left: {} Right: " \
               "{} Parent: {}".format(
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
    JOB_TREE = OrderedDict()
    batches = glob("batches/isONbatch_*.cer")
    batch_ids = [int(re.search('batches/isONbatch_(.*)\\.cer$', x).group(1))
                 for x in batches]
    LEVELS = OrderedDict()
    LEVELS[0] = []
    for Id, bf in sorted(zip(batch_ids, batches), key=lambda x: x[0]):
        n = Node(
            Id,
            "clusters/isONcluster_{}.cer".format(Id),
            None,
            None,
            None,
            0)
        n.Done = True
        JOB_TREE[Id] = n
        LEVELS[0].append(n)
    level = 0
    max_id = LEVELS[0][-1].Id
    while len(LEVELS[level]) != 1:  # Final level will be link
        next_level = level + 1
        LEVELS[next_level] = []
        for l_, r in grouper(2, LEVELS[level]):
            if r is None:  # End of a level
                LEVELS[level].pop()  # remove last node?
                l_.Level += 1  # ncrement level
                LEVELS[next_level].append(l_)  # Add the left to the next level
                continue
            max_id += 1
            new_batch = "clusters/isONcluster_{}.cer".format(max_id)
            new_node = Node(max_id, new_batch, l_, r, None, next_level)
            l_.Parent = new_node
            r.Parent = new_node
            r.RightSide = True
            LEVELS[next_level].append(new_node)
            JOB_TREE[max_id] = new_node
        level = next_level
    ROOT = JOB_TREE[len(JOB_TREE) - 1].Id
    JOB_TREE[ROOT].RightSide = True

    return JOB_TREE, LEVELS


def main():
    """Entry point."""
    Path('clusters').mkdir(exist_ok=True)
    job_tree, levels = build_job_tree()

    init_template = 'isONclust2 cluster -x {} -v -Q -l batches/' \
                    'isONbatch_{}.cer -o clusters/isONcluster_{}.cer {}; ' \
                    'sync;\n'
    template = 'isONclust2 cluster -x {} -v -Q -l clusters/isONcluster_{}' \
               '.cer -r clusters/isONcluster_{}.cer -o clusters/isONcluster' \
               '_{}.cer {}; sync\n'

    for nr, l in levels.items():
        jobs_out = 'jobs_level_{}.sh'.format(nr)
        with open(jobs_out, 'w') as fh:
            for n in l:
                purge = "-z" if n.RightSide else ""
                if nr == 0 or n.Left is None or n.Right is None:
                    jr = init_template.format('sahlin', n.Id, n.Id, purge)
                    fh.write(jr)
                else:
                    jr = template.format('sahlin', n.Left.Id,
                                         n.Right.Id, n.Id, purge)
                    fh.write(jr)
        # Run a level in parallel
        cmd = "parallel < {}".format(jobs_out)
        sub.call(cmd, shell=True)
    sub.call("ln -s `realpath clusters/isONcluster_{}"
             ".cer` isONcluster_ROOT.cer".format(n.Id), shell=True)


if __name__ == '__main__':
    # The cwd should be the process dir that contains 'batches/'
    main()
