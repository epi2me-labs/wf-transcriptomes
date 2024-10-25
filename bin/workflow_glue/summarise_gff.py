"""Get summary statistics from GFF file."""
from collections import Counter
from pathlib import Path
import pickle

import gffutils
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("summ_gff")
    parser.add_argument(
        "gff",
        help="Report output file",
        type=Path)
    parser.add_argument(
        "sample_id",
        help="Output TSV file path")
    parser.add_argument(
        "out",
        default="gff_summary.tsv",
        help="Output TSV file path",
        type=Path)

    return parser


def main(args):
    """Entry point."""
    db = gffutils.create_db(
        str(args.gff), dbfn=':memory:', force=True, keep_order=True,
        merge_strategy='merge', sort_attribute_values=True
    )

    num_transcripts = db.count_features_of_type('transcript')
    num_genes = db.count_features_of_type('gene')

    transcript_lens = []
    exons_per_transcript = Counter()
    isoforms_per_gene = Counter()

    for gene in db.features_of_type('gene'):

        n_isos = len(list(db.children(gene, featuretype='transcript')))
        isoforms_per_gene[n_isos] += 1

        for transcript in db.children(
                gene, featuretype='transcript', order_by='start'):
            tr_len = 0
            exons = list(db.children(transcript, featuretype='exon'))
            if len(exons) == 0:
                continue
            exons_per_transcript[len(exons)] += 1
            for ex in exons:
                tr_len += abs(ex.end - ex.start)

            transcript_lens.append(tr_len)

    results = {
        'sample_id': args.sample_id,
        'summaries': {
            'Total genes': [num_genes],
            'Total transcripts': [num_transcripts],
            'Max trans. len': max(transcript_lens),
            'Min trans. len': min(transcript_lens)
        },
        'transcript_lengths': transcript_lens,
        'exons_per_transcript': exons_per_transcript,
        'isoforms_per_gene': isoforms_per_gene
    }

    with open(args.out, 'wb') as fh:
        pickle.dump(results, fh)
