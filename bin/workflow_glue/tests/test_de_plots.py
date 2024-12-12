"""Test assign_barcodes."""
from pathlib import Path

import pytest
from workflow_glue.de_plots import get_translations


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return Path(request.config.getoption("--test_data")) / "workflow_glue"


@pytest.mark.parametrize(
    'annotation_file,expected',
    [
        [
            'MSTRG.11088.gtf',
            dict(gid_to_gene_name={
                'ENSG00000236051.7': 'MYCBP2-AS1',
                'ENSG00000283208.2': 'ENSG00000283208',
                'ENSG00000102805.16': 'CLN5',
                'MSTRG.11088': 'MSTRG.11088'
            },
                txid_to_gene_name={
                    'ENST00000636183.2': 'CLN5',
                    'ENST00000636780.2': 'CLN5',
                    'ENST00000638147.2': 'ENSG00000283208',
                    'ENST00000637192.1': 'ENSG00000283208',
                    'ENST00000636737.1': 'MYCBP2-AS1',
                    'ENST00000450627.6': 'MYCBP2-AS1',
                    'MSTRG.11088.2': 'MSTRG.11088'
                },
                txid_to_gene_id={
                    'ENST00000636183.2': 'ENSG00000102805.16',
                    'MSTRG.11088.2': 'MSTRG.11088',
                    'ENST00000636780.2': 'ENSG00000102805.16',
                    'ENST00000638147.2': 'ENSG00000283208.2',
                    'ENST00000637192.1': 'ENSG00000283208.2',
                    'ENST00000636737.1': 'ENSG00000236051.7',
                    'ENST00000450627.6': 'ENSG00000236051.7'
                })

        ],
        # Small test to check that GFF3 works
        [
            'MSTRG.11088.gff3',
            dict(gid_to_gene_name={
                "ENSG00000290825.1": "DDX11L2",
                "ENSG00000236397.3": "DDX11L2"
            },
                txid_to_gene_name={
                    "ENST00000456328.2": "DDX11L2",
                    "ENST00000437401.1": "DDX11L2"
                },
                txid_to_gene_id={
                    'ENST00000437401.1': 'ENSG00000236397.3',
                    'ENST00000456328.2': 'ENSG00000290825.1'
                })
        ]
    ]
)
def test_get_translations(test_data, annotation_file, expected):
    """Test that correct feature identifiers are extracted from the annotation.

    `stringtie --merge` can sometimes generate gene models that may span multiple
    reference genes. Possibly related issue:
    https://github.com/gpertea/stringtie/issues/217
    This can lead to the original genes and transcripts being assigned to that
    incorrectly-merged gene model. The test data contains such a gene model generated
    from `stringtie --merge` but actually consists of multiple different genes.


    """
    input_gtf = test_data / annotation_file
    txid_to_gene_name, txid_to_gene_id, gid_to_gene_name = get_translations(input_gtf)

    assert expected['gid_to_gene_name'] == gid_to_gene_name
    assert expected['txid_to_gene_name'] == txid_to_gene_name
    assert expected['txid_to_gene_id'] == txid_to_gene_id
