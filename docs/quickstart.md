## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop),
[Singularity](https://sylabs.io/singularity/) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Each method is automated out-of-the-box provided
either docker, singularity or conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).


**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-isoforms --help
```

to see the options for the workflow.



**Example execution of a workflow for reference-based transcript assembly**

This example uses a synthetic SIRV dataset, so we need to tell minimap2 about the non-canonical splice junctions with 
--minimap2_opts '-uf --splice-flank=no'
```
OUTPUT=~/output;
nextflow run wf-isoforms/ --fastq test_data/fastq  --ref_genome test_data/SIRV_150601a.fasta --ref_annotation test_data/SIRV_isofroms.gtf
--minimap2_opts '-uf --splice-flank=no' --out_dir outdir -w workspace_dir -profile conda -resume
```

```
# To evaluate the workflow on a larger Drosophila dataset
./evaluation/run_evaluation_dmel.sh outdir
```

**Example workflow for denovo transcript assembly**
```
OUTPUT=~/output
nextflow run . --fastq test_data/fastq --denovo --ref_genome test_data/SIRV_150601a.fasta  -profile local --out_dir ${OUTPUT} -w ${OUTPUT}/workspace \
--sample sample_id -resume
```
A full list of options can be seen in nextflow_schema.json. Below are some commonly used ones.

- Threshold for including isoforms into interactive table `transcript_table_cov_thresh = 50`
- Run the denovo pipeline `denovo = true` (default false)
- To run the workflow with direct RNA reads `--direct_rna` (skips the pychopper step).


Pychopper and minimap2 can take options via `minimap2_opts` and `pychopper_opts`, for example:


- When using the SIRV synthetic test data  
  - `minimap2_opts = '-uf --splice-flank=no'`
- pychopper needs to know which cDNA synthesis kit used
  - SQK-PCS109: use `pychopper_opts = '-k PCS109'` (default)
  - SQK-PCS110: use `pychopper_opts = '-k PCS110'`
- pychopper can use one of two available backends for identifying primers in the raw reads
  - nhmmscan `pychopper opts = '-m phmm'` 
  - edlib `pychopper opts = '-m edlib'`

__Note__: edlib is set by default in the config as it's quite a lot faster. However it may be less sensitive than nhmmscan. 
  
