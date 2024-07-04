
These are instructions to install and run the workflow on command line.
You can also access the workflow via the
[EPI2ME Desktop application](https://labs.epi2me.io/downloads/).

The workflow uses [Nextflow](https://www.nextflow.io/) to manage
compute and software resources,
therefore Nextflow will need to be
installed before attempting to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop
or [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html)
to provide isolation of the required software.
Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.
This is controlled by the
[`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles)
parameter as exemplified below.

It is not required to clone or download the git repository
in order to run the workflow.
More information on running EPI2ME workflows can
be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow.
This will pull the repository in to the assets folder of
Nextflow and provide a list of all parameters
available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-transcriptomes --help
```
To update a workflow to the latest version on the command line use
the following command:
```
nextflow pull epi2me-labs/wf-transcriptomes
```
A demo dataset is provided for testing of the workflow.
It can be downloaded and unpacked using the following commands:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-transcriptomes/wf-transcriptomes-demo.tar.gz
tar -xzvf wf-transcriptomes-demo.tar.gz
```
The workflow can then be run with the downloaded demo data using:
```
nextflow run epi2me-labs/wf-transcriptomes \
	--de_analysis \
	--direct_rna \
	--fastq 'wf-transcriptomes-demo/differential_expression_fastq' \
	--jaffal_annotation 'genCode22' \
	--jaffal_genome 'hg38_chr20' \
	--jaffal_refBase 'wf-transcriptomes-demo/chr20' \
	--minimap2_index_opts '-k15' \
	--ref_annotation 'wf-transcriptomes-demo/gencode.v22.annotation.chr20.gtf' \
	--ref_genome 'wf-transcriptomes-demo/hg38_chr20.fa' \
	--sample_sheet 'wf-transcriptomes-demo/sample_sheet.csv' \
	-profile standard
```
For further information about running a workflow on
the command line see https://labs.epi2me.io/wfquickstart/
