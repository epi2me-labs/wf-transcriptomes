*Does the workflow support de novo assembly?* - Currently the workflow does not have a *de novo* mode.

*Why is the IGV panel not showing?* - The workflow expects either an uncompressed or [`bgzip`](https://www.htslib.org/doc/bgzip.html)-compressed reference. If the user provides a reference compressed not with `bgzip`, the workflow will run to completion, but won't be able to generate the necessary indexes to visualize the outputs in IGV.

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-transcriptomes/issues) page or start a discussion on the [community](https://community.nanoporetech.com/).