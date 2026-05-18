### Does the workflow support both cDNA and direct RNA?

Yes. Use `--direct_rna` for direct RNA data. cDNA is the default mode.

### Do I need both `--ref_genome` and `--ref_annotation` in fixed-annotation mode?

Yes. `bambu` still uses the genome together with the imported annotation.

### Does the workflow create one shared transcriptome or one per sample?

Both. The joint cohort model is the primary result for reporting and DE/DTU, and
the workflow also emits independent per-sample transcriptomes.

### Can I run DE/DTU without transcript discovery?

Yes. Use `--transcriptome_mode fixed_annotation` together with `--de_analysis`.

### What changed from the previous workflow version?

The current workflow uses `bambu`, optional `SQANTI3`, `DESeq2`, and `DEXSeq`.
The main transcriptome result is now one shared `bambu` model built from all
samples together, with separate per-sample `bambu` outputs published alongside
it. The most important differences are summarised below.

#### Why do the results now mention `bambu` and `SQANTI3` instead of StringTie or GffCompare?

The workflow now uses a different set of transcript analysis tools. It builds
its transcript models with `bambu`, optionally classifies them with `SQANTI3`,
and performs DE/DTU from the cohort `bambu` outputs.

#### Why is `--ref_annotation` now required even in fixed-annotation mode?

`bambu` still requires the annotation together with the genome in both
`discover` and `fixed_annotation` modes. Fixed-annotation mode means
annotation-driven quantification, not annotation-only execution.

#### Can I still use `--ref_transcriptome`?

No. `--ref_transcriptome` has been removed from the workflow interface.

Short old-to-new example:

```text
Previous workflow version:
  --transcriptome_source precomputed --ref_transcriptome transcripts.fa

Current workflow:
  --transcriptome_mode fixed_annotation \
  --ref_genome genome.fa \
  --ref_annotation annotation.gtf
```

#### Why do I now get both cohort and per-sample transcriptomes?

The workflow now treats both as important outputs. The shared cohort model is
the main transcriptome used for reporting and optional DE/DTU, while the
per-sample transcriptomes are provided for looking at each sample separately.

#### Why are DE results under `de_analysis/<contrast>/`?

The workflow now writes one subdirectory per comparison instead of publishing
one flat DE result set. This makes the output folder clearer when more than one
comparison is present.

Short old-to-new example:

```text
Previous workflow version:
  de_analysis/results_dge.tsv

Current workflow:
  de_analysis/<contrast>/results_dge.tsv
```

#### Can I still run fixed-annotation quantification without transcript discovery?

Yes. Use `--transcriptome_mode fixed_annotation` together with `--ref_genome`,
`--ref_annotation`, and any optional DE/DTU settings.

#### What should I expect to differ in report contents and output filenames?

Expect the report and output folder to emphasise:

+ the joint cohort `bambu` transcriptome under `cohort/`
+ the per-sample `bambu` transcriptomes under `samples/<alias>/`
+ optional `SQANTI3` results under cohort and per-sample directories
+ contrast-specific DE/DTU outputs under `de_analysis/<contrast>/`

If your question is not answered here, please report issues or suggestions on
the [GitHub issues](https://github.com/epi2me-labs/wf-transcriptomes/issues)
page or start a discussion on the
[community](https://community.nanoporetech.com/).
