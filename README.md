# Polygenic liability to bipolar disorder and the plasma proteome (UK Biobank)

> **Work in progress.** A manuscript based on this project is in preparation with [Dr Rona Strawbridge](https://www.gla.ac.uk/schools/healthwellbeing/staff/ronastrawbridge/)
> (University of Glasgow). Specific results are omitted here until publication. No licence is currently
> granted; please contact the authors before any reuse.

This project comes from my MSc dissertation in Medical Genetics and Genomics at the University of Glasgow
(June–August 2025, supervised by [Dr Rona Strawbridge](https://www.gla.ac.uk/schools/healthwellbeing/staff/ronastrawbridge/)). It asks whether polygenic liability to bipolar disorder
is associated with circulating plasma protein levels in UK Biobank, and whether any such associations remain
after adjustment for demographic, technical, lifestyle, cardiometabolic and psychiatric factors.

A two-page overview is here:
[Research Project Overview (PDF)](https://github.com/bleulievre/bleulievre/blob/main/portfolio/Research_Project_Overview.pdf).

## Contents

| Folder | What it is |
| --- | --- |
| [`thesis_code/`](thesis_code/) | A public copy of the R script I wrote for the dissertation, in the order the analysis was done |
| [`pipeline/`](pipeline/) | A separate re-implementation of the analysis, written after the MSc with AI coding assistance. It is driven by configuration files, has unit and end-to-end tests, and runs on simulated data |

### `thesis_code/UKB.R`

The dissertation analysis: data cleaning, protein missingness profiling, Spearman screening, four stages of
multivariable regression with Bonferroni thresholds, a sensitivity analysis that excludes participants with
any history of mental illness or psychiatric medication, and Gene Ontology enrichment with clusterProfiler.
The section headers follow the weeks of the project.

This public copy differs from the script used for the dissertation in three ways:
- comments are translated into English;
- comments that stated intermediate results are replaced with neutral ones;
- four lines that hard-coded the final protein names now take them from the model results instead (marked
  `PUBLIC COPY`).

Apart from these lines, the R code is unchanged. It needs the UK Biobank extract used in the dissertation, so
it does not run without approved data access.

### `pipeline/`

After the MSc, the dissertation analysis was re-implemented as a tested pipeline, with AI coding assistance,
in preparation for the manuscript:
- exclusions, models and thresholds are set in `config/analysis_plan.yaml`, which is still a draft;
- it can run the original staged procedure (`R/legacy.R`), following the dissertation code step by step, as
  well as draft model specifications for the manuscript;
- it is designed to write aggregate tables only and includes a check that blocks the export of participant
  identifiers. This is not a substitute for UK Biobank's own output checks;
- it ships with a simulator (`synthetic/make_synthetic.R`), so it has been run and tested on simulated data
  only.

See [`pipeline/README.md`](pipeline/README.md) for details.

```
cd pipeline
Rscript synthetic/make_synthetic.R
Rscript run_all.R --synthetic --layout glasgow
Rscript tests/run_tests.R
```

## Data

No participant-level UK Biobank data, protein identities, effect estimates or findings are included. Access
to UK Biobank data is controlled by UK Biobank and limited to authorised researchers working on approved
projects (https://www.ukbiobank.ac.uk/). The polygenic scores used were LDpred scores from published
Psychiatric Genomics Consortium summary statistics, provided by the research group.

## Contact

Yunqi (Vicky) Liu — vickyliu1998424@gmail.com
