
# Specific info about the King County build
- The `builds.yaml` for the variant specific King County build can be found under `my_profiles/variant_build/` along with the associated config files. This can be run by using the Nextstrain CLI command as follows: `nextstrain build . --configfile my_profiles/variant_build/builds.yaml -p`
- In order to incorporate `Puma` metadata, a custom Snakemake rule (`download_metadata_with_puma.smk`) was added to the build which uses `scan_seq_puma_sourcename.tsv` obtained from a Metabase query and contains strain name and associated Puma and then adds it to the downloaded metadata from the Nextstrain S3 bucket using `scripts/puma_kc_seq_merge.py`


#Everything below is from the upstream Nextrain/ncov repository

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/nextstrain/ncov)](https://github.com/nextstrain/ncov/releases)
[![See recent changes](https://img.shields.io/badge/changelog-See%20recent%20changes-blue)](https://docs.nextstrain.org/projects/ncov/en/latest/reference/change_log.html)


# About

This repository analyzes viral genomes using [Nextstrain](https://nextstrain.org) to understand how SARS-CoV-2, the virus that is responsible for the COVID-19 pandemic, evolves and spreads.

We maintain a number of publicly-available builds, visible at [nextstrain.org/ncov](https://nextstrain.org/ncov).

[See our change log for details about backwards-incompatible or breaking changes to the workflow](https://docs.nextstrain.org/projects/ncov/en/latest/reference/change_log.html).

# Resources

## Use Nextstrain to analyze your SARS-CoV-2 data

**We've written a comprehensive guide to get you up and running in <1 hr. Click on the below links to follow it. It covers:**

0. [**Introduction**](https://nextstrain.github.io/ncov/index) _(Start here)_
1. [**Setup and installation**](https://nextstrain.github.io/ncov/setup)
2. [**Preparing your data**](https://nextstrain.github.io/ncov/data-prep)
3. [**Orientation: analysis workflow**](https://nextstrain.github.io/ncov/orientation-workflow)
4. [**Orientation: which files should I touch?**](https://nextstrain.github.io/ncov/orientation-files)
5. [**Running & troubleshooting**](https://nextstrain.github.io/ncov/running)
6. [**Customizing your analysis**](https://nextstrain.github.io/ncov/customizing-analysis) (see also: [reference for all configuration parameters](https://nextstrain.github.io/ncov/configuration))
7. [**Customizing your visualization**](https://nextstrain.github.io/ncov/customizing-visualization)
8. [**Options for visualizing and sharing results**](https://nextstrain.github.io/ncov/sharing) (including working with sensitive metadata)
9. [**Interpreting your results**](https://nextstrain.github.io/ncov/interpretation)
10. [**Writing a narrative to highlight key findings**](https://nextstrain.github.io/ncov/narratives)
11. [**Running the pipeline starting with multiple inputs**](https://nextstrain.github.io/ncov/multiple_inputs)

## Reference guides

  - [Metadata field definitions](https://docs.nextstrain.org/projects/ncov/en/latest/reference/metadata-fields.html)
  - [Workflow configuration parameters](https://docs.nextstrain.org/projects/ncov/en/latest/reference/configuration.html)

## Download formatted datasets

The hCoV-19 / SARS-CoV-2 genomes were generously shared via GISAID. We gratefully acknowledge the Authors, Originating and Submitting laboratories of the genetic sequence and metadata made available through GISAID on which this research is based.

In order to download the GISAID data to run the analysis yourself, please see [this guide](https://docs.nextstrain.org/projects/ncov/en/latest/analysis/data-prep.html).
> Please note that `data/metadata.tsv` is no longer included as part of this repo. However, we provide continually-updated, pre-formatted metadata & fasta files for download through GISAID.

## Read previous Situation Reports
We issued weekly Situation Reports for the first ~5 months of the pandemic. You can find the Reports and their translations [here](https://nextstrain.org/ncov-sit-reps).

## FAQs

- Can't find your sequences in Nextstrain? Check [here](./docs/data_faq.md) for common reasons why your sequences may not be appearing.
You can also use [clades.nextstrain.org](https://clades.nextstrain.org/) to perform some basic quality control on your sequences. If they are flagged by this tool, they will likely be excluded by our pipeline.
- For information about how clades are defined, and the currently named clades, please see [here](./docs/naming_clades.md). To assign clades to your own sequences, you can use our clade assignment tool at [clades.nextstrain.org](https://clades.nextstrain.org/).

## Bioinformatics notes

Site numbering and genome structure uses [Wuhan-Hu-1/2019](https://www.ncbi.nlm.nih.gov/nuccore/MN908947) as reference. The phylogeny is rooted relative to early samples from Wuhan. Temporal resolution assumes a nucleotide substitution rate of [8 &times; 10^-4 subs per site per year](http://virological.org/t/phylodynamic-analysis-176-genomes-6-mar-2020/356). There were SNPs present in the nCoV samples in the first and last few bases of the alignment that were masked as likely sequencing artifacts.

# Contributing

We welcome contributions from the community! Please note that we strictly adhere to the [Contributor Covenant Code of Conduct](https://github.com/nextstrain/.github/blob/master/CODE_OF_CONDUCT.md).

### Contributing to software or documentation
Please see our [Contributor Guide](https://github.com/nextstrain/.github/blob/master/CONTRIBUTING.md) to get started!

### Contributing data
**Please note that we automatically pick up any SARS-CoV-2 data that is submitted to GISAID.**

If you're a lab and you'd like to get started sequencing, please see:
* [Protocols from the ARTIC network](https://www.protocols.io/groups/artic/publications)
* [Funding opportunities for sequencing efforts](https://twitter.com/firefoxx66/status/1242147905768751106)
* Or, if these don't meet your needs, [get in touch](mailto:hello@nextstrain.org)

---

# Get in touch

To report a bug, error, or feature request, please [open an isssue](https://github.com/nextstrain/ncov/issues).

For questions, head over to the [discussion board](https://discussion.nextstrain.org/); we're happy to help!
