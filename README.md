# Platinum Pedigree 

This repository includes the code, data, and documentation for analyzing large pedigrees, such as CEPH-1463.

For easy haplotype mapping across a pedigree see: [README](HAPLOTYPING.md)

## Pipelines

This folder includes Snakemake files used to generate the truth sets, covering both small variants and structural variants (SVs).

## Analyses

This folder contains README files detailing the various analyses performed:
 - **Pedigree phasing filtering**: Concordance analysis code.
 - **De novo mutation discovery**: Analysis for Porubsky et al., 2024.
 - **High Confidence Regions (HCR)**: Definition and selection of HCR for the truth set.

## Code

 - **Rust**: Source code for the toolkit (an ensemble of binaries).
 - **Python**: Source code for various analyses, including SV merging.
