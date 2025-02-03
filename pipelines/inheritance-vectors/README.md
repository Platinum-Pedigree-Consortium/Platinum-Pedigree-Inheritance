# Inheritance vectors

## Workflow steps

1. Prepare SNPs (`prepare_snps.py`)

This step assigns haplotype identity for each SNP in the cohort vcf, by taking the identity from the grandparents, as defined in the pedigree file. Haplotypes are named arbitrarily as A, B, C, D (paternal grandfather - A, paternal grandmother - B, maternal grandfather - C, maternal grandmother - D).

It is run as follows:

```bash
python prepare_snps.py \
    --cohort-calls {input.cohort_vcf} \
    --pedigree {input.cohort_ped} \
    --dad-sample {params.dad} \
    --mom-sample {params.mom} \
    --subset-children {params.children} \
    --filtered {output.filtered} \
    > {output.bed}
```

It ouputs a BED file with the following structure:

```
| #CHROM | start | end  | REF | ALT | called_parent | grandparent | phase | children_calls |
|--------|-------|------|-----|-----|---------------|-------------|-------|----------------|
| chr1   | 4736  | 4736 | C   | A   | NA12886       | NA12878     | B     | 200102         |
| chr1   | 5324  | 5324 | C   | A   | NA12886       | NA12878     | B     | 200102         |
| chr1   | 6137  | 6137 | G   | A   | NA12886       | NA12878     | B     | 200102         |
```

Assigned haplotype is contained in the `phase` column, the children in which this variant is observed are listed in the `children_calls` column.

Output bed files are split into seperate files for variants occuring in each parent, these bed files are then split by chromosome and into windows along each chromosome to allow parallelisation by window.

2. Viterbi (`viterbi.py`)

The HMM is defined according to the "snp_punishment" & "change_punishment" parameters. These create the transmission and emission matrices, which describe the probability of a SNP calling error or a change in inheritance (recombination) respectively. The viterbi algorithm is run per window, calculating the most likely haplotype blocks in each window.

```bash
python3 viterbi.py \
    --input {split_dir} \
    --file-prefix {file_prefix} \
    --parents-list {parents} \
    --children {children} \
    --male-children {params.male_children} \
    --transmission-matrix {t_mat} \
    --emission-matrix {e_mat} \
    --test-outdir {outdir} \
    --punishment "{snp_punishment},{change_punishment}" \
    --output {viterbi}
```

Per window, per parent viterbi file is output to `output/windows/viterbi_{parent}_{chrom}.0{split}/{file_prefix}.{chrom}_{parent}.viterbi_df.txt`

| Row | CHROM | start   | end     | REF | ALT | called_parent | grandparent | phase | children_calls            | x_seq_opt | x_seq_opt_state                   |
|-----|-------|---------|---------|-----|-----|---------------|-------------|-------|---------------------------|-----------|-----------------------------------|
| 788 | chr1  | 2194656 | 2194656 | C   | T   | NA12879       | NA12877     | C     | 200084;200086;200087      | 40.0      | [('200084', '200086', '200087')] |
| 789 | chr1  | 2195012 | 2195012 | G   | A   | NA12879       | NA12877     | C     | 200084;200086;200087      | 40.0      | [('200084', '200086', '200087')] |
| 790 | chr1  | 2199387 | 2199387 | C   | G   | NA12879       | NA12878     | D     | 200081;200082;200085      | 40.0      | [('200084', '200086', '200087')] |
| 791 | chr1  | 2241520 | 2241520 | A   | G   | NA12879       | NA12877     | C     | 200084;200086;200087      | 40.0      | [('200084', '200086', '200087')] |

This file outputs the most likely inheriting children for the first grandparental haplotype (A or C) at each genomic position in the window, contained in the `x_seq_opt_state` column. 

Optimum states are then collapsed into continuous haplotype blocks, and windowed files are combined into a single output inheritance vectors file (`output/viterbi/{vcf_prefix}.inht_vectors.csv`):

| CHROM | start    | end      | NA12879 | NA12881 | NA12882 | NA12883 | NA12884 | NA12885 | NA12886 | NA12887 | NA12877 | NA12878 |
|-------|----------|----------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|
| chr1  | 5324     | 996948   | AC      | BD      | BC      | AC      | AD      | BC      | BD      | BC      | AB      | CD      |
| chr1  | 1027469  | 2648764  | AC      | BD      | BC      | AC      | AD      | AC      | BD      | BC      | AB      | CD      |
| chr1  | 2690284  | 4613470  | AC      | BD      | BC      | AC      | AD      | AC      | AD      | BC      | AB      | CD      |

## How to run the snakemake workflow

Set up and activate your snakemake environment using conda to run the workflow:

```bash
conda create -n snakemake -f envs/snakemake.yaml
conda activate snakemake
```

Run a test dataset using just chromosome 22:

```bash
VCF="data/chr22.dv.all.vcf.gz"

snakemake \
  -s inhHMM.smk \
  --config \
  "input_vcf=${VCF}" \
  "mom=NA12878" \
  "dad=NA12877" \
  "children=NA12879,NA12881,NA12882,NA12883,NA12884,NA12885,NA12886,NA12887" \
  "male_children=NA12882,NA12883,NA12884,NA12886" \
  "chromosomes=chr22" \
  --use-conda \
  -p
```

Final inheritance vectors will be output to `output/viterbi/{vcf_prefix}.inht_vectors.csv`

To run on the entire genome just remove the `chromosome` parameter from the above command, and update your vcf path, ie:

```bash
VCF="data/cohort.vcf.gz"

snakemake \
  -s inhHMM.smk \
  --config \
  "input_vcf=${VCF}" \
  "mom=NA12878" \
  "dad=NA12877" \
  "children=NA12879,NA12881,NA12882,NA12883,NA12884,NA12885,NA12886,NA12887" \
  "male_children=NA12882,NA12883,NA12884,NA12886" \
  --use-conda \
  -p
```