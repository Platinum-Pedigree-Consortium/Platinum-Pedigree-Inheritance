# Methods for generating High Confidence Regions (HCR) for G2/G3


### Get HiFi and ONT Depth against GRCh38 (per sample per tech)
```
 mosdepth -t 2 --no-per-base --fast-mode -b 100 -Q 10 {sample_name} {bam}
```

### Get 10-fold or better coverage (per sample per tech)
```
 zcat {sample_name}.regions.bed.gz | perl -lane 'print if $F[-1] > 9' > {sample_name}.10-fold.bed
```

### Get 10 fold coverage across all samples (per tech)
```
 bedtools multiinter -i *10-fold.bed.gz | perl -lane 'print if $F[3] == 10' > 10-fold-total.bed
```

### Merge (Union of two techs)

```
 cat ont/10-fold-total.bed hifi/10-fold-total.bed | bedtools sort -i - | bedtools merge -i - > HCR.bed
```