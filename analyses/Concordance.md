# Pedigree filtering and phasing using inheritance vectors

 The concordance program accepts a VCF and the inheritance vectors in the data folder. The concordance code generates all possible combinations of valid genotypes over a haplotype (region) and tests if the observed genotypes at a site match a valid combination. Sites that pass the inheritance filtering are kept if they contain a genotype call for every individual (excluding any site with nocalls). Using the inheritance vectors we can directly phase genotypes for bi-allelic and multi-allelic sites. There is a bijective mapping of observed genotypes in the VCF with the haplotype structure with the exception of all heterozygous genotypes, where the phase is ambiguous. We exclude all heterozygous sites from the truth set. The command used to filter and phase the variants is:

building:
```
1. checkout source code
2. cd code/rust ; cargo build --release
3. built code should be found in the generated folder: target/release/concordance
```


usage:
```
concordance \
--father NA12877 \
--mother NA12878 \
--vcf data/ceph_oa_test.vcf.gz \
--inheritance ceph.GRCh38.viterbi.oa.csv.gz \
--prefix test > std.out
```

runtime for test data:
```
real    0m0.103s
user    0m0.073s
sys     0m0.012s
```


The open access inheritance vectors are distributed with this repo.

