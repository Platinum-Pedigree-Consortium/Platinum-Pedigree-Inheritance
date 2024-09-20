# Structural Variant Analysis 

After haplotype inheritance filtering we sequentially merged the four following callsets:
1. Sawfish - HiFi
2. Sniffles - ONT
3. PGGB - assembly
4. PAV - assembly

The merging pipeline requires python 3.10.6, it's easiest to setup a conda environment.

### Commands used to merge the different callers
```
 mkdir stripped_vcfs
 python strip_vcf.py -i merging_config.json -o stripped_vcfs
 python interval_tree_merge.py --diff_threshold 30 -i  stripped_vcfs/sawfish.strip.pedfilt.vcf stripped_vcfs/sniffles_ont.strip.pedfilt.vcf stripped_vcfs/pggb.strip.pedfilt.vcf stripped_vcfs/pav.strip.pedfilt.vcf -o merged_hg38.svs.vcf
```

