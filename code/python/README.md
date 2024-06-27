# SV merging

## VCF normalization

Prior to SV merging it is important to handle the differences of SV callers by normalizing the VCF files. Most notably this includes dropping INFO/FORMAT fields that are not shared across different callers, reformatting fields, and prefixing SV caller id to variant IDs. This can be done using `strip_vcf.py`, given an input path containing VCF files (*.pedfilt.vcf files) and an output path. Which will write a normalized VCF file (*.strip.pedfilt.vcf) for each input.

Supported SV callers: sawfish, sniffles, pav, pggb, pbsv

```
python strip_vcf.py -i ../sv_vcfs/in/ -o ../sv_vcfs/out/
```

## VCF merging

After normalization the VCF files can be merged using `interval_tree_merge.py` given the path of normalized VCFs, an output file, and the merging order. Optional CLI `--flank_len` (default = 200), allows for setting the size of flanks when building intervals.

The below command (assuming each SV VCF file is present), first merge pbsv into sawfish, then sniffles in that result, and finally do the same with PAV.

```
python interval_tree_merge.py --sv_path ../sv_vcfs/out/ -o merged.vcf --sv_order Sawfish,pbsv,Sniffles,PAV
```