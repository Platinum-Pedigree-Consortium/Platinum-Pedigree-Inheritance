## Verification


### Verify simple insertion case

```vcf
chr1	30515560	.	N	<X:1>	30	PASS	END=30515560;AN=20;NS=20;NA=2;ALEN=0,144;AC=4;VS=>s1138;VE=>s1139;AWALK=*,<s89506	GT:GT0	0:0	0:0	1:1	0:0	1:1	0:0	0:0	0:0	1:1	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	1:1	0:0
```

Ref has `*` as AWALK, and the samples match in the bam at `code/rust/notes/igv.xml`.

The lengths are 144 or 145, which matches expected `ALEN`.

### Verify simple deletion case


## Remaining issues

1. Properly represent indels per VCF specification.
