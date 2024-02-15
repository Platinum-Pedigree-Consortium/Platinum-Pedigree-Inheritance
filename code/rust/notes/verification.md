## Verification


### Verify simple insertion case

```vcf
chr1	30515560	.	N	<X:1>	30	PASS	END=30515560;AN=20;NS=20;NA=2;ALEN=0,144;AC=4;VS=>s1138;VE=>s1139;AWALK=*,<s89506	GT:GT0	0:0	0:0	1:1	0:0	1:1	0:0	0:0	0:0	1:1	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	1:1	0:0
```

Ref has `*` as AWALK, and the samples match in the bam at `code/rust/notes/igv.xml`.

The lengths are 144 or 145, which matches expected `ALEN`.

### Verify simple deletion case

```vcf
chr1	120812954	.	N	<X:1>	30	PASS	END=120812954;AN=8;NS=8;NA=2;ALEN=6311,0;AC=2;VS=>s2129;VE=>s2130;AWALK=<s92056,*	GT:GT0	.	1:1	.	0:0	0:0	.	0:0	.	.	.	0:0	1:1	0:0	.	.	.	.	.	0:0	.
```

Breakpoint I don't trust, at least not without understanding the sequence.
There are many reads punching into the delete region, but the coverage drop is significant for NA12878-hap2.

Could not verify this one, but it's enormous:
```
chr1	16655140	.	N	<X:1>	30	PASS	END=16655140;AN=20;NS=20;NA=2;ALEN=88126,0;AC=1;VS=>s895;VE=>s896;AWALK=<s94649,*	GT:GT0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	0:0	1:1
```



## Remaining issues

1. Properly represent indels per VCF specification.
