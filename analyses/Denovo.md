
## Running rust tool to find denovo small variant candidates for HiFi DeepVariant calls
```
 denovo --ped CEPH1463.ped --vcf CEPH-1463.joint.GRCh38.deepvariant.glnexus.vcf.gz --dp 10 --gq 20 --founders NA12878,NA12877 > deepvariant_results.txt
```

## Running rust tool to find denovo small variant candidates from ONT (R9) Clair calls
```
 denovo --ped CEPH1463.ped --vcf CEPH1463.clair3.glnexus.vcf.gz --dp 10 --gq 20 --founders NA12878,NA12877 --unknown 15 > clair_ont_results.txt
```

## Running rust tool to find denovo small variant candidates from Illumina Dragen calls
 denovo --ped CEPH1463.ped --vcf palladium.v4.2.4.grc38.multisample.joint_genotyped.hard-filtered.reheader.vcf.gz --dp 10 --gq 20 --founders NA12878,NA12877 > dragen_results.txt