## V1.1 truthset generation

### requirements
- vcflib
- bedtools
- aardvark
- snakemake

### input data

The input data can be found on Amazon S3:

```
s3://platinum-pedigree-data/truthset_v1.1/small-variant-input
```

For more information on public datasets see:

[platinum pedigree](https://github.com/Platinum-Pedigree-Consortium/Platinum-Pedigree-Datasets)

### Running the pipeline

1. Getting the data in place:
```
 aws s3 cp --no-sign-request s3://platinum-pedigree-data/truthset_v1.1/small-variant-input input --recusive 
 ```

3. Edit config json
- You need the path to the GRCh38 reference.
- You need to make sure that input data folder is in place and relative to the config paths.

2. Running the snakemake
```
 snakemake --cores 20 --configfile truthset-v1_1.json  -s truthset-v1_1.snakemake --verbose -p
```