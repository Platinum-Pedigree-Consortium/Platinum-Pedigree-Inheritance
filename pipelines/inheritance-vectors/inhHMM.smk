import pandas as pd
import os

chromosomes=["chr" + str(chrom) for chrom in  list(range(1,23))] + ["chrX"]

parents=[config["dad"], config["mom"]]
children=config["children"].split(",")
male_children=config["male_children"].split(",")

split_dict={
    "chr1": [*range(1, 9)],
    "chr2": [*range(1, 9)],
    "chr3": [*range(1, 8)],
    "chr4": [*range(1, 8)],
    "chr5": [*range(1, 7)],
    "chr6": [*range(1, 8)],
    "chr7": [*range(1, 7)],
    "chr8": [*range(1, 7)],
    "chr9": [*range(1, 6)],
    "chr10": [*range(1, 6)],
    "chr11": [*range(1, 6)],
    "chr12": [*range(1, 6)],
    "chr13": [*range(1, 5)],
    "chr14": [*range(1, 4)],
    "chr15": [*range(1, 4)],
    "chr16": [*range(1, 4)],
    "chr17": [*range(1, 4)],
    "chr18": [*range(1, 4)],
    "chr19": [*range(1, 3)],
    "chr20": [*range(1, 3)],
    "chr21": [*range(1, 3)],
    "chr22": [*range(1, 3)],
    "chrX": [*range(1, 4)]
}
chrom_keys, split_values = zip(*split_dict.items())

def chromosome_split_list(chromosomes, split_dict):
    outlist=[]
    for chrom in chromosomes:
        items=split_dict[chrom]
        for split in items:
            outlist.append(f"{chrom}.0{split}")
    return outlist



try:
    config["chromosomes"]
    chromosomes = [config["chromosomes"]]
except NameError:
    pass

outlist = chromosome_split_list(chromosomes, split_dict)

# viterbi params
snp_punishment="0.05"
change_punishment="0.8"

wildcard_constraints:
    parent="|".join(parents),
    chrom="|".join(chromosomes)


file_prefix = ".".join(os.path.basename(config["input_vcf"]).split('.')[:-2])

rule all:
    input:
        expand("output/split/{file_prefix}.{chroms}.{parents}.bed", file_prefix=file_prefix, chroms=chromosomes, parents=parents),
        expand("output/windows/viterbi/{file_prefix}.{chrom}.viterbi.csv", chrom=chromosomes, file_prefix=file_prefix),
        f"output/viterbi/{file_prefix}.inht_vectors.csv"

rule prepare_snps:
    params:
        mom=config["mom"],
        dad=config["dad"],
        children=",".join(children),
    input:
        script="../../code/inheritance_vectors/prepare_snps.py",
        cohort_ped="data/CEPH1463.inht_vectors.ped",
        cohort_vcf=config["input_vcf"]
    output:
        filtered=f"output/prepare_snps/filtered/{file_prefix}.filtered.vcf",
        bed=f"output/prepare_snps/{file_prefix}.bed"
    conda:
        "envs/vcfpy.yaml"
    shell:
        """
        python {input.script} \\
            --cohort-calls {input.cohort_vcf} \\
            --pedigree {input.cohort_ped} \\
            --dad-sample {params.dad} \\
            --mom-sample {params.mom} \\
            --subset-children {params.children} \\
            --filtered {output.filtered} \\
            > {output.bed}
        """

rule split_by_chrom:
    params:
        dad=config["dad"],
        mom=config["mom"]
    input:
        sites=f"output/prepare_snps/{file_prefix}.bed"
    output:
        dad=f"output/split/{file_prefix}.{{chrom}}.{config["dad"]}.bed",
        mom=f"output/split/{file_prefix}.{{chrom}}.{config["mom"]}.bed",
    conda:
        "envs/gnugrep.yaml"
    shell:
        """
        head -n1 {input.sites} > {output.dad}
        grep -P '^{wildcards.chrom}\t' {input.sites} | grep '{params.dad}' >> {output.dad}

        head -n1 {input.sites} > {output.mom}
        grep -P '^{wildcards.chrom}\t' {input.sites} | grep '{params.mom}' >> {output.mom}
        """


rule split_sites_per_chrom:
    params:
        splits=lambda wildcards: split_dict[(wildcards.chrom)][-1]
    input:
        sites="output/split/{file_prefix}.{chrom}.{parent}.bed"
    output:
        temp_sites="output/split/{file_prefix}.{chrom}.{parent}.temp",
        split_dir=directory("output/split/{file_prefix}.{chrom}.{parent}"),
        #split_files=expand("output/split/{{chrom}}.{{parent}}/{{chrom}}.{{parent}}.0{split}", split = list(range(1, params.splits)))
    conda:
        "envs/coreutils.yaml"
    shell:
        """
        tail -n +2 {input.sites} > {output.temp_sites}
        split -n l/{params.splits} --numeric-suffixes=1 {output.temp_sites} "{wildcards.chrom}.{wildcards.parent}."
        mkdir -p {output.split_dir}
        find "." -maxdepth 1 -name "{wildcards.chrom}.{wildcards.parent}.*" -exec mv -t "{output.split_dir}/" {{}} \;
        for file in "{output.split_dir}/*"; do
            sed -i '1s/^/CHROM\\tstart\\tend\\tREF\\tALT\\tcalled_parent\\tgrandparent\\tphase\\tchildren_calls\\n/' $file
        done
        """

rule viterbi_window:
    params:
        test_outdir="output/windows/viterbi_{parent}_{chrom}.0{split}",
        snp_punishment=snp_punishment,
        change_punishment=change_punishment,
        children=";".join(children),
        male_children=";".join(male_children),
        parents=",".join(parents)
    input:
        split_dir=directory("output/split/{file_prefix}.{chrom}.{parent}"),
        #sites="output/split/{chrom}.{parent}/{chrom}.{parent}.0{split}",
        script="../../code/inheritance_vectors/viterbi.py"
    output:
        t_mat="output/windows/viterbi_{parent}_{chrom}.0{split}/{file_prefix}.{chrom}_{parent}.t_matrix.txt",
        e_mat="output/windows/viterbi_{parent}_{chrom}.0{split}/{file_prefix}.{chrom}_{parent}.e_matrix.txt",
        viterbi="output/windows/viterbi_{parent}_{chrom}.0{split}/{file_prefix}.{chrom}_{parent}.out.tsv",
        out_df="output/windows/viterbi_{parent}_{chrom}.0{split}/{file_prefix}.{chrom}_{parent}.viterbi_df.txt"
    conda:
        "envs/inht_vectors.yaml"
    shell:
        """
        python3 {input.script} \\
            --input "{input.split_dir}/{wildcards.chrom}.{wildcards.parent}.0{wildcards.split}" \\
            --file-prefix {wildcards.file_prefix} \\
            --parents-list "{params.parents}" \\
            --children "{params.children}" \\
            --male-children "{params.male_children}" \\
            --transmission-matrix {output.t_mat} \\
            --emission-matrix {output.e_mat} \\
            --test-outdir {params.test_outdir} \\
            --punishment "{params.snp_punishment},{params.change_punishment}" \\
            --output {output.viterbi}
        """

rule gather_window_results:
    input:
        chr1="output/windows/viterbi_{parent}_{chrom}.01/{file_prefix}.{chrom}_{parent}.viterbi_df.txt",
        out_files=expand("output/windows/viterbi_{{parent}}_{splits_ending}/{{file_prefix}}.{{chrom}}_{{parent}}.out.tsv", splits_ending=outlist)
    output:
        "output/windows/viterbi/{file_prefix}.{chrom}.{parent}.viterbi.tsv"
    conda:
        "envs/gnugrep.yaml"
    shell:
        """
        c1grep() {{ grep "$@" || test $? = 1; }}
        head -n 1 {input.chr1} > {output}

        OUT_DIRS=($(dirname {input.out_files}))
        for i in ${{OUT_DIRS[@]}}; do 
            cat ${{i}}/chr*_NA*.viterbi_df.txt | c1grep -v "CHROM" | c1grep -P "{wildcards.chrom}\\t" >> {output}
        done
        """

rule viterbi_to_inht_vectors:
    params:
        parents=",".join(parents),
        children=";".join(children),
        male_children=";".join(male_children)
    input:
        script="../../code/inheritance_vectors/viterbi_to_inht_vectors.py",
        dad_sites=lambda wildcards: f"output/windows/viterbi/{wildcards.file_prefix}.{wildcards.chrom}.{config["dad"]}.viterbi.tsv",
        mom_sites=lambda wildcards: f"output/windows/viterbi/{wildcards.file_prefix}.{wildcards.chrom}.{config["mom"]}.viterbi.tsv"
    output:
        inht_vectors="output/windows/viterbi/{file_prefix}.{chrom}.viterbi.csv",
        chr_summary="output/windows/viterbi/{file_prefix}.{chrom}.summary.tsv",
        gaps=expand("output/windows/viterbi/{{file_prefix}}.{{chrom}}.{child}.tsv", child=children)
    conda:
        "envs/pandas.yaml"
    shell:
        """
        python3 {input.script} \\
            --input-dad {input.dad_sites} \\
            --input-mom {input.mom_sites} \\
            --children "{params.children}" \\
            --male-children "{params.male_children}" \\
            --parents "{params.parents}" \\
            --output-inht {output.inht_vectors} \\
            --output-gaps "output/windows/viterbi/{wildcards.file_prefix}.{wildcards.chrom}" \\
            --output-analysis "output/windows/viterbi/{wildcards.file_prefix}.{wildcards.chrom}.test_analysis.txt" \\
            --output-summary {output.chr_summary}
        """

rule combine_vectors:
    input:
        expand("output/windows/viterbi/{file_prefix}.{chrom}.viterbi.csv", chrom=chromosomes, file_prefix=file_prefix)
    output:
        f"output/viterbi/{file_prefix}.inht_vectors.csv"
    shell:
        """
        head -n 1 {input[0]} > {output}
        cat {input} | grep -v "CHROM" >> {output}
        """
