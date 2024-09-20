#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import vcfpy


#Define usage for command and subcommands
parser = argparse.ArgumentParser(
    description="Script to annotate child samples with parental inheritance",
    usage="inheritance_vectors_prep.py <command> [<args>]")
parser.add_argument(
    "-c", "--cohort-calls", required=True, help="Cohort VCF",
    type=str)
parser.add_argument(
    "-p", "--pedigree", required=True, help="Pedigree file",
    type=str)
parser.add_argument(
    "-d", "--dad-sample", required=True, help="Dad sample id",
    type=str)
parser.add_argument(
    "-m", "--mom-sample", required=True, help="Mom sample id",
    type=str)
parser.add_argument(
    "-s", "--subset-children", required=False, help="Optional. Subset of children from pedigree to consider.",
    type=str)
parser.add_argument(
    "-g", "--subset-gparents", required=False, help="Optional. Subset of gparents from pedigree to consider.",
    type=str)
parser.add_argument(
    "-f", "--filtered", required=False, help="Variant file",
    type=str)

### SNP filtering
###

def filter_variant(record, mom_sample, dad_sample, annotate_filtered_vars, filtered, children, grandparents, autosome):
    if not record.call_for_sample[mom_sample].is_variant and not record.call_for_sample[dad_sample].is_variant:
        annotate_filtered_vars(record, filtered, "parent_neither")
        return None
    # elif record.call_for_sample[mom_sample].is_variant and record.call_for_sample[dad_sample].is_variant:
    #     annotate_filtered_vars(record, filtered, "parent_both")
    #    return None
    elif record.call_for_sample[mom_sample].is_variant:
        called_parent=mom_sample
        called_parent_gt=record.call_for_sample[mom_sample].gt_alleles
        other_parent_sample=dad_sample
    else:
        called_parent=dad_sample
        called_parent_gt=record.call_for_sample[dad_sample].gt_alleles
        other_parent_sample=mom_sample
    parent_gt=record.call_for_sample[called_parent].gt_alleles
    # Filtering variants not relevant for transmission vectors
    ## SNV
    if not record.is_snv():
        annotate_filtered_vars(record, filtered, "variant_indel")
        return None
    ## Called in parent
    if not record.call_for_sample[called_parent].is_variant:
        annotate_filtered_vars(record, filtered, "variant_filtered")
        return None
    # Filter any sites where we have no call aka "./."
    if "./." in [record.call_for_sample[sample].data['GT'] for sample in children + [dad_sample] + [mom_sample] + grandparents]:
        annotate_filtered_vars(record, filtered, "variant_no_call")
        return None
    ## Remove - Hom variant in autosome
    if record.CHROM in autosome and not record.call_for_sample[called_parent].is_het:
        annotate_filtered_vars(record, filtered, "variant_autosome_not_het")
        return None
    ## Remove - Hom variant in mom and in chrX
    elif called_parent == mom_sample and record.CHROM == "chrX" and not record.call_for_sample[called_parent].is_het:
        annotate_filtered_vars(record, filtered, "variant_mom_chrX_hom")
        return None
    ## Remove - Dad variant in chrY
    elif record.CHROM == "chrY":
        annotate_filtered_vars(record, filtered, "variant_dad_chrY")
        return None

        ## Not called in other parent
    if record.call_for_sample[other_parent_sample].is_variant:
        annotate_filtered_vars(record, filtered, "variant_called_other_parent")
        return None
    ## Not multi-allelic
    if not len(record.ALT) == 1:
        annotate_filtered_vars(record, filtered, "variant_multialleleic")
        return None
    return record, called_parent, other_parent_sample

def annotate_filtered_vars(record, filtered, reason):
    # header
    #"chrom", "start", 'REF', "ALT", "parent, ""calls", "genotype", "depth", "filter"
    filtered_line=[record.CHROM, record.POS, record.REF, reason]
    filtered.write("\t".join(str(x) for x in filtered_line) + "\n")

    
class parent_info:
    """Required info regarding mom or dad"""
    # If all parental haplotypes are known (ie have all four grandparents)
    known_haplotypes={
        "A": "B",
        "B": "A",
        "C": "D",
        "D": "C"
    }
    # If not all parental haplotypes are known (ie don't have all grandparents)
    inferred_haplotypes={
        "A": "B_inferred",
        "B": "A_inferred",
        "C": "D_inferred",
        "D": "C_inferred",
        "hap1": "hap2",
        "hap2": "hap1",
        "hap3": "hap4",
        "hap4": "hap3"
    }
    def __init__(self, parent_sample, parent, grandparent_dict):
        self.sample = parent_sample
        self.parent = parent
        self.grandparents = grandparent_dict
        # self.other_parent = parent_dict[parent]
        # All known
        known_haplotypes = self.known_haplotypes
        inferred_haplotypes = self.inferred_haplotypes
        if grandparent_dict["dad"] and grandparent_dict["mom"]:
            if parent == "dad":
                self.haplotype1 = "A"
                self.haplotype2 = known_haplotypes["A"]
            else:
                self.haplotype1 = "C"
                self.haplotype2 = known_haplotypes["C"]
        elif grandparent_dict["dad"] and not grandparent_dict["mom"]:
            if parent == "dad":
                self.haplotype1 = "A"
                self.haplotype2 = inferred_haplotypes["A"]
            else:
                self.haplotype1 = "C"
                self.haplotype2 = inferred_haplotypes["C"]
        elif not grandparent_dict["dad"] and grandparent_dict["mom"]:
            if parent == "dad":
                self.haplotype2 = "B"
                self.haplotype1 = inferred_haplotypes["B"]
            else:
                self.haplotype2 = "D"
                self.haplotype1 = inferred_haplotypes["D"]
        else:
            if parent == "dad":
                self.haplotype1 = "hap1"
                self.haplotype2 = inferred_haplotypes["hap1"]
            else:
                self.haplotype1 = "hap3"
                self.haplotype2 = inferred_haplotypes["hap3"]

                

def assign_haplotype_snp(record, parent_info, grandparent_dict):
    if len([x for x in list(grandparent_dict.values()) if x is not None]) == 2:
        if record.call_for_sample[grandparent_dict['dad']].is_variant and record.call_for_sample[grandparent_dict['mom']].is_variant:
            grandparent ="/".join([grandparent_dict['dad'], grandparent_dict['mom']])
            haplotype ="/".join([parent_info.haplotype1, parent_info.haplotype2])
        elif record.call_for_sample[grandparent_dict['dad']].is_variant:
            grandparent = grandparent_dict['dad']
            haplotype = parent_info.haplotype1
        elif record.call_for_sample[grandparent_dict['mom']].is_variant:
            grandparent = grandparent_dict['mom']
            haplotype = parent_info.haplotype2
        else:
            grandparent = "unknown"
            haplotype = "unknown"
    elif len([x for x in list(grandparent_dict.values()) if x is not None]) == 1:
        if record.call_for_sample[grandparent_dict['dad']].is_variant:
            grandparent = grandparent_dict['dad']
            haplotype = parent_info.haplotype1
        elif record.call_for_sample[grandparent_dict['mom']].is_variant:
            grandparent = grandparent_dict['mom']
            haplotype = parent_info.haplotype2
        else:
            grandparent = "unknown"
            haplotype = "unknown"
    else:
        grandparent = "unknown"
        haplotype = "unknown"
    return grandparent, haplotype
        

def main():
    args = parser.parse_args()
    import vcfpy

    if hasattr(args, "func"):
        args.func(args)
    
    parents={
        "dad": args.dad_sample,
        "mom": args.mom_sample
    }

    pedigree=pd.read_csv(args.pedigree, sep='\t').astype(str)

    children_dict={
        "mom": list(pedigree.loc[pedigree["mom"] == parents["mom"]].to_dict()['child'].values()),
        "dad": list(pedigree.loc[pedigree["dad"] == parents["dad"]].to_dict()['child'].values())
    }

    # If subset_children parameter is used
    # Check all subset children are children of specified parents
    if args.subset_children:
        subset_children=set(args.subset_children.split(","))
        if subset_children.issubset(set(children_dict["mom"])) and subset_children.issubset(set(children_dict["dad"])):
            children_dict["mom"]=list(subset_children.intersection(set(children_dict["mom"])))
            children_dict["dad"]=list(subset_children.intersection(set(children_dict["dad"])))
        else:
            sys.exit(f"Subset children specified; {*subset_children,}, not children of mom/dad")

    children=list(set(children_dict["mom"]).intersection(set(children_dict["dad"])))
    # Sort children to match ordered dict from vcfpy
    children.sort()

    mom_parents={
        "dad": pedigree.loc[pedigree["child"] == str(parents["mom"])].to_dict('records')[0]['dad'],
        "mom": pedigree.loc[pedigree["child"] == str(parents["mom"])].to_dict('records')[0]['mom']
        }
    dad_parents={
        "dad": pedigree.loc[pedigree["child"] == str(parents["dad"])].to_dict('records')[0]['dad'],
        "mom": pedigree.loc[pedigree["child"] == str(parents["dad"])].to_dict('records')[0]['mom']
        }
    
    # Subset grandparents 
    if args.subset_gparents:
        grandparents=list(set(args.subset_gparents.split(",")))
        for k, v in mom_parents.items():
            if v not in grandparents:
                mom_parents[k] = None
        for k, v in dad_parents.items():
            if v not in grandparents:
                dad_parents[k] = None
    grandparents=list(mom_parents.values()) + list(dad_parents.values())
    

    
    # Depending on missing gparents
    ## Use a function for this
    
    ## gparent sample -> haplotype
    
    reader = vcfpy.Reader.from_path(args.cohort_calls)
    autosome = ["chr" + str(chrom) for chrom in  list(range(1,23))]
    dad = parent_info(parents["dad"], "dad", dad_parents)
    mom = parent_info(parents["mom"], "mom", mom_parents)
    header=[
        "#CHROM", "start", "end", 'REF', "ALT", "called_parent", "grandparent", "phase",
        #*[child + "_child" for child in children],
        "children_calls"#,
        #"transmission_vector_run", "run_id"
        ]
    print("\t".join(header))

    
    filtered = open(args.filtered, 'w', encoding="utf-8")
    filtered_header=header=[
        "#CHROM", "start", "end", 'REF', "ALT", "calls", "genotype", "depth", "filter"
        ]
    filtered.write("\t".join(filtered_header))

    mom_transmission_vector_run, dad_transmission_vector_run, run_id = (0, 0, 0)

    for record in reader:
        # Filter variants - what is returned here? #TODO
        filtered_variant = filter_variant(record, parents["mom"], parents["dad"], annotate_filtered_vars, filtered, children, grandparents, autosome)
        if filtered_variant is not None:
            record, called_parent, other_parent_sample = filtered_variant
        else:
            continue
        
        called_in=[]
        parents_called_in=[]
        gparents_called_in=[]
        #
        child_called_haplotypes=[]
        child_called_depth=[]
        haplotypes=[]
        for called in record.calls:
            if called.sample in children:
                haplotypes.append(called.is_variant)
                if called.is_variant is True:
                    called_in.append(called.sample)
                    child_called_haplotypes.append(called.data['GT'])
                    child_called_depth.append(called.data['DP'])
            elif called.sample in [parents["dad"]] + [parents["mom"]]:
                if called.is_variant is True:
                    parents_called_in.append(called.sample)
            elif called.sample in grandparents:
                if called.is_variant is True:
                    gparents_called_in.append(called.sample)
                else:
                    continue
            else:
                continue
        
        # Assign haplotypes
        if called_parent == parents["mom"]:
            grandparent, haplotype = assign_haplotype_snp(record, mom, mom_parents)
        elif called_parent == parents["dad"]:
            grandparent, haplotype = assign_haplotype_snp(record, dad, dad_parents)

        line = [
            record.CHROM, record.POS, record.POS, record.REF, record.ALT[0].value,
            called_parent, grandparent, haplotype
            ]
        line += [";".join(called_in)]


        print("\t".join(map(str, line)))
    
    filtered.close()


if __name__ == "__main__":
    main()