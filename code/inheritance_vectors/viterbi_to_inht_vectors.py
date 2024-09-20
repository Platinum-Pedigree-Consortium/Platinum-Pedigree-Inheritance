#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import numpy as np
import itertools
import pickle
from ast import literal_eval

parser = argparse.ArgumentParser(
    description="Script to convert viterbi to inht vectors",
    usage="viterbi_to_inht_vectors.py <command> [<args>]")
parser.add_argument(
    "-d", "--input-dad", required=True, help="States sites file",
    type=str)
parser.add_argument(
    "-m", "--input-mom", required=True, help="States sites file",
    type=str)
parser.add_argument(
    "-c", "--children", required=True, help="Semi colon seperated list of children sample ids",
    type=str)
parser.add_argument(
    "-mc", "--male-children", required=True, help="Semi colon seperated list of male children sample ids",
    type=str)
parser.add_argument(
    "-p", "--parents", required=True, help="Comma seperated list of parent sample ids [dad],[mom] (eg NA12877,NA12878)",
    type=str)
parser.add_argument(
    "-o", "--output-inht", required=True, help="out recomb bed file in inht vectors format",
    type=str)
parser.add_argument(
    "-g", "--output-gaps", required=True, help="out gaps, prefix, per child",
    type=str)
parser.add_argument(
    "-a", "--output-analysis", required=True, help="out snps joined to recombs for simplifying analysis",
    type=str)
parser.add_argument(
    "-s", "--output-summary", required=True, help="Summary file",
    type=str)


def filter_values(start, end, values):
    return [value for value in values if start <= value <= end]

def get_recombinations(file_path, parents_dict, sample_id, autosome, children):
    #
    parent=parents_dict[sample_id]
    
    viterbi_sites=pd.read_csv(file_path, sep="\t")
    
    # check unique chrom
    chrom = viterbi_sites['CHROM'].unique()
    if len(list(chrom)) > 1:
        raise Exception("More than one chromosome found in dataframe")
    
    viterbi_sites['previous_same'] = viterbi_sites.x_seq_opt_state.eq(viterbi_sites.x_seq_opt_state.shift())
    viterbi_sites['next_same'] = viterbi_sites.x_seq_opt_state.eq(viterbi_sites.x_seq_opt_state.shift(-1))
    viterbi_sites['x_seq_opt_state'] = viterbi_sites['x_seq_opt_state'].str.strip("[]").apply(literal_eval)
    
    recombs = viterbi_sites.loc[~viterbi_sites['previous_same'] == viterbi_sites['next_same']]
    recombs['first_or_last'] = False
    recombs.loc[recombs.index.isin([recombs.index[0], recombs.index[-1]]), 'first_or_last'] = True
    
    
    recombs['recombination'] = ~(recombs['previous_same'] & recombs['next_same']) & ~recombs['first_or_last']
    recombs['previous_state'] = recombs['x_seq_opt_state'].shift(1)
    recombs['difference_state'] = recombs.apply(lambda row: list(set(row['x_seq_opt_state']) ^ set(row['previous_state'])) if row['previous_state'] is not None else list(), axis=1)
    recombs['next_start'] = recombs['start'].shift(-1)
    recombs.drop(["REF", "ALT"], axis=1, inplace=True)
    
    recombs['rle'] = recombs['x_seq_opt_state'].ne(recombs['x_seq_opt_state'].shift()).cumsum()
    haplotypes = recombs.groupby('rle').first()
    haplotypes['end'] = haplotypes['next_start']
    haplotypes['x_seq_opt_state'] = haplotypes['x_seq_opt_state'].apply(list)
    
    haplotypes.reset_index(inplace=True)
    
    haplotypes_pivot = haplotypes[['rle', 'x_seq_opt_state']].explode('x_seq_opt_state').pivot_table(index="rle", columns =["x_seq_opt_state"], aggfunc=lambda x: 1, fill_value=0)
    
    
    # Check if any children are missing (ie no recombinations)
    # This should only happen if child is B or D for the entire chromosome (ie not listed)
    if len(set(children) ^ set(haplotypes_pivot.columns)) > 0:
        missing_children = set(children) ^ set(haplotypes_pivot.columns)
        for missing_child in missing_children:
            haplotypes_pivot[missing_child] = 0
    
    haplotypes_pivot.sort_index(axis=1, inplace=True)
    
    haplotypes = haplotypes[['CHROM', 'start', 'end', 'rle']].set_index('rle').join(haplotypes_pivot)
    
    
    if chrom in autosome:
        if parent == "dad":
            haplotypes.fillna(0, inplace=True)
            haplotypes.replace({1: "A", 0: "B"}, inplace=True)
        else:
            haplotypes.fillna(0, inplace=True)
            haplotypes.replace({1: "C", 0: "D"}, inplace=True)
    if chrom == "chrX":
        if parent == "dad":
            haplotypes.fillna(0, inplace=True)
            haplotypes.replace({1: "B", 0: "B"}, inplace=True)
        else:
            haplotypes.fillna(0, inplace=True)
            haplotypes.replace({1: "C", 0: "D"}, inplace=True)
    haplotypes.reset_index(inplace=True)
    haplotypes=haplotypes.set_index(['CHROM', 'start', 'end']).add_suffix(f"_{parent}").reset_index()
    
    # Add supporting SNPs, we will use this at the end
    supporting_snps = haplotypes.apply(lambda row: filter_values(row['start'], row['end'], list(viterbi_sites['start'])), axis=1)
    haplotypes['supporting_snps'] = supporting_snps
    return chrom, haplotypes

    
def output_haplotypes_summary(combined, children, chrom, output_summary_file):
    print("test this function")
    header = "\t".join(["chrom", "child", "haplotypes", "recombinations", "gap_size", "A_inh", "B_inh", "C_inh", "D_inh", "dad_recombs", "mom_recombs"])
    
    with open(output_summary_file, "w") as outfile:
        outfile.write(header + "\n")
        for child in children:
            child_df = combined[['start', 'end', child]]
            child_df['rle'] = combined[child].ne(combined[child].shift()).cumsum()
            
            child_df_grouped = child_df.groupby('rle')
            last_values = child_df_grouped['end'].last()
            first_values_next = child_df_grouped['start'].first().shift(-1)
            differences = first_values_next - last_values
            gap_size = int(differences.sum())
            num_haplotypes = max(child_df['rle'])
            num_recombs = max(child_df['rle']) - 1
            
            # Haplotypes bases inherited
            child_df['length'] = child_df['end'] - child_df['start']
            A_length = int(child_df[child_df[child].str.contains("A")]["length"].sum())
            B_length = int(child_df[child_df[child].str.contains("B")]["length"].sum())
            C_length = int(child_df[child_df[child].str.contains("C")]["length"].sum())
            D_length = int(child_df[child_df[child].str.contains("D")]["length"].sum())
            
            # Dad haplotypes
            dad_haps = child_df[['start', 'end', child]]
            dad_haps[child] = dad_haps[child].str[0]
            dad_haps['rle'] = dad_haps[child].ne(dad_haps[child].shift()).cumsum()
            num_dad_recombs = max(dad_haps['rle']) - 1
            
            # Mom haplotypes
            mom_haps = child_df[['start', 'end', child]]
            mom_haps[child] = mom_haps[child].str[1]
            mom_haps['rle'] = mom_haps[child].ne(mom_haps[child].shift()).cumsum()
            num_mom_recombs = max(mom_haps['rle']) - 1

            line = "\t".join(str(x) for x in [chrom[0], child, num_haplotypes, num_recombs, gap_size, A_length, B_length, C_length, D_length, num_dad_recombs, num_mom_recombs])
            outfile.write(line + "\n")

def dif(a, b):
    return [i for i in range(len(a)) if a[i] != b[i]]

def get_hap_difference(haplotype_1, haplotype_2):
    if pd.isna(haplotype_1) or pd.isna(haplotype_2):
        return ''
    else:
        difference = list(set(list(haplotype_1)) ^ set(list(haplotype_2)))
        if len(difference) == 2:
            return difference[1]
        elif len(difference) == 4:
            return "".join([difference[1], difference[3]])
        else:
            return None
        
    
def merge_shifted(row, df):
    merged = ''
    for col_idx, col in enumerate(df.columns):
        if col.endswith('recombination'):
            merged += str(row[col])
            if col_idx < len(df.columns) - 1 and len(str(row[col])) > 0:
                merged += ";"
    return merged

def output_gaps_analysis(chrom, combined, output_prefix, children, parents_dict):
    
    for col in combined[children]:
        combined[f"{col}_shift"] = combined[col].shift()
        combined[f"{col}_diff"] = combined.apply(lambda x: get_hap_difference(x[col], x[f"{col}_shift"]), axis=1)
        
        conditions = [(combined[f"{col}_diff"] == "A") | (combined[f"{col}_diff"] == "B"), (combined[f"{col}_diff"] == "C") | (combined[f"{col}_diff"] == "D"), (combined[f"{col}_diff"] == "AC") | (combined[f"{col}_diff"] == "AD") | (combined[f"{col}_diff"] == "BC") | (combined[f"{col}_diff"] == "BD") , (combined[f"{col}_diff"].isnull())]
        choices = [f"{col}_NA12877", f"{col}_NA12878", f"{col}_both", ""]
        combined[f"{col}_recombination"] = np.select(conditions, choices, default="NA")
    
    combined["recombination"] = combined.apply(merge_shifted, args=(combined,), axis=1)
    combined["rle_recomb"] = combined["recombination"].ne(combined["recombination"].shift(-1)).cumsum()
    combined["size"] = (combined.groupby('rle_recomb').transform('size'))
    combined["dbl_recomb"] = np.where(combined["size"] > 1, "dbl", "single")
    
    
    for child in children:
        
        # Keep one child column
        child_df = combined[['CHROM', 'start', 'end', child, "recombination", "rle_recomb", "dbl_recomb"]]
        # RLE to identify differences in inheritance
        child_df['rle'] = child_df[child].ne(child_df[child].shift()).cumsum()
        child_df_grouped = child_df.groupby('rle')

        test = child_df_grouped[child].apply(lambda x: pd.DataFrame(list(x.str.split('').dropna()))).iloc[:, 1:-1].groupby('rle').first()
        test['rle_dad'] = test[1].ne(test[1].shift()).cumsum()
        test['rle_mom'] = test[2].ne(test[2].shift()).cumsum()
        
        test['rle_dad_shift'] = test['rle_dad'].shift(-1)
        test['rle_mom_shift'] = test['rle_mom'].shift(-1)
        test['recomb_parent_dad'] = np.where(test['rle_dad']!=test['rle_dad_shift'], 'NA12877', "")
        test['recomb_parent_mom'] = np.where(test['rle_mom']!=test['rle_mom_shift'], 'NA12878', "")
        test['recomb_parent'] = test['recomb_parent_dad'] + test['recomb_parent_mom']
        gap_start = list(child_df_grouped['end'].last())
        gap_end = list(child_df_grouped['start'].first())[1:]
        gap_end.append(np.nan)
        haps = list(child_df_grouped[child].last())
        dbls = child_df_grouped["dbl_recomb"].last()
        
        
        test = {
            "CHROM": [chrom] * len(gap_start),
            "start": gap_start,
            "end": gap_end,
            "haplotype": haps,
            "child": child,
            "diff": test['recomb_parent'],
            "dbl": dbls
        }
        test_df = pd.DataFrame(test)

        test_df.reset_index(inplace=True)
        test_df.drop(test_df.tail(1).index,inplace=True)
        test_df[['start', 'end']] = test_df[['start', 'end']].astype(int)
        test_df['CHROM'] = test_df['CHROM'].str[0]
        test_df = test_df.drop('rle', axis=1)
        # Keep only gaps with differences
        # Label difference (mom or dad id)
        # Label gap change (a->b, c->d etc)
        
        outfile = f"{output_prefix}.{child}.tsv"
        test_df.to_csv(outfile, sep='\t', index=False)
    
    
    
    
    


def main():
    # Load viterbi output
    args = parser.parse_args()
    dad_id = args.parents.split(",")[0]
    mom_id = args.parents.split(",")[1]
    
    children=args.children.split(";")
    male_children=args.male_children.split(";")
    
    parents_dict={
        dad_id: "dad",
        mom_id: "mom"
    }
    
    autosome=["chr" + str(chrom) for chrom in  list(range(1,23))]
    
    
    chrom, mom_recombs=get_recombinations(args.input_mom, parents_dict, mom_id, autosome, children)
    if chrom in autosome:
        dad_chrom, dad_recombs=get_recombinations(args.input_dad, parents_dict, dad_id, autosome, children)
    elif chrom =="chrX":
        dad_recombs = pd.DataFrame(columns=['CHROM', 'start', 'end'] + [child + "_dad" for child in children])
    
    # print(mom_recombs)
    # print(dad_recombs)
    
    if chrom in autosome:
        combined=pd.concat([dad_recombs, mom_recombs], ignore_index=True, sort=True).sort_values(['start'])
        combined=combined.fillna(method='ffill')
        combined=combined.fillna(method='bfill')
        # combined[children]=combined[children].fillna(method='ffill')
        # combined[children]=combined[children].fillna(method='bfill')
        # print("Combined 1:")
        # print(combined)
        for child in children:
            child_columns = combined[combined.columns[pd.Series(combined.columns).str.startswith(child)]]
            #print(child_columns)
            combined[child] = child_columns.iloc[:, 0] + child_columns.iloc[:, 1]
    elif chrom == "chrX":
        combined=pd.concat([dad_recombs, mom_recombs], ignore_index=True, sort=True).sort_values(['start'])
        for child in children:
            child_columns = combined[combined.columns[pd.Series(combined.columns).str.startswith(child)]]
            #print(child_columns)
            if child in male_children:
                combined[child] = child_columns.iloc[:, 1] + child_columns.iloc[:, 1]
            else:
                combined[child] = "B" + child_columns.iloc[:, 1]
                
    
    # Add parent columns
    if chrom in autosome:
        combined[list(parents_dict.keys())[list(parents_dict.values()).index("dad")]] = "AB"
        combined[list(parents_dict.keys())[list(parents_dict.values()).index("mom")]] = "CD"
    elif chrom == "chrX":
        
        combined[list(parents_dict.keys())[list(parents_dict.values()).index("dad")]] = "BB"
        combined[list(parents_dict.keys())[list(parents_dict.values()).index("mom")]] = "CD"
        
    
    # Fix duplicate columns at start of each chrom
    # Group by vectors - rle
    combined['paste_children'] = combined[children].agg("".join, axis=1)
    combined['rle'] = combined['paste_children'].ne(combined['paste_children'].shift()).cumsum()
    first_dict = dict.fromkeys(children + ["CHROM", "start", "supporting_snps"] + list(parents_dict.keys()), 'first')
    last_dict = dict.fromkeys(['end'], 'last')
    agg_dict = {**first_dict, **last_dict}
    combined=combined.groupby('rle').agg(agg_dict).reset_index(drop=True)
    
    
    # Fix overlapping vectors, using original SNP support
    combined['next_start'] = combined['start'].shift(-1)
    for index, row in combined.iterrows():
        if row['next_start'] and row['end'] > row['next_start']:
            max_value = max(i for i in row['supporting_snps'] if i <= row['next_start'])
            combined.at[index, 'end'] = max_value 
    #print(combined)
    
    # End can be less than start for same vector, how is this happening
    
    
    
    out_df = combined[['CHROM', 'start', 'end'] + children + list(parents_dict.keys())]
    out_df['end'] = out_df['end'].astype(int)
    
    out_df.to_csv(args.output_inht, index=False)
    
    output_haplotypes_summary(combined, children, chrom, args.output_summary)
    
    output_gaps_analysis(chrom, combined, args.output_gaps, children, parents_dict)

if __name__ == "__main__":
    main()