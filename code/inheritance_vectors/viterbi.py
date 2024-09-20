#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
import numpy as np
import itertools


#Define usage for command and subcommands
parser = argparse.ArgumentParser(
    description="Script to run viterbi algo on inheritance sites",
    usage="inht_vectors_viterbi.py <command> [<args>]")
parser.add_argument(
    "-i", "--input", required=True, help="States sites file",
    type=str)
parser.add_argument(
    "-c", "--children", required=True, help="Children, comma separated list",
    type=str)
parser.add_argument(
    "-m", "--male-children", required=True, help="Male children, comma separated list",
    type=str)
parser.add_argument(
    "-l", "--parents-list", required=True, help="Comma seperated list of parent sample ids [dad],[mom] (eg NA12877,NA12878)",
    type=str)
parser.add_argument(
    "-t", "--transmission-matrix", required=True, help="Out t mat",
    type=str)
parser.add_argument(
    "-e", "--emission-matrix", required=True, help="Out e mat",
    type=str)
parser.add_argument(
    "-d", "--test-outdir", required=True, help="out",
    type=str)
parser.add_argument(
    "-p", "--punishment", required=True, help="Comma seperated list of punishment values for snp error, and state change respectively",
    type=str)
parser.add_argument(
    "-o", "--output", required=True, help="out",
    type=str)
parser.add_argument(
    "-f", "--file-prefix", required=True, help="out prefix",
    type=str)



# Convert states to be in context of 1st haplotype, make common state #TODO
def convert_state(state, haplotype, all_states):
    if haplotype in ["B", "D", "hap2", "hap4"]:
        common_state = set(all_states) ^ set(state)
    else:
        common_state = state
    return common_state

def state_difference(state1, state2, all_states):
    names1 = state1.split(',')
    names2 = state2.split(',')
    
    # Check if states are complements
    if set(names1).isdisjoint(names2):
        return 1 if state1 != state2 else 0
    
    
    num_differences = len(set(names1) ^ set(names2))
    num_differences_complement = len(set(set(names1) ^ set(all_states)) ^ set(names2))
    if num_differences_complement < num_differences:
        num_differences = num_differences_complement
    
    return num_differences

# Function to create transition likelihood matrix with adjustable parameters
def create_transition_matrix(states, stay_same_likelihood, change_punishment):
    """_summary_

    Args:
        states (_type_): _description_
        stay_same_likelihood (_type_): _description_
        error_probability (_type_): _description_

    Returns:
        _type_: _description_
    """
    #states = convert_state(states, haplotype, all_states)
    num_states = len(states)
    transition_matrix = np.zeros([num_states, num_states])
    
    for i in range(num_states):
        for j in range(num_states):

            dissimilarity = len(set(states[i]) ^ set(states[j]))
            
            if dissimilarity == 0:
                likelihood = stay_same_likelihood
            else:
                likelihood = (stay_same_likelihood - change_punishment) ** dissimilarity
            
            transition_matrix[i][j] = likelihood
    
    return transition_matrix

# Function to calculate the likelihood of observed states given hidden states
def create_emission_matrix(hidden_states, observed_states, error_punishment):
    """_summary_

    Args:
        hidden_states (_list_): _description_
        observed_states (_list_): _description_
        error_probability (_float_): _description_

    Returns:
        _type_: _description_
    """
    num_hidden_states = len(hidden_states)
    num_observed_states = len(observed_states)
    
    # Initialize emission matrix
    emission_matrix = np.zeros([num_hidden_states, num_observed_states])
    
    # Define error probability (presence when absent)
    
    for i, hidden_state in enumerate(hidden_states):
        for j, observed_state in enumerate(observed_states):
            # Calculate the number of differences between observed and hidden states
            observed_set = set(observed_state)
            hidden_set = set(hidden_state)
            
            # Calculate likelihood based on the number of differences
            if observed_set == hidden_set:
                likelihood = 0.99
            else:
                likelihood = (0.99 - error_punishment) ** 1
            
            emission_matrix[i][j] = likelihood
    
    return emission_matrix
    

def viterbi(y, A, B, pi, state_index, haplotype_list, all_states, chromosome, male_children, args):
    """
        viterbi algorithm
        :param y: observation sequence
        :param A: the transition matrix
        :param B: the emission matrix
        :param pi: the initial probability distribution
        :param state_index: Dictionary to translate between the state and its index row/columns in the e and t matrices
        https://medium.com/@zhe.feng0018/coding-viterbi-algorithm-for-hmm-from-scratch-ca59c9203964
        
        Returns:
        :x_seq_opt: optimal sequence
    """
    
    # Alter t mat and e mat for chrX
    ## 
    parents_dict={
        "dad": args.parents_list.split(",")[0],
        "mom": args.parents_list.split(",")[1]
    }
    #
    
    N = B.shape[0]
    x_seq = np.zeros([N, 0])

    V = B[:, state_index[y[0]]] * pi

    # forward to compute the optimal value function V
    for position, y_ in enumerate(y[1:]):
        
        
        #set up state indexing
        observed_state = y_
        
        haplotype = haplotype_list[position + 1]
        converted_state = list(convert_state(observed_state, haplotype, all_states))
        converted_state.sort()
        if not converted_state:
            converted_state = tuple()
        
        # Change from state to it's index
        obs_y_ = state_index[tuple(y_)]
        y_ = state_index[tuple(converted_state)]
        
        _V = np.tile(B[:, y_], reps=[N, 1]).T * A.T * np.tile(V, reps=[N, 1])

        x_ind = np.argmax(_V, axis=1)
        x_seq = np.hstack([x_seq, np.c_[x_ind]])
        V = _V[np.arange(N), x_ind]
    x_T = np.argmax(V)
    

    # backward to fetch optimal sequence
    x_seq_opt, i = np.zeros(x_seq.shape[1]+1), x_seq.shape[1]
    prev_ind = x_T
    while i >= 0:
        x_seq_opt[i] = prev_ind
        i -= 1
        prev_ind = x_seq[int(prev_ind), i]
    
    return x_seq_opt
    
def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def search_dict(dict, search_index):
    return [state for state, index in dict.items() if index == search_index]
    
def main():
    args = parser.parse_args()
    test_outdir = args.test_outdir
    prefix = args.file_prefix
    
    all_individuals=args.children.split(";")
    male_children=args.male_children.split(";")
    
    possible_states=list(powerset(all_individuals))
    state_index = {k: v for v, k in enumerate(possible_states)}
    
    snp_error_punishment = float(args.punishment.split(",")[0])
    change_punishment = float(args.punishment.split(",")[1])
    
    # Create the emission matrix
    emission_matrix = create_emission_matrix(possible_states, possible_states, snp_error_punishment)

    # Create the transition likelihood matrix with custom parameters
    transition_matrix = create_transition_matrix(possible_states, 0.95, change_punishment)

    #print("Transition Likelihood Matrix:")
    #print(transition_matrix)
    transition_matrix_df = pd.DataFrame(transition_matrix.astype(float))
    transition_matrix_df.to_csv(args.transmission_matrix, sep ="\t", header = possible_states, index = possible_states)
    
    #print("Emission matrix:")
    #print(emission_matrix)
    emission_matrix_df = pd.DataFrame(emission_matrix.astype(float))
    emission_matrix_df.to_csv(args.emission_matrix, sep ="\t", header = possible_states, index = possible_states)
    
    sites = pd.read_csv(args.input, sep="\t")
    parents = list(sites["called_parent"].unique())
    
    chromosomes = ["chr" + str(chrom) for chrom in  list(range(1,23))] + ["chrX"]
    chromosomes = list(set.intersection(set(chromosomes), set(list(sites["CHROM"].unique()))))

    with open(args.output, "w") as f:
        for chrom in chromosomes:
        # Run viterbi per parent, per chromosome
            for parent in parents:
                
                filtered_df = sites.query('CHROM == @chrom & called_parent == @parent & phase in ["A", "B", "C", "D", "hap1", "hap2", "hap3", "hap4"]').dropna()
                
                observed_sequence = pd.Series(sites.query('CHROM == @chrom & called_parent == @parent & phase in ["A", "B", "C", "D", "hap1", "hap2", "hap3", "hap4"]')[['children_calls', 'phase']].dropna()['children_calls']).to_list()
                
                haplotype_list = pd.Series(sites.query('CHROM == @chrom & called_parent == @parent & phase in ["A", "B", "C", "D", "hap1", "hap2", "hap3", "hap4"]')[['children_calls', 'phase']].dropna()['phase']).to_list()
                observed_sequence = [tuple(x.split(";")) for x in observed_sequence]
                optimal_sequence = viterbi(observed_sequence, transition_matrix, emission_matrix, (0.5 ** 8), state_index=state_index, haplotype_list=haplotype_list, all_states=all_individuals, chromosome=chrom, male_children=male_children, args=args)
                
                filtered_df['x_seq_opt']=list(optimal_sequence)
                
                
                states_list=[]
                for i in optimal_sequence:
                    states_list.append(search_dict(state_index, i))
                
                filtered_df['x_seq_opt_state']=states_list
                
                filtered_df.to_csv(f"{test_outdir}/{prefix}.{chrom}_{parent}.viterbi_df.txt", sep = "\t")
                # Convert index back to inheritance state
                f.write(f"chrom: {chrom}, parent: {parent}\n")
                f.write(f"optimal sequence:\n")
                np.savetxt(f, optimal_sequence, delimiter = "\t", fmt='%f')
            

if __name__ == "__main__":
    main()