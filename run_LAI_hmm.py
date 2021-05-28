#script to take fps for input (emissions, genotype) /output filepaths and run all the HMM functions
#also outputs the probability matrix

import random as rn
import numpy as np
import pandas as pd
import sys
import math

import datetime
from LAI_hmm_scriptFINAL import HMMOptimalPathLAI, HammingDist, makeSNPseq, standardizeIndices
from LAI_hmm_script_progress import HMMOptimalPathLAI_progress


def HMM_implementation_wrapper(popA, popB, genotype_fp, emission_df, recomb, output_fp, chrm):
    #read in sample genotype to test -- assumes current dir is the team3 dir
    genotype_df = pd.read_csv(genotype_fp, sep='\t',header = 0)

    #cut down both dfs to only be positions shared in both
    fin_emission_df, fin_genotype_df = standardizeIndices(emission_df, genotype_df, "POS")

    #create SNPseq
    a1_idx = list(fin_genotype_df.columns).index('A1')
    a2_idx = list(fin_genotype_df.columns).index('A2')
    snps = makeSNPseq(fin_genotype_df,a1_idx,a2_idx)

    #create the hmm object
    hmm = HMMOptimalPathLAI_progress(popA, popB, fin_emission_df, recomb, snps)
    hmm.get_transition_matrix()

    #perform the Viterbi algorithm and reconstruct the most probable path through states
    hmm.get_optimal_path()
    prob_mat = hmm.optimal_matrix
    #print("Final probability matrix")
    #print(prob_mat)

    hmm.reconstruct_path()

    #convert the path to the desired output format (pass the string to use for CHR col)
    path_df = hmm.output_path(chrm)
    path_df.to_csv(output_fp,sep='\t',header=True,index=False)



if __name__ == '__main__':
    input = sys.stdin.read()
    data = list(map(str, input.split()))
    emission_fp = data[0]
    genotype_fp = data[1]
    output_fp = data[2]

    #define names of the two test populations
    test_PopA = "0"
    test_PopB = "1"

    #read in emissions data (simulated df)
    emission_df1 = pd.read_csv(emission_fp, sep='\t',header=0)

    #renaming the columns (only necessary for simulated data)
    colnames = ["POS","0_A","0_C","0_G","0_T","1_A","1_C","1_G","1_T"]
    emission_df1.columns = colnames
    #emission_df1.head()

    #read in sample genotype to test
    #assumes current dir is the team3 dir
    genotype_df1 = pd.read_csv(genotype_fp, sep='\t',header = 0)
    #genotype_df1.head()

    HMM_implementation_wrapper(test_PopA, test_PopB, genotype_fp1, emission_df1, 0.1, output_fp, "21")
