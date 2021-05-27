#!/opt/conda/bin/python
from collections import OrderedDict
import pandas as pd
import os
import numpy as np

chr_maf_file = "/home/lbruce/teams/CSE284_SP21_A00/team3/chromosome_21_files/chr21_genotypes_afr_eur_allelefreqs.csv"
#chr_maf_file = "/home/lbruce/teams/CSE284_SP21_A00/team3/chromosome_14_files/chr14_genotypes_afr_eur_allelefreqs.csv"
#chr_maf_file = "/home/lbruce/teams/CSE284_SP21_A00/team3/chromosome_14_files/chr14_genotypes_afr_eur_allelefreqs_HEAD.csv"
genotype_maf_df = pd.read_csv(chr_maf_file)
genotype_maf_df = genotype_maf_df[["CHR", "POS", "SNP","A2", "A1", "MAF_AFR", "MAF_EUR"]]
#genotype_maf_df.rename(columns={"A2": "REF", "A1": "ALT"}, inplace=True)

# Keep rows only if the length of A1 or A2 is == 1
genotype_maf_df = genotype_maf_df[genotype_maf_df['A1'].apply(lambda x: len(x) == 1)]
genotype_maf_df = genotype_maf_df[genotype_maf_df['A2'].apply(lambda x: len(x) == 1)]

def calculate_bp_prob(row, genotype_maf_df):

    # Start with all bps having a 0 probability
    # and only update if bp is in A1 or A2
    afr_prob_dict = OrderedDict([('A',0.000001), ('C',0.000001), ('G',0.000001), ('T',0.000001)])
    eur_prob_dict = OrderedDict([('A',0.000001), ('C',0.000001), ('G',0.000001), ('T',0.000001)])

    # A1 = Alternate allele, P = MAF
    # A2 = Reference allele, P = 1-MAF

    # Set reference alelle probability
    afr_prob_dict[row['A2']] = 1 - float(row['MAF_AFR'])
    eur_prob_dict[row['A2']] = 1 - float(row['MAF_EUR'])

    # Set alternate alelle probability
    afr_prob_dict[row['A1']] = float(row['MAF_AFR'])
    eur_prob_dict[row['A1']] = float(row['MAF_EUR'])
    

    # If any values equal 0, replace with 0.00000001
    for key,value in afr_prob_dict.items():
        if value < 0.00001:
            afr_prob_dict[key] = 0.000001
    for key,value in eur_prob_dict.items():
        if value < 0.00001:
            eur_prob_dict[key] = 0.000001
            
    #print(afr_prob_dict)
    #print(eur_prob_dict)
    
    return(afr_prob_dict, eur_prob_dict)

    


# Loop through all SNPs and calculate the Values
i = 0
AFR_A =[]
AFR_C =[]
AFR_G =[]
AFR_T =[]
EUR_A =[]
EUR_C =[]
EUR_G =[]
EUR_T =[]
for index, row in genotype_maf_df.iterrows():
    #print(row)
        
    (afr_prob_dict, eur_prob_dict) = calculate_bp_prob(row, genotype_maf_df)
    AFR_A.append(afr_prob_dict['A'])
    AFR_C.append(afr_prob_dict['C'])
    AFR_G.append(afr_prob_dict['G'])
    AFR_T.append(afr_prob_dict['T'])

    EUR_A.append(eur_prob_dict['A'])
    EUR_C.append(eur_prob_dict['C'])
    EUR_G.append(eur_prob_dict['G'])
    EUR_T.append(eur_prob_dict['T'])

    #break

#print(AFR_A)
#print(len(AFR_A))

# Store the BP lists in the dataframe and save to file
genotype_maf_df['AFR_A'] = AFR_A
genotype_maf_df['AFR_C'] = AFR_C
genotype_maf_df['AFR_G'] = AFR_G
genotype_maf_df['AFR_T'] = AFR_T

genotype_maf_df['EUR_A'] = EUR_A
genotype_maf_df['EUR_C'] = EUR_C
genotype_maf_df['EUR_G'] = EUR_G
genotype_maf_df['EUR_T'] = EUR_T

output_file = "/home/lbruce/teams/CSE284_SP21_A00/team3/chromosome_21_files/chr21_genotypes_afr_eur_allelefreqs.bybp.csv"
#output_file = "/home/lbruce/teams/CSE284_SP21_A00/team3/chromosome_14_files/chr14_genotypes_afr_eur_allelefreqs.bybp.csv"

genotype_maf_df.to_csv(output_file, index=False)

