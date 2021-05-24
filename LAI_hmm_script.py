import random as rn
import numpy as np
import pandas as pd


def HammingDist(str1,str2):
    """
    This function computes the Hamming distance (edit distance) between 2 strings.
    """
    dist = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    return dist



def makeSNPseq(df,col1,col2):
    """
    This function takes a dataframe containing a single individual's SNP genotypes
    at various positions. It combines each row into one haplotype string (2nt) and
    returns these as a list.
    """
    snps = []
    for i in range(len(df)):
        snps.append(df.iloc[i,col1]+df.iloc[i,col2])

    return snps



def standardizeIndices(df1,df2,col_name):
    """
    This function cuts down df1 to only contain values in the column, col_name that exist in
    df2[,col_name]. It also checks that the two dataframes now have equal values in the two
    columns. In practice, this is used to cut down an emission df to perfectly match an
    individual's SNPs (in a genotype df) based on POS.
    """
    df2_pos = list(df2.loc[:,col_name])
    fin_df1 = df1[df1[col_name].isin(df2_pos)]
    print(df2_pos == list(fin_df1.loc[:,col_name]))
    return fin_df1



class HMMOptimalPathLAI:
    def __init__(self, PopA, PopB, emission_df, recombination_freq, SNPseq):
        """
        - States = ["AA", "AB", "BA", "BB"] --> states = {0:"AA",1:"AB",2:"BA",3:"BB"}
        - Emissions = we will calculate these as we go based on the emission_df (which contains
          the prob of each nt for popA and popB)
        - Alphabet = ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
        - Transitions = will build matrix using an input recombination frequency in separate function
        """
        #read in states, make a pop code dictionary, and construct an indexing state_dict {index:state}
        self.states = {0: "AA", 1: "AB", 2:"BA", 3:"BB"}
        self.n_states = len(self.states.keys())
        self.pop_codes = {'A':PopA, 'B':PopB}

        #read in emissions_df and recombination frequency and initialize the transition matrix
        self.emissions_df = emission_df
        self.recomb_freq = recombination_freq
        self.transitions = {0:[], 1:[], 2:[], 3:[]}

        #read in the SNPseq as a list
        self.sequence = SNPseq
        self.alphabet = ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]

        #pick a random starting state (0,1,2,3) and initialize the matrices for the Viterbi algorithm
        self.current_state = rn.randint(0, self.n_states-1)
        self.optimal_matrix = np.zeros((self.n_states, len(SNPseq)))
        self.backtrack = np.zeros((self.n_states, len(SNPseq)))
        self.final_path = []


    def get_transition_matrix(self):
        """
        This function creates the transition matrix based on the recombination rate given
        during initialization. By calculating the Hamming Distance (edit distance) between
        populations (AA,AB,BA,BB), we can see how many ancestries are changed in each
        transition and the use this and the recombination rate to assign a transition rate.
        """
        for i in range(self.n_states):
            state1 = self.states[i]
            for j in range(self.n_states):
                state2 = self.states[j]
                dist = HammingDist(state1,state2)
                if dist == 0:
                    self.transitions[i].append((1-self.recomb_freq)**2)
                elif dist == 1:
                    self.transitions[i].append(self.recomb_freq*(1-self.recomb_freq))
                elif dist == 2:
                    self.transitions[i].append(self.recomb_freq**2)


    def get_emissions_matrix(self,ind):
        """
        This function takes a SNP position index (works bc we have matched indexing for genotype
        and emissions), and uses the emission_df to output a 4x16 emission probability matrix
        (rows = 4 ancestry pop combos, cols = 16 genotype combos). The genotype combos are ordered
        lexicographically: ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC",
        "TG","TT"].
        """
        #go through the columns with pop_nt freqs and make the emissions matrix
        emissions_mat = np.zeros((4,16),dtype=float)
        for i in range(self.n_states):
            pop1 = self.pop_codes[self.states[i][0]]
            pop2 = self.pop_codes[self.states[i][1]]
            for j in range(len(self.alphabet)):
                nt1 = self.alphabet[j][0]
                nt2 = self.alphabet[j][1]
                prob1 = self.emissions_df.loc[ind,pop1+"_"+nt1]
                prob2 = self.emissions_df.loc[ind,pop2+"_"+nt2]
                emissions_mat[i,j] = prob1*prob2

        return emissions_mat


    def get_optimal_path(self):
        """
        This function performs the Viterbi algorithm using the input sequence of SNPs.
        For each SNP, a position (index) specific emissions matrix is calculated and then used
        to calculate the highest possible probability of being in each state for each SNP.
        A backtrack matrix is stored which will facilitate reconstruction of the most probable
        path through states, and thus the haplotypes of a chromosome.
        """
        #go through every SNP in the input sequence (i=index, c=SNP string)
        for i, c in enumerate(self.sequence):
            #pull out the lexicographic indexing of the seen 2nt SNP (idx)
            idx = self.alphabet.index(c)
            #create the emissions matrix for this SNP with the enumerated index
            snp_emissions = self.get_emissions_matrix(i)

            #go through all states and calculate max probability based on emission and transition probs
            for j in range(self.n_states):
                #get the probability of emitting c (SNP with index idx) when in state j
                emit_prob = snp_emissions[j][idx]

                #if this is the first SNP, then just set all positions in the matrix to the emissions prob
                #(this assumes equal prob of start in any state, which is fine here)
                if i == 0:
                    self.optimal_matrix[j, i] =  emit_prob
                    self.backtrack[j, i] = -1

                #otherwise go through all possible previous states and keep track of transition
                #and emission probabilities from these states to the current one
                else:
                    probs = []
                    for previous in range(self.n_states):
                        transition_prob = self.transitions[previous][j]
                        probs.append(emit_prob * transition_prob * self.optimal_matrix[previous, i-1])
                    #select the max probability and record backtrack
                    max_prob = max(probs)
                    back = probs.index(max_prob)
                    self.optimal_matrix[j, i] = max_prob
                    self.backtrack[j, i] = back


    def reconstruct_path(self):
        """
        This function uses the backtrack list stored from the Viterbi algorithm to reconstruct
        the optimal path of states through the input sequence (SNPs).
        """
        #find the ending state (state with highest prob for final SNP) and add it to the path
        end_values = [self.optimal_matrix[v, len(self.sequence) - 1] for v in range(self.n_states)]
        max_prob = max(end_values)
        end_state = end_values.index(max_prob)
        path = [end_state]
        previous_state = end_state

        #now walk backwards through the backtrack list until reaching the start
        for i in range(len(self.sequence)-1, -1, -1):
            #print(previous_state, i)
            path.append(int(self.backtrack[previous_state, i]))
            previous_state = int(self.backtrack[previous_state, i])

        #reverse path and output with state names
        path.reverse()
        self.final_path = [self.states[p] for p in path if p != -1]


    def output_path(self,chr_num):
        """
        This function translates the final path to a dataframe of the desired output
        format: cols = CHR, POS, POP1, POP2. It then returns this dataframe so it can
        be viewed and saved externally.
        """
        chrom = [chr_num] * len(self.final_path)
        pos = list(self.emissions_df.loc[:,"POS"])

        #reconstruct the two populations from the final path using the pop_codes dict
        hap1 = []
        hap2 = []
        for gt in self.final_path:
            pop1 = gt[0]
            pop2 = gt[1]
            hap1.append(self.pop_codes[pop1])
            hap2.append(self.pop_codes[pop2])

        df_dict = {'CHR':chrom, 'POS':pos, 'POP1':hap1, 'POP2':hap2}
        results_df = pd.DataFrame(df_dict)
        return(results_df)
