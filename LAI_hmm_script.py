import random as rn
import numpy as np


class HMMOptimalPathLAI:
    def __init__(self, states, emissions, emission_df, SNPseq, alphabet, PopA, PopB):
        """
        - States = ["AA", "AB", "BA", "BB"] --> state_dict = {0:"AA",1:"AB",2:"BA",3:"BB"}
        - Emissions = we will calculate these as we go based on the emission_df (which contains
          the prob of each nt for popA and popB)
        - Alphabet = ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"]
        -
        """
        #read in states, make a pop code dictionary, and construct an indexing state_dict {index:state}
        self.n_states = len(states)
        self.pop_codes = {'A':PopA, 'B':PopB}
        self.state_dict = {}
        for s in range(self.n_states):
            self.state_dict[s] = states[s]

        #read in emission options and construct an indexing emission_dict {index:symbol}
        ### NOTE: may not need the emissions_dict, considering we have the alphabet list ########
        self.n_emissions = len(emissions)  #########
        self.emissions_dict = {}
        self.emissions_df = emission_df
        for s in range(self.n_emissions):
            self.emissions_dict[s] = emissions[s]

        #pick a random starting state (0,1,2,3) and initialize the matrices for the Viterbi algorithm
        self.current_state = rn.randint(0, self.n_states-1)
        self.optimal_matrix = np.zeros((self.n_states, len(SNPseq)))
        self.backtrack = np.zeros((self.n_states, len(SNPseq)))

        ### NOTE: still discussing how to represent this information ########
        #read in the SNPseq into a list
        self.text = SNPseq
        self.alphabet = emissions


    ### TO DO: test that this works with easy data!
    def get_emissions_matrix(self,position):
        """
        This function takes a SNP position (on chr 21), and uses the emission_df
        to output a 4x16 emission probability matrix (rows = 4 ancestry pop combos,
        cols = 16 genotype combos). The genotype combos are ordered lexicographically:
        ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"].
        """
        #first grab the index of position in emission_df
        ind = list(self.emission_df.loc[:,"POS"]).index(position) ######## might change this

        #now go through the columns with pop_nt freqs and make the emissions matrix
        emissions_mat = np.zeros((4,16),dtype=float)
        for i in range(self.n_states):
            pop1 = self.pop_codes[self.states[i][0]]
            pop2 = self.pop_codes[self.states[i][1]]
            for j in range(len(self.alphabet)):
                nt1 = self.alphabet[j][0]
                nt2 = self.alphabet[j][1]
                prob1 = self.emission_df.loc[ind,pop1+"_"+nt1]
                prob2 = self.emission_df.loc[ind,pop2+"_"+nt2]
                emissions_mat[i,j] = prob1*prob2

        return emissions_mat

    ### TO DO: update this for new DSs
    def get_optimal_path(self):
        for i, c in enumerate(self.text): ###
            idx = self.alphabet.index(c)
            for j in range(self.n_states):
                emit_prob = self.emissions_dict[j][idx] ###
                if i == 0:
                    self.optimal_matrix[j, i] =  emit_prob
                    self.backtrack[j, i] = -1
                else:
                    probs = []
                    for previous in range(self.n_states):
                        transition_prob = self.state_dict[previous][j]
                        probs.append(emit_prob * transition_prob * self.optimal_matrix[previous, i-1])
                    max_prob = max(probs)
                    back = probs.index(max_prob) ### think carefully about indexing here
                    self.optimal_matrix[j, i] = max_prob
                    self.backtrack[j, i] = back

    ### TO DO: update this for new DSs
    def reconstruct_path(self, state_names):
        end_values = [self.optimal_matrix[v, len(self.text) - 1] for v in range(self.n_states)]
        max_prob = max(end_values)
        end_state = end_values.index(max_prob)
        path = [end_state]
        previous_state = end_state
        for i in range(len(self.text)-1, -1, -1): ###?
            #print(previous_state, i)
            path.append(int(self.backtrack[previous_state, i]))
            previous_state = int(self.backtrack[previous_state, i])
        path.reverse()
        output = [state_names[p] for p in path if p != -1]
        return output


### big changes necessary
    # - text --> SNPseq (with SNP position info too)
    # - can hard code alphabet to be all 2nt combos?
    # - think more about how we want to represent state_names and state_probs
def parse_hmm_data(data):
    text = data[0] ###
    alphabet = data[2].split(" ") ###
    state_names = data[4].split(" ") ###
    n_states = len(state_names)
    state_probs = []
    for i in range(7, 7+n_states):
        state_probs.append([float(x) for x in data[i].split("\t")[1:]])
    emit_probs = []
    for j in range(9+n_states, 9+n_states*2):
        emit_probs.append([float(x) for x in data[j].split("\t")[1:]])
    return state_probs, emit_probs, text, alphabet, state_names
