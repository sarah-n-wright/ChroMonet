import random as rn
import numpy as np

### TO DO:
    # - Refresh on genome representation, change indexing methods as necessary
    #   (just cut down emission_df and then use the index in SNPseq?)
    # - Adapt optimal path and reconstruct path for any changes and test!


def HammingDist(str1,str2):
    dist = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            dist += 1
    return dist



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


    #function which creates an emissions matrix for a specific SNP (based on index in SNPseq) ###(move this down)
    ### Need to adjust use of position depending on how we do our indexing (could be row # instead) ##########
    def get_emissions_matrix(self,position):
        """
        This function takes a SNP position (on chr 21 ####), and uses the emission_df
        to output a 4x16 emission probability matrix (rows = 4 ancestry pop combos,
        cols = 16 genotype combos). The genotype combos are ordered lexicographically:
        ["AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"].
        """
        #first grab the index of position in emission_df
        ind = list(self.emissions_df.loc[:,"POS"]).index(position) ######## might change this

        #now go through the columns with pop_nt freqs and make the emissions matrix
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
