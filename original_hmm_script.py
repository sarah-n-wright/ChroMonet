import random as rn
import numpy as np

class HMMOptimalPath:
    def __init__(self, states, emissions, text, alphabet):
        self.n_states = len(states)
        self.n_emissions = len(emissions)
        self.state_dict = {}
        self.emissions_dict = {}
        for s in range(self.n_states):
            self.state_dict[s] = states[s]
            self.emissions_dict[s] = emissions[s]
        self.current_state = rn.randint(0, self.n_states-1)
        self.optimal_matrix = np.zeros((self.n_states, len(text)))
        self.backtrack = np.zeros((self.n_states, len(text)))
        self.text = text
        self.alphabet = alphabet


    def get_optimal_path(self):
        for i, c in enumerate(self.text):
            idx = self.alphabet.index(c)
            for j in range(self.n_states):
                emit_prob = self.emissions_dict[j][idx]
                if i == 0:
                    self.optimal_matrix[j, i] =  emit_prob
                    self.backtrack[j, i] = -1
                else:
                    probs = []
                    for previous in range(self.n_states):
                        transition_prob = self.state_dict[previous][j]
                        probs.append(emit_prob * transition_prob * self.optimal_matrix[previous, i-1])
                    max_prob = max(probs)
                    back = probs.index(max_prob)
                    self.optimal_matrix[j, i] = max_prob
                    self.backtrack[j, i] = back

    def reconstruct_path(self, state_names):
        end_values = [self.optimal_matrix[v, len(self.text) - 1] for v in range(self.n_states)]
        max_prob = max(end_values)
        end_state = end_values.index(max_prob)
        path = [end_state]
        previous_state = end_state
        for i in range(len(self.text)-1, -1, -1):
            #print(previous_state, i)
            path.append(int(self.backtrack[previous_state, i]))
            previous_state = int(self.backtrack[previous_state, i])
        path.reverse()
        output = [state_names[p] for p in path if p != -1]
        return output


def parse_hmm_data(data):
    text = data[0]
    alphabet = data[2].split(" ")
    state_names = data[4].split(" ")
    n_states = len(state_names)
    state_probs = []
    for i in range(7, 7+n_states):
        state_probs.append([float(x) for x in data[i].split("\t")[1:]])
    emit_probs = []
    for j in range(9+n_states, 9+n_states*2):
        emit_probs.append([float(x) for x in data[j].split("\t")[1:]])
    return state_probs, emit_probs, text, alphabet, state_names


if __name__ == '__main__':
    text = "xyxzzxyxyy"
    alphabet = ['x', 'y', 'z']
    st = [[0.641, 0.359], [0.729, 0.271]]
    em = [[0.177, 0.691, 0.192], [0.097, 0.42, 0.483]]
    hmm = HMMOptimalPath(st, em, text, alphabet)
    hmm.get_optimal_path()
    path = hmm.reconstruct_path(["A", "B"])
    assert "".join(path) == "AAABBAAAAA"