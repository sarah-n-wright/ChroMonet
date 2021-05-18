import random as rn
import pandas as pd
from os import path
import re as re

# TODO Modify data table for ref data
# TODO Update methods to work with new ref data frame
# TODO Save


def make_fake_frequencies(ancestries, positions, seed=0, save=False):
    """
    Generates a fake set of allele frequencies. Two alleles for each position are chosen at random to be ref/alt.
    The reference allele frequency is a random number between 0 and 1.
    :param ancestries: Either integer number of ancestries, or list of specific ancestries
    :param positions: Either integer number of positions, or list of positions
    :param seed: Random seen if wanted. Default=no random seed
    :param save: Should data be saved to file? Default=False. If True data will be saved as
        Data/simData_N[#ancestries]_P[#positions]_seed[seed_value].tsv
    :return: Data frame with position and frequencies for each nucleotide for each population
    """
    if seed != 0:
        rn.seed(seed)
    if type(ancestries) == int:
        ancestries = [int(x) for x in range(ancestries)]
    if type(positions) == int:
        positions = range(positions)
    nts = ['A', 'G', 'C', 'T']
    data = pd.DataFrame({"POS": [], "POP": [], "A": [], "C": [], "G": [], "T": []})
    for pos in positions:
        alleles = rn.sample(nts, 2)  # randomly select the two alleles
        zero_alleles = [x for x in nts if x not in alleles]
        for pop in ancestries:
            freq = rn.random()  # randomly set the ref allele frequency
            data = data.append(pd.Series({"POS": int(pos), "POP": int(pop), alleles[0]: round(freq, 4),
                                          alleles[1]: round(1-freq, 4), zero_alleles[0]: 0,
                                          zero_alleles[1]: 0}),
                               ignore_index=True)
    data["POS"] = data["POS"].astype(int)
    if data["POP"].dtype == float:
        data["POP"] = data["POP"].astype(int)
    if save:
        data.to_csv("Data/simData_N"+str(len(ancestries))+"_P"+str(len(positions))+"_seed" + str(seed)+".tsv", sep="\t",
                    index=False)
    return data


def freq_one_row_per_position(data):
    """
    Reformats ancestry representation of frequency data to one row per position
    :param data: frequency data with 6 columns : POS, POP, A, C, G, T
    :return: frequency data with one row for each position
    """
    unique_pops = data.POP.unique()
    pop_data = []
    for pop in unique_pops:
        pop_data.append(data.loc[(data.POP == pop), ("POS", "A", "C", "G", "T")])
    out_data = pop_data[0]
    for i in range(1, len(pop_data)):
        out_data = out_data.merge(pop_data[i], on="POS", suffixes=[unique_pops[i-1], unique_pops[i]])
    return out_data


def freq_collapse_pop_to_column(data, pop_is_suffix=True):
    """
    Reformats population representation of frequency data to one column for each nucleotide
    :param data: Frequency table with one row per position
    :param pop_is_suffix: Nucleotide column names are e.g. A_EUR. If false, columns are e.g EUR_A
    :return: Data with populations collapsed into a column, with one column per nucleotide
    """
    if pop_is_suffix:
        pops = set([re.sub("^[ACTG]", "", col) for col in data.columns[1:]])
    else:
        pops = set([re.sub("[ACTG]$", "", col) for col in data.columns[1:]])
    out_list = []
    for pop in pops:
        cols = [col for col in data.columns if re.search(pop, col)]
        d_pop = data.filter(items=["POS"]+cols)
        cols = [re.sub(pop, "", col) for col in cols]
        d_pop.columns = ["POS"] + cols
        d_pop["POP"] = pop
        out_list.append(d_pop)
    out_data = out_list[0]
    for i in range(1, len(pops)):
        out_data = pd.concat([out_data, out_list[i]])
    out_data.sort_values(by=["POS", "POP"], inplace=True)
    return out_data


def load_fake_frequencies(num_ancestries, num_positions, seed_num=0):
    """
    Checks if the dataset has been generated and loads simulated data based on parameters
    :param num_ancestries: Number of ancestries in the simulated data
    :param num_positions: Number of positions in the simulated data
    :param seed_num: Seed number used to generate data (0 if no seed used)
    :return: If file exists, returns dataframe of simulated data
    """
    filename = "Data/simData_N"+str(num_ancestries)+"_P"+str(num_positions)+"_seed"+str(seed_num)+".tsv"
    if path.exists(filename):
        data = pd.read_csv(filename, index_col=False,  sep="\t")
        return data
    else:
        print("File:", filename, "does not exist, please generate it first.")
        return None


def simulate_chrom(ancestry, data):
    """
    Simulates a chromosome based on given allele frequencies by using a random number generator to select the allele
    :param ancestry: Single ancestry identifier to match values in data
    :param data: Simulated data set to use for simulation
    :return: A string representing a single chromosome for an individual of input ancestry
    """
    if data.shape[1] != 6:
        data = freq_collapse_pop_to_column(data)
    data_sub = data.loc[(data.POP == ancestry), :]
    vals = data_sub.apply(lambda x: x.iloc[2+sum([rn.random() > x["ref_freq"]])], axis=1)
    data_sub = data_sub.assign(select=vals)
    G = "".join(data_sub.select)
    return G


def simulate_random_genomes(ancestry_pairs, data, N=1):
    """
    Simulate a set of genomes for given ancestries based on allele frequencies
    :param ancestry_pairs: List lists/tuples of ancestry pairs
    :param data: Simulated datasets
    :param N: Number of repeats
    :return: N x len(ancestry_pairs) diploid genomes as a dictionary where:
        keys = "[ancestry1][ancestry2].[repeat]"
        values = list of length len(Genome) of diploid genotypes e.g. ["AA", "TA", "GG" .....]
    """
    if data.shape[1] != 6:
        data = freq_collapse_pop_to_column(data)
    genomes = {}
    for n in range(N):
        for pair in ancestry_pairs:
            g1 = simulate_chrom(pair[0], data)
            g2 = simulate_chrom(pair[1], data)
            idd = "".join([str(p) for p in pair] + [".", str(n)])
            genos = [g1[i] + g2[i] for i in range(len(g1))]
            genomes[idd] = genos
    return genomes


def simulate_perfect_genomes(ancestries, data):
    """
    Simulate genomes where the most common allele for a given ancestry is always chosen
    :param ancestries: List of all available ancestries, all combinations will be made by the function
    :param data: Simulated dataset
    :return: A dictionary of simulate genomes for each combination of ancestries:
        keys = [ancestry1, ancestry2]
        values = list of length len(Genome) of diploid genotypes e.g. ["AA", "TA", "GG" .....]
    """

    if type(data) == str:
        data = pd.read_csv(data, sep="\t", index=None)
    if data.shape[1] != 6:
        data = freq_collapse_pop_to_column(data)
    G = []
    for indiv in ancestries:
        data_sub = data.loc[data.POP == indiv]
        vals = data_sub.apply(lambda x: x.iloc[2+sum([x.ref_freq < x.alt_freq])], axis=1)
        data_sub = data_sub.assign(select=vals)
        G.append("".join(data_sub.select.values))
    diploid = {}
    for i, indiv in enumerate(ancestries):
        for j, indiv2 in enumerate(ancestries):
            diploid[(indiv, indiv2)] = [G[i][k] + G[j][k] for k in range(len(G[i]))]
    return diploid


if __name__ == "__main__":
    data = make_fake_frequencies(2, 5, seed=0, save=False)
    x = cast_population_freq_table(data)
    d = melt_population_freq_table(x)
    x.head()