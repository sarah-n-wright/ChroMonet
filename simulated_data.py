import random as rn
import pandas as pd
from os import path

# TODO save structure for genomes after deciding on data structures


def make_fake_frequencies(ancestries, positions, seed=0, save=False):
    """
    Generates a fake set of allele frequencies. Two alleles for each position are chosen at random to be ref/alt.
    The reference allele frequency is a random number between 0 and 1.
    :param ancestries: Either integer number of ancestries, or list of specific ancestries
    :param positions: Either integer number of positions, or list of positions
    :param seed: Random seen if wanted. Default=no random seed
    :param save: Should data be saved to file? Default=False. If True data will be saved as
        Data/simData_N[#ancestries]_P[#positions]_seed[seed_value].tsv
    :return: Data frame with position, ancestry, alleles and frequency. Each row is a single ancestry:position pair
    with both alleles
    """
    if seed != 0:
        rn.seed(seed)
    if type(ancestries) == int:
        ancestries = range(ancestries)
    if type(positions) == int:
        positions = range(positions)
    nts = ['A', 'G', 'C', 'T']
    data = pd.DataFrame({"POS": [], "POP": [], "REF": [], "ALT": [], "ref_freq": []})
    for pos in positions:
        alleles = rn.sample(nts, 2)  # randomly select the two alleles
        for pop in ancestries:
            freq = rn.random()  # randomly set the ref allele frequency
            data = data.append(pd.Series({"POS": pos,  "POP": pop, "REF": alleles[0],
                                          "ALT": alleles[1], "ref_freq": round(freq, 4)}),
                               ignore_index=True)
    data["alt_freq"] = 1 - data.ref_freq  # calculate the alt frequency
    if save:
        data.to_csv("Data/simData_N"+str(len(ancestries))+"_P"+str(len(positions))+"_seed" + str(seed)+".tsv", sep="\t",
                    index=False)
    return data


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

