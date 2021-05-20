import random as rn
import pandas as pd
from os import path
import re as re
import numpy as np

# TODO make it possible to save simulated data outside of Git


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
    data = freq_one_row_per_position(data)
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
    vals = data_sub.apply(lambda x: select_allele(x), axis=1)
    data_sub = data_sub.assign(select=vals)
    G = "".join(data_sub.select)
    return G


def select_allele(row):
    """
    Subroutine of simulate_chrom() to select an allele based on frequencies.
    :param row: row of data frame
    :return: A single selected allele
    """
    probs = np.array([row["A"], row["C"], row["G"], row["T"]])
    cumulative_probs = np.cumsum(probs)
    allele = ["A", "C", "G", "T"][np.argmax(cumulative_probs > rn.random())]
    return allele


def simulate_random_genome(ancestry_pair, data, save=False):
    """

    :param ancestry_pair:
    :param data:
    :param save:
    :return:
    """
    if data.shape[1] != 6:
        data = freq_collapse_pop_to_column(data)
    g1 = simulate_chrom(ancestry_pair[0], data)
    g2 = simulate_chrom(ancestry_pair[1], data)
    G = [g1[i]+g2[i] for i in range(len(g1))]
    genome = genome_to_table(ancestry_pair, G, data)
    if save:
        n_pos = data.POS.nunique()
        genome.to_csv("Data/simGenome_"+str(n_pos)+"_"+str(ancestry_pair[0])+"_"+str(ancestry_pair[1])+".tsv", sep="\t",
                      index=None)
    return G


def simulate_perfect_genomes(ancestries, data, save=False):
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
        vals = data_sub.filter(items=["A", "C", "G", "T"]).idxmax(axis=1)
        data_sub = data_sub.assign(select=vals)
        G.append(data_sub.select.values)
    diploid = {}
    for i, indiv in enumerate(ancestries):
        for j, indiv2 in enumerate(ancestries):
            diploid[(indiv, indiv2)] = [G[i][k] + G[j][k] for k in range(len(G[i]))]

    if save:
        save_perfect_genomes(diploid, data)
    return diploid


def genome_to_table(ancestry_pair, genome, data):
    """

    :param ancestry_pair:
    :param genome:
    :param data:
    :return:
    """
    genos = genome
    c1 = [x[0] for x in genos]
    c2 = [x[1] for x in genos]
    genome_table = pd.DataFrame({"POS": data.POS.unique(), "A1": c1, "A2": c2, "POP1": ancestry_pair[0],
                                 "POP2": ancestry_pair[1]})
    return genome_table


def table_to_genome(table):
    """

    :param table:
    :return:
    """
    if type(table) == str:
        table = pd.read_csv(table, sep = "\t", index_col=False)
    c1 = table["A1"]
    c2 = table["A2"]
    G = [c1[i] + c2[i] for i in range(len(c1))]
    return G


def save_perfect_genomes(all_genomes, data):
    """

    :param all_genomes:
    :param data:
    :return:
    """
    for pair in all_genomes.keys():
        tab = genome_to_table(pair, all_genomes[pair], data)
        n_pos = tab.POS.nunique()
        tab.to_csv("Data/simGenome_perfect_"+str(n_pos)+"_"+str(pair[0])+"_"+str(pair[1])+".tsv", sep="\t",
                   index=None)


def simulate_admixed_chrom(ancestry_pair, data, n_recomb):
    """
    Performs specified number of recombinations at randomly selected locations from two single ancestry chromosomes
    :param ancestry_pair: Pair of ancestries to admix
    :param data: Frequency data
    :param n_recomb: Number of recombination events
    :return: A string representing a single admixed chromosome (selected randomly from the two)
    """
    if data.shape[1] != 6:
        data = freq_collapse_pop_to_column(data)
    c1 = simulate_chrom(ancestry_pair[0], data)
    pop1 = [ancestry_pair[0] for i in range(len(c1))]
    c2 = simulate_chrom(ancestry_pair[1], data)
    pop2 = [ancestry_pair[1] for i in range(len(c2))]
    n_pos = data.POS.nunique()
    crossovers = [rn.randint(0, n_pos-1) for _ in range(n_recomb)]
    crossovers.sort()
    for cross in crossovers:
        c1, c2 = crossover(c1, c2, cross)
        pop1, pop2 = crossover(pop1, pop2, cross)
    chroms = [c1, c2]
    pops = [pop1, pop2]
    choice = rn.randint(0, 1)
    return chroms[choice], pops[choice]


def simulate_admixed_genome(ancestry_pair, data, n_recomb, sample_id="", save=False):
    """
    :param ancestry_pair: Names of the two ancestries to admix
    :param data: Allele frequency table
    :param n_recomb: Number of recombination events
    :param sample_id: Optional suffix for save file
    :param save: Should the generated genome be saved
    :return: Genome in format ["AA", "TA", ....]
    """
    a1, pop1 = simulate_admixed_chrom(ancestry_pair, data, n_recomb)
    a2, pop2 = simulate_admixed_chrom(ancestry_pair, data, n_recomb)
    genome_table = pd.DataFrame({"POS": data.POS.unique(), "A1": list(a1), "A2": list(a2), "POP1": pop1,
                                 "POP2": pop2})
    if save:
        filename = "Data/simAdmixedGenome_"+str(len(data))+"_"+"_".join(ancestry_pair)+"_Rx"+str(n_recomb)+"_"+sample_id+".tsv"
        genome_table.to_csv(filename, index=False, sep="\t")
    return table_to_genome(genome_table)


def crossover(a, b, cross_pos):
    """
    Performs a suffix exchange between two objects of equal length of the suffix from index cross_pos to end
    :param a: List or string to be crossed
    :param b: List or string to be crossed
    :param cross_pos: Position at which to perform crossover
    :return: both chromosomes following crossover
    """
    a_pref = a[:cross_pos]
    a_suff = a[cross_pos:]
    b_pref = b[:cross_pos]
    b_suff = b[cross_pos:]
    a_out = a_pref + b_suff
    b_out = b_pref + a_suff
    return a_out, b_out



if __name__ == "__main__":
    data = make_fake_frequencies(["0", "1"], 10, 520, save=False)
    x, y = simulate_admixed_genome(["0", "1"], data, 2)

    perf = simulate_perfect_genomes(["0", "1"], data, save=False)
    save_genome(G, data, "affga")