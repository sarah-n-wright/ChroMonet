import random as rn
import pandas as pd
import numpy as np


def get_fake_frequencies(ancestries, positions, seed=0, save=False):
    """
    Generates a fake set of allele frequencies. Two alleles for each position are chosen at random to be ref/alt.
    The reference allele frequency is a random number between 0 and 1.
    :param ancestries: Either integer number of ancestries, or list of specific ancestries
    :param positions: Either integer number of positions, or list of positions
    :param seed: Random seen if wanted. Default=no random seed
    :param save: Should data be saved to file? Default=False
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
        alleles = rn.sample(nts, 2)
        for pop in positions:
            freq = rn.random()
            data = data.append(pd.Series({"POS": pos,  "POP": pop, "REF": alleles[0],
                                          "ALT": alleles[1], "ref_freq": round(freq, 4)}),
                               ignore_index=True)
    data["alt_freq"] = 1 - data.ref_freq
    if save:
        data.to_csv("Data/simData_N"+str(len(ancestries))+"_P"+str(len(positions))+"_seed" + str(seed)+".tsv", sep="\t",
                    index=False)
    return data



if __name__ == "__main__:":
    x = get_fake_frequencies(2, 20, 20, save=True)
    x.head()

