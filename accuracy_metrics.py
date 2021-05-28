import numpy as np
import pandas as pd


def accuracy(predict_geno, truth_geno):
    """
    Calculate the overall accuracy of predictions. Both chromosomes are considered independently.
    :param predict_geno: Predicted genotype table
    :param truth_geno: Truth genotype table
    :return: Percent of positions correct across all positions in both chromosomes
    """
    n_pos = len(predict_geno)
    p1 = predict_geno.POP1.astype(str).to_numpy()
    p2 = predict_geno.POP2.astype(str).to_numpy()
    t1 = truth_geno.POP1.astype(str).to_numpy()
    t2 = truth_geno.POP2.astype(str).to_numpy()
    #print(type(t2))
    correct1 = np.equal(p1, t1).sum()
    correct2 = np.equal(p2, t2).sum()
    accuracy = (correct1 + correct2) / (2 * n_pos)
    return accuracy


def pivot_accuracy(predict_geno, truth_geno):
    """
    Accuracy metrics calculated from comparing the number and position of pivots between two genomes
    :param predict_geno: Predicted genome table
    :param truth_geno: Truth genome table
    :return: 1. percent_correct : DO NOT USE
            2. skew: ratio of predicted pivots compared to expected number of pivots
    """
    n_pos = len(predict_geno)
    p1 = predict_geno.POP1.astype(str).to_numpy()
    p2 = predict_geno.POP2.astype(str).to_numpy()
    t1 = truth_geno.POP1.astype(str).to_numpy()
    t2 = truth_geno.POP2.astype(str).to_numpy()
    pivots1 = find_pivots(np.equal(p1, t1))
    pivots2 = find_pivots(np.equal(p2, t2))
    # TODO some way to scale the number of correct pivots?
    # TODO incorrect pivots vs number of pivots??
    # Should it be compared to truth or predicted? Percent wrong OR
    # number wrong vs expected number?
    # Should it matter how close it is to the actual point? Nah that's what accuracy is for?
    total1 = find_pivots(np.equal(p1, p1[0])) + 1 * (p1[0] != t1[0])
    total2 = find_pivots(np.equal(p2, p2[0])) + 1 * (p2[0] != t2[0])
    # Get total percent of correct switches in the predicted genome
    if total1 + total2 == 0:
        percent_correct = 1
    else:
        percent_correct = 1 - (pivots1 + pivots2) / (total1 + total2)
    # Get skew
    expected1 = find_pivots(np.equal(t1, t1[0]))
    expected2 = find_pivots(np.equal(t2, t2[0]))
    if expected1 + expected2 == 0:
        skew = (total1 + total2)/2
    else:
        skew = (total1 + total2) / (expected1 + expected2)
    return percent_correct, skew


def find_pivots(bool_array):
    """
    Given a boolean array, count the number of times it switches True -> False and False -> True assuming start in True
    e.g. T,T,T,T,T,T,F,F,F,F,F,T,T,T,T = 2 pivots
    e.g. F,F,F,F,T,T,T,T,T,T,F,F,F,F,F = 3 pivots (counting initial pivot to False)
    e.g. T,T,T,T,T,T,T,T = 0 pivots
    e.g. F,F,F,F,F,F,F,F = 1 pivot
    :param bool_array: boolean array of any length
    :return: The number of pivots
    """
    idx = 0
    count = 0
    while idx < len(bool_array)-1:
        falses = np.where(bool_array[idx:] == False)
        if len(falses[0]) == 0:
            break
        else:
            count += 1
            idx += falses[0][0]
        trues = np.where(bool_array[idx:] == True)
        if len(trues[0]) == 0:
            break
        else:
            count += 1
            idx += trues[0][0]
    return count


def compare_genome_files(path1, path2):
    """
    :param path1: Path to predicted genome file
    :param path2: Path to truth genome file
    :return: Three skew metrics:
        1. bp_accuracy: percent of base pairs assigned correct ancestry
        2. correct_crossovers: DO NOT USE
        3. skew: Does the predict accuracy have fewer or more crossovers than expected.
            Skew < 1 = fewer, skew > 1 = more
    """
    a = pd.read_csv(path1, sep="\t", index_col=False)
    b = pd.read_csv(path2, sep="\t", index_col=False)
    bp_accuracy = accuracy(a, b)
    correct_crossovers, skew = pivot_accuracy(a, b)
    return bp_accuracy, correct_crossovers, skew


if __name__ == "__main__":
    p1 = "Data/simAdmixedGenome_100_0_1_Rx2_1.tsv"
    p2 = "Data/simGenome_100_0_1.tsv"
    x, y, z, = compare_genome_files(p2, p1)
    print(x)