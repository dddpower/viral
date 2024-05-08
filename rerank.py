"""
# reranking code for the result
## get the probability of RNA(nucleotide, codon) mutation
## basic formula is aX + bY + cZ

|Probability||||
|---|---|---|---|
A to C	1.1309cZ882976023E-05
A to G	3.79592984444666E-05
A to T	2.17601073885478E-05
C to A	4.82364467163978E-05
C to G	1.43123713007997E-05
C to T	0.000615014695634072
G to A	0.000111490944193118
G to C	7.0550851521742E-05
G to T	0.000311342109361459
T to A	7.57183086181891E-06
T to C	5.12463411619656E-05
T to G	6.88348260165355E-06

A to C	0.000480370596901429
A to G	0.00161226521002784
A to T	0.000924228464347168
C to A	0.000214773685034829
C to G	6.37261020475316E-05
C to T	0.00273836448419408
G to A	0.000681726831066903
G to C	0.000431392960074668
G to T	0.00190374448013487
T to A	0.00034768052837326
T to C	0.00235311053529112
T to G	0.000316073207612054
"""

from itertools import permutations
from Bio.Data import CodonTable
import pandas as pd
import argparse
import json
from enum import Enum


class Virus(Enum):
    COV = 1
    FLU = 2


def get_aa_from_position(codon_df, position):
    return codon_df.iloc[position]["1AA"]


def get_codon_from_position(codon_df, position):
    return codon_df.iloc[position]["Codon"]


# receive amino_acid tranlate it to 3char codon
def to_codons(ammino_acid):
    std_codon_table = CodonTable.unambiguous_dna_by_id[1]
    codons = [
        codon
        for codon, aa in std_codon_table.forward_table.items()
        if aa == ammino_acid
    ]
    assert len(codons) > 0, f"ammino_acid = {ammino_acid}"
    # print(f'{ammino_acid} is mapped to {codons}')
    return codons


"""
Example
prob of from wt(codon) to mutaion(ammino acid) 
like AGG => 'C', AGG => "D", ....
"""


def get_codon_aa_mutation_independent_probability(wt_codon, mut_aa, p_table):
    # 3 char for each codon
    # probability of mutation from wildtype based on p_table
    def get_codon_independent_probability(wt_codon, mut_codon, p_table):
        p = 1
        for ch1, ch2 in zip(wt_codon, mut_codon):
            p *= p_table[ch1 + ch2]
        return p

    mut_codons = to_codons(mut_aa)
    prob_list = [
        get_codon_independent_probability(wt_codon, mut_codon, p_table)
        for mut_codon in mut_codons
    ]
    assert (
        len(prob_list) > 0
    ), f"changed from {wt_codon}, to {mut_codons}, mut aa is {mut_aa}"
    return sum(prob_list)


def wt_codon(position, wt_codon_df):
    try:
        codon = wt_codon_df[wt_codon_df["pos"] == position]["Codon"].values[0]
    except Exception as e:
        print(f"{e}. Error position is {position}")
    return codon


def load_dict(file):
    return json.load(open(file))


def process_viral_result_table(virus_type: Virus, weight_func):
    h1n1_df = pd.read_csv(
        "results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", delimiter="\t"
    )
    cov_df = pd.read_csv(
        "results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", delimiter="\t"
    )
    cov_wt_codon_df = pd.read_csv("COV_wildtype_codon.csv")
    h1n1_wt_codon_df = pd.read_csv("H1N1_wildtype_codon.csv")

    rank_df, p_table, wt_codon_df = (
        (cov_df, load_dict("SARS-COV2-static-probability.json"), cov_wt_codon_df)
        if virus_type == Virus.COV
        else (h1n1_df, load_dict("H1N1-static-probability.json"), h1n1_wt_codon_df)
    )

    # move position -1 to sync with original table index
    wt_codon_df["pos"] -= 1
    # get rid of junklike alphabet
    rank_df = rank_df[~rank_df["wt"].str.contains("X|B|Z|J|U", case=False, na=False)]
    rank_df = rank_df[~rank_df["mut"].str.contains("X|B|Z|J|U", case=False, na=False)]

    # translate aa to codon
    rank_df["codon"] = rank_df.apply(
        lambda row: wt_codon(row["pos"], wt_codon_df), axis=1
    )

    # get wild type codon to mutation probability
    rank_df["codon_prob"] = rank_df.apply(
        lambda row: weight_func(wt_codon(row["pos"], wt_codon_df), row["mut"], p_table),
        axis=1,
    )

    # make rank columns
    rank_df["grammar_rank"] = rank_df["prob"].rank(method="min", ascending=False)
    rank_df["semantic_rank"] = rank_df["change"].rank(method="min", ascending=False)
    rank_df["codon_mut_rank"] = rank_df["codon_prob"].rank(
        method="min", ascending=False
    )
    rank_df["rank_sum"] = (
        rank_df["grammar_rank"] + rank_df["semantic_rank"] + rank_df["codon_mut_rank"]
    )
    rank_df["total_rank"] = rank_df["rank_sum"].rank(method="min", ascending=True)
    return rank_df


# need to be careful of boundary exception
def get_codon_aa_mutation_dependent_weight(codon_df, starting_position, mut_aa):
    """
    Get probability-like weight value related to chance from wildtype codon to
    some specific codon (mutation).
    Source contains prefix + wildtype + postfix (5 nucliotide chars)
    """

    def get_dependent_weights(source, mut_codon, reduction_func) -> list:
        """
        Current state contains from 0 to 4 position.
        Permutation number matches for
        0 => [0:3], mut_codon[0]
        1 => [1:4], mut_codon[1]
        2 => [2:5], mut_codon[2].
        """
        mut_dependent_prob_table = pd.read_excel(
            "12276_2021_658_MOESM2_ESM.xlsx", sheet_name=3
        )

        def calculate_weight(current, remain_permutation_number):
            # target_str format: "XX>XX"
            def get_mutation_dependent_probability(target_str):
                if target_str[1] == target_str[3]:
                    return 1
                sum = mut_dependent_prob_table["SARS-CoV-2"].sum()
                val = mut_dependent_prob_table[
                    mut_dependent_prob_table["Substition type"] == target_str
                ]["SARS-CoV-2"].values[0]
                return val / sum

            if not remain_permutation_number:
                return 1
            head = remain_permutation_number[0]
            keyword = (
                current[head]
                + current[head + 1]
                + ">"
                + mut_codon[head]
                + current[head + 2]
            )
            p = get_mutation_dependent_probability(keyword)
            current[head + 1] = mut_codon[head]
            return p * calculate_weight(current, remain_permutation_number[head + 1 :])

        return reduction_func(
            calculate_weight(source, numbers) for numbers in permutations([0, 1, 2], 3)
        )

    assert (
        get_aa_from_position(codon_df, starting_position) != mut_aa
    ), f"starting position = {starting_position}, mut_aa = {mut_aa}"
    mut_codons = to_codons(mut_aa)

    pre = starting_position - 1
    post = starting_position + 1
    prefix = get_codon_from_position(codon_df, pre)[2]
    postfix = get_codon_from_position(codon_df, post)[0]
    wt_codon = get_codon_from_position(codon_df, starting_position)
    probs = [
        sum(get_dependent_weights(prefix + wt_codon + postfix, mut_codon))
        for mut_codon in mut_codons
    ]
    return max(probs)


# cov_df = pd.read_csv("results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", delimiter='\t')
#
# flu_df = pd.read_csv("results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", delimiter='\t')
# rank_df = cov_df
# # get rid of start and last position
# max_pos = rank_df['pos'].max()
# rank_df = rank_df[(rank_df['pos'] != 0) & (rank_df['pos'] != max_pos)]
#
# # get rid of junklike alphabet
# rank_df = rank_df[~rank_df['wt'].str.contains('X|B|Z|J|U', case=False, na=False)]
# rank_df = rank_df[~rank_df['mut'].str.contains('X|B|Z|J|U', case=False, na=False)]
#
#
#
# # translate aa to codon
# rank_df['codon'] = rank_df.apply(lambda row: wt_codon(row['pos'], cov_wt_codon_df), axis=1)
# # get wild type codon to mutation probability
# rank_df['codon_prob'] = rank_df.apply(
#     lambda row: aa_mutation_dependent_probability(cov_wt_codon_df, row['pos'], row['mut']), axis=1)
# rank_df
#
#
#
# # make rank columns
# rank_df['grammar_rank'] = rank_df['prob'].rank(method='min', ascending=False)
# rank_df['semantic_rank'] = rank_df['change'].rank(method='min', ascending=False)
# rank_df['codon_mut_rank'] = rank_df['codon_prob'].rank(method='min', ascending=False)
# rank_df['rank_sum'] = rank_df['grammar_rank'] + rank_df['codon_mut_rank'] + rank_df['semantic_rank']
# # rank_df['rank_sum'] = rank_df['codon_mut_rank']
# rank_df['total_rank'] = rank_df['rank_sum'].rank(method='min', ascending=True)
# rank_df
#
#
# df2 = rank_df[rank_df['is_escape'] == True]
# mean, std = df2['total_rank'].mean(), df2['total_rank'].std()
# print(f'total number of cases are {len(df)}')
# print(f'for escape mutants priority rank applying independent prob ways, mean = {mean}, std = {std}')
# df2.to_csv("dependent_prob_df.csv")


def get_codon_mutation_probability(Virus=None, indendent=True):
    return


if __name__ == "__main__":
    pass
