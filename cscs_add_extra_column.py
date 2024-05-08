from Bio.Data import CodonTable
import pandas as pd
from itertools import permutations, combinations
import json
import argparse
from enum import Enum, auto

all_amino_acids = "ARNDCEQGHILKMFPSTWYV"


# receive amino_acid and translate it to 3char codon
def to_codons(ammino_acid):
    std_codon_table = CodonTable.unambiguous_dna_by_id[1]
    codons = [
        codon
        for codon, aa in std_codon_table.forward_table.items()
        if aa == ammino_acid
    ]
    assert len(codons) > 0, f"ammino_acid = {ammino_acid}"
    return codons


# class Virus(Enum):
#     COV = auto()
#     FLU = auto()
#     HIV = auto()


# Wildtype position table (Header: pos,Codon,1AA)
wt_positions = pd.read_csv("cov_wildtype_codon.csv")

# Probability table
with open("cov_static_probability.json", "r") as json_fs:
    pos_indep_prob_table = json.load(json_fs)

# Load position-dependent probability table
input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
pos_dep_prob_table = pd.read_excel(input_file_name, sheet_name=3)

def codon_from_aa_pos(source_position):
    return wt_positions[wt_positions["pos"] == source_position]["Codon"].values[0]


"""Build_independent table with funcions written below."""


def pos_indep_prob(aa_pos, target: str):
    """wildtype codon to mutant amino_acid independent probability"""
    source_codon = codon_from_aa_pos(aa_pos)
    target_codons = to_codons(target)
    target_codons = filter(lambda x: x != source_codon, target_codons)

    def sum_iter():
        for target_codon in target_codons:
            p = 1
            for s, t in zip(source_codon, target_codon):
                p *= pos_indep_prob_table[s + t]
            yield p

    return sum(sum_iter())


freq_sum = pos_dep_prob_table["SARS-CoV-2"].sum()
def pos_dep_weight(source_position, mut_aa):
    # must ignore first 'M'
    def pos_dep_codon_weight(source_position, mut_codon):
        """Compute the position-dependent mutation weight"""
        def codon_weight(current: str, remain_permutation_number):
            def single_prob(keyword):
                val = pos_dep_prob_table[
                    pos_dep_prob_table["Substition type"] == keyword
                ]["SARS-CoV-2"].values[0]
                return val / freq_sum

            if not remain_permutation_number:
                return 1

            head = remain_permutation_number[0]
            keyword = (
                current[head - 1]
                + current[head]
                + ">"
                + mut_codon[head - 1]
                + current[head + 1]
            )
            p = single_prob(keyword)
            # p = 1 if keyword[1] == keyword[3] else single_prob(keyword)
            new_current = current[:head] + mut_codon[head - 1] + current[head + 1 :]
            result = p * codon_weight(
                new_current, remain_permutation_number[1:]
            )
            return 0 if result == 1 else result

        pre = codon_from_aa_pos(source_position - 1)[2]
        codon = codon_from_aa_pos(source_position)
        end = codon_from_aa_pos(source_position + 1)[0]
        # diff = map(lambda x: x[0] != x[1], zip(codon, mut_codon))
        return sum(
            codon_weight(pre + codon + end, _) for _ in permutations(_ for _ in range(1, 4) if codon[_ - 1] != mut_codon[_ - 1])
        )

    return sum(
        map(lambda x: pos_dep_codon_weight(source_position, x), to_codons(mut_aa))
    )


"""test functions"""


# from 2 ~ until aa_pos
def test_pos_dep_weight(aa_pos):
    ...

def test_pos_indep_prob(position):
    """test function for  c2a_indep_prob"""
    return map(lambda x: pos_indep_prob(position, x), (_ for _ in all_amino_acids))


def make_extra_column():
    """add extra rank columns to existing (cscs) table"""
    cscs_path = "results/cov/semantics/analyze_semantics_cov_bilstm_512.txt"
    cscs = pd.read_csv(cscs_path, sep="\t")
    cscs["pos"] += 1

    """remove unusable elements"""
    # 1. get rid of row that matches condition "pos == 1"
    cscs = cscs[cscs["pos"] != 1]

    # 2. get rid of unusable alphabet
    cscs = cscs[
        ~cscs["wt"].str.contains("X|B|Z|J|U", case=False)
        & ~cscs["mut"].str.contains("X|B|Z|J|U", case=False)
    ]
    cscs["pos_dep_weight"] = cscs.apply(
        lambda row: pos_dep_weight(row["pos"], row["mut"]), axis=1
    )
    cscs["pos_indep_prob"] = cscs.apply(
        lambda row: pos_indep_prob(row["pos"], row["mut"]), axis=1
    )

    # make rank columns
    cscs["grammar rank"] = cscs["prob"].rank(method="min", ascending=False)
    cscs["semantic rank"] = cscs["change"].rank(method="min", ascending=False)
    cscs["pos_indep_rank"] = cscs["pos_indep_prob"].rank(method="min", ascending=False)
    cscs["pos_dep_rank"] = cscs["pos_dep_weight"].rank(method="min", ascending=False)

    base_list = ["grammar rank", "semantic rank"]
    for x, y in combinations(base_list + ["pos_indep_rank"], 2):
        cscs[x + "+" + y] = cscs[x] + cscs[y]
    for x, y in combinations(base_list + ["pos_dep_rank"], 2):
        cscs[x + "+" + y] = cscs[x] + cscs[y]

    name1 = "grammar rank"
    name2 = "semantic rank"
    name3 = "pos_indep_rank"
    name4 = "pos_dep_rank"
    sum_rank0 = "+".join([name1, name2])
    sum_rank1 = "+".join([name1, name2, name3])
    sum_rank2 = "+".join([name1, name2, name4])
    cscs[sum_rank0] = cscs[name1] + cscs[name2]
    cscs[sum_rank1] = cscs[name1] + cscs[name2] + cscs[name3]
    cscs[sum_rank2] = cscs[name1] + cscs[name2] + cscs[name4]
    cscs[sum_rank0 + " rank"] = cscs[sum_rank0].rank(method="min", ascending=True)
    cscs[sum_rank1 + " rank"] = cscs[sum_rank1].rank(method="min", ascending=True)
    cscs[sum_rank2 + " rank"] = cscs[sum_rank2].rank(method="min", ascending=True)
    # cscs["rank_sum"] = (
    #     cscs["grammar_rank"] + cscs["semantic_rank"] + cscs["codon_mut_rank"]
    # )
    # cscs["total_rank"] = cscs["rank_sum"].rank(method="min", ascending=True)
    cscs.to_csv("0507.csv", index=False)


if __name__ == "__main__":
    make_extra_column()