from Bio.Data import CodonTable
import pandas as pd
from itertools import permutations, product
import json
import argparse
from enum import Enum, auto

# receive amino_acid and translate it to 3char codon
def to_codons(ammino_acid):
    std_codon_table = CodonTable.unambiguous_dna_by_id[1]
    codons = [codon for codon, aa in std_codon_table.forward_table.items() \
        if aa == ammino_acid]
    assert len(codons) > 0, f'ammino_acid = {ammino_acid}'
    return codons



# class Virus(Enum):
#     COV = auto()
#     FLU = auto()
#     HIV = auto()

if __name__ =="__main__":
    print("need to input virus type when excuting the program")

# Wildtype position table (Header: pos,Codon,1AA)
wt_positions = pd.read_csv("cov_wildtype_codon.csv")

# Probability table
with open("cov_static_probability.json", "r") as json_fs:
    pos_indep_prob_table = json.load(json_fs)

# Load position-dependent probability table
input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
pos_dep_prob_table = pd.read_excel(input_file_name, sheet_name=3)

def codon_from_aa_pos(source_position):
    return wt_positions[wt_positions['pos'] == source_position]['Codon'].values[0]


"""Let's build_independent table"""
def c2a_indep_prob(source_position, target:str):
    """wildtype codon to mutant amino_acid independent probability"""
    source_codon = codon_from_aa_pos(source_position)
    target_codons = to_codons(target)
    target_codons = filter(lambda x: x != source_codon, target_codons)
    def sum_iter():
        for target_codon in target_codons:
            p = 1
            for s, t in zip(source_codon, target_codon):
                p *= pos_indep_prob_table[s + t]
            yield p
        
    return sum(sum_iter())


def test_c2a_indep_prob(position):
    """test function for  c2a_indep_prob"""
    all_amino_acids = "ARNDCEQGHILKMFPSTWYV"
    return map(lambda x: c2a_indep_prob(position, x), (_ for _ in all_amino_acids))


# must ignore first 'M'
def pos_dep_codon_weight(aa_source_position, mut_codon):
    """Compute the position-dependent mutation weight"""
    def codon_weight(current: str , remain_permutation_number):
        def single_prob(keyword):
            sum = pos_dep_prob_table["SARS-CoV-2"].sum()
            val = pos_dep_prob_table[
                pos_dep_prob_table["Substition type"] == keyword
            ]["SARS-CoV-2"].values[0]
            return val / sum

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
        p = 1 if keyword[1] == keyword[3] else single_prob(keyword)
        new_current = current[:head] + mut_codon[head - 1] + current[head + 1:]
        result = p * codon_weight(new_current, remain_permutation_number[head + 1 :])
        return 0 if result == 1 else result

    pre = codon_from_aa_pos(aa_source_position - 1)[2]
    codon = codon_from_aa_pos(aa_source_position)
    end = codon_from_aa_pos(aa_source_position + 1)[0]
    return sum(codon_weight(pre + codon + end, _) for _ in permutations(range(1, 4)))




# from 2 ~ until aa_pos
def test_pos_dep_codon_weight(aa_pos):
    def all_codons():
        codon_table = CodonTable.unambiguous_dna_by_id[1]
        return codon_table.forward_table.keys()


    for _ in range(2, aa_pos + 1):
        yield list(map(lambda x: pos_dep_codon_weight(aa_pos, x), all_codons()))
