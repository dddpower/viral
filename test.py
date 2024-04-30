from Bio.Data import CodonTable
import pandas as pd
from itertools import permutations
import json

# receive amino_acid and translate it to 3char codon
def to_codons(ammino_acid):
    std_codon_table = CodonTable.unambiguous_dna_by_id[1]
    codons = [codon for codon, aa in std_codon_table.forward_table.items() \
        if aa == ammino_acid]
    assert len(codons) > 0, f'ammino_acid = {ammino_acid}'
    # print(f'{ammino_acid} is mapped to {codons}')
    return codons
# 
# # testing to_codons()
# def list_all_to_codons():
#     aas = list("ARNDCEQGHILKMFPSTWYV")
#     return list(map(to_codons, aas))

def single_codon_mutation_independent_prob(
    virus_type, wildtype_position, target_codon
) -> float:
    # flu_cscs_df = pd.read_csv("results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", delimiter='\t')
    # cov_cscs_df = pd.read_csv("results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", delimiter='\t')
    wildtype_file, mut_prob_file = ("cov_wildtype_codon.csv", "cov_static_probability.json")
    wildtype_df = pd.read_csv(wildtype_file)
    mut_prob_table = json.load(open(mut_prob_file))
    p = 1
    for ch1, ch2 in zip(wildtype_df, target_codon):
        p *= mut_prob_table[ch1 + ch2]
    return p

def codon_aa_mutation_independent_prob(
    wildtype_codon_table, prob_table, wildtype_position, target_aa
):
    def get_prob_from_table(source: str, dest: str, table: dict):
        return table[source + dest]

    to_codons(target_aa)
    return sum(
        map(
            lambda x: single_codon_mutation_independent_prob(wt_codon, x, p_table),
            aa_to_codons(mut_aa),
        )
    )




cov_wildtype_codon_info = pd.read_csv("cov_wildtype_codon.csv")
input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
cov_prob_table = pd.read_excel(input_file_name, sheet_name=3)
# must ignore first 'M'
def position_dependent_codon_weight(aa_source_position, mut_codon, prob_table):

    def codon_from_aa_pos(aa_pos, wildtype_codon_info):
        return wildtype_codon_info[wildtype_codon_info["pos"] == aa_pos]["Codon"].values[0]


    def codon_weight(current: str , remain_permutation_number):


        def single_prob(target_str, prob_table):
            # source and mutant are same. return identity element for multiplication
            if target_str[1] == target_str[3]:
                return 1

            sum = prob_table["SARS-CoV-2"].sum()
            val = prob_table[
                prob_table["Substition type"] == target_str
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
        p = single_prob(keyword, prob_table)
        new_current = current[:head] + mut_codon[head - 1] + current[head + 1:]
        result = p * codon_weight(new_current, remain_permutation_number[head + 1 :])
        return 0 if result == 1 else result

    pre = codon_from_aa_pos(aa_source_position - 1, cov_wildtype_codon_info)[2]
    codon = codon_from_aa_pos(aa_source_position, cov_wildtype_codon_info)
    end = codon_from_aa_pos(aa_source_position + 1, cov_wildtype_codon_info)[0]
    return sum(codon_weight(pre + codon + end, _) for _ in permutations(range(1, 4)))


def all_codons():
    codon_table = CodonTable.unambiguous_dna_by_id[1]
    return codon_table.forward_table.keys()


# from 2 ~ until aa_pos
def test_position_dependent_codon_weight(aa_pos):
    for _ in range(2, aa_pos + 1):
        yield list(map(lambda x: position_dependent_codon_weight(aa_pos, x, cov_prob_table), get_all_codons()))

if __name__ == "__main__":
    print("entered main scope")