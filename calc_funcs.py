"""
Define calc_funtion for calculating codon mutation weight values here.
Signature of calc_func should be (virus_type: Virus, wildtype_position: int, target_codon: str) -> weight value: float
"""

"""
applied position independent mutation probability.
"""
import pandas as pd
def get_codon_independent_probability(virus_type, wildtype_position, target_codon) -> float:
    # flu_cscs_df = pd.read_csv("results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", delimiter='\t')
    # cov_cscs_df = pd.read_csv("results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", delimiter='\t')
    wild_type_codon_file = "COV_wildtype_codon.csv" if virus_type == 
    cov_wt_codon_df = pd.read_csv("COV_wildtype_codon.csv")
    flu_wt_codon_df = pd.read_csv("H1N1_wildtype_codon.csv")
    wildtype_codon_df
        p = 1
        for ch1, ch2 in zip(wt_codon, mut_codon):
            p *= p_table[ch1 + ch2]
        return p

def get_codon_aa_mutation_independent_probability(wt_codon, mut_aa, p_table):
    mut_codons = to_codons(mut_aa)
    prob_list = [get_codon_independent_probability(wt_codon, mut_codon, p_table) for mut_codon in mut_codons]
    assert len(prob_list) > 0, f"changed from {wt_codon}, to {mut_codons}, mut aa is {mut_aa}"
    return sum(prob_list)
