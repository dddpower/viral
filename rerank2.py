from itertools import permutations
from Bio.Data import CodonTable
from collections.abc import Callable
import pandas as pd
import argparse
import json
from enum import Enum, auto


class Virus(Enum):
    COV = auto()
    FLU = auto()
    HIV = auto()


# receive amino_acid tranlate it to 3char codon
def aa_to_codons(ammino_acid):
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
Define calc_funtion for calculating codon mutation weight values here.
Signature of calc_func should be (virus_type: Virus, wildtype_position: int, target_codon: str) -> weight value: float
"""


"""
applied position independent mutation probability.
"""

def single_codon_mutation_independent_prob(
    virus_type, wildtype_position, target_codon
) -> float:
    # flu_cscs_df = pd.read_csv("results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", delimiter='\t')
    # cov_cscs_df = pd.read_csv("results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", delimiter='\t')
    wildtype_file, mut_prob_file = (
        ("cov_wildtype_codon.csv", "cov_static_probability.json")
        if virus_type == Virus.COV
        else ("flu_wildtype_codon.csv", "flu_static_probability.json")
    )
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

    aa_to_codons(target_aa)
    return sum(
        map(
            lambda x: single_codon_mutation_independent_prob(wt_codon, x, p_table),
            aa_to_codons(mut_aa),
        )
    )


def single_codon_mutation_conditional_weight(
    virus_type, wildtype_position, target_codon
) -> float:
    wildtype_file, mut_prob_file = (
        "cov_wildtype_codon.csv",
        "cov_static_probability.json",
    )

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


# need to be careful to handle boundary exception
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


def aa_to_codons(ammino_acid) -> tuple:
    std_codon_table = CodonTable.unambiguous_dna_by_id[1]
    codons = (
        codon
        for codon, aa in std_codon_table.forward_table.items()
        if aa == ammino_acid
    )
    assert len(codons) > 0, f"ammino_acid = {ammino_acid}"
    return codons


# we will make static file in the end.
def get_codon_mutation_probability(
    virus_type: Virus,
    source_position: int,
    target_codon: str,
    calc_func: Callable[[Virus, int, str], float],
) -> float:
    return calc_func(virus_type, source_position, target_codon)


def get_aa_mutation_weight(
    aa: str, virus_type: Virus, source_position, target_aa, calc_func
) -> float:
    return sum(
        map(
            lambda x: get_codon_mutation_probability(
                virus_type, source_position, x, calc_func
            ),
            aa_to_codons(target_aa),
        )
    )


def add_codon_mutation_rank(virus: Virus):
    return
