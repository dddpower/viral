import pandas as pd
from itertools import permutations
input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
prob_table = pd.read_excel(input_file_name, sheet_name=3)
wildtype_codon_info = pd.read_csv("cov_wildtype_codon.csv")
def get_position_dependent_codon_weight(source_position, mut_codon) -> float:
    wildtype_codon_info[source_position]
    # print(f'{get_position_dependent_codon_prob.__name__} called')

    """
    Current state contains from 0 to 4 position.
    Permutation number matches for
    0 => [0:3], mut_codon[0]
    1 => [1:4], mut_codon[1]
    2 => [2:5], mut_codon[2].
    """
    input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
    prob_table = pd.read_excel(input_file_name, sheet_name=3)

    def get_codon_weight(current, remain_permutation_number):
        # target_str format: "XX>XX"
        def get_single_prob(target_str):

            # source and mutant are same. return identity element for multiplication
            if rget_str[1] == target_str[3]:
                return 1
            sum = prob_table["SARS-CoV-2"].sum()
            val = prob_table[
                prob_table["Substition type"] == target_str
            ]["SARS-CoV-2"].values[0]
            return val / sum

        # identity element for multiplication
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
        p = get_single_prob(keyword)
        current[head + 1] = mut_codon[head]
        return p * get_codon_weight(current, remain_permutation_number[head + 1 :])

    source_codon = wildtype_codon_info[source_position:source_position + 3]
    return sum(
        get_codon_weight(source_codon, numbers) for numbers in permutations(range(0, 3), 3)
    )