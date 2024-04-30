import pandas as pd
import argparse
from itertools import permutations

def add_nucleotide_mutant_rank_column(
    new_column_name,
    cscs_table: pd.DataFrame,
    wildtype_codon_info,
    conditional_prob_table,
    nonconditional_prob_table
) -> pd.DataFrame:
    def _args():
        return (
            new_column_name,
            cscs_table,
            wildtype_codon_info,
            conditional_prob_table,
            nonconditional_prob_table
        )

    # Independent, dependent mechanism are quite different but returns same type(float).



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="parse arguments for adding extra column to existing dataframe"
    )
    parser.add_argument("column_name", type=str, help="column name to be added")
    parser.add_argument("cscs_table_path")
    parser.add_argument(
        "wildtype_info",
        help="path of wildtype codon positional information table",
        type=str,
    )
    parser.add_argument(
        "non_cond_prob",
        type=str,
        help="path of position-independent mutation probability table",
    )
    parser.add_argument(
        "cond_prob",
        type=str,
        help="path of position-dependent mutation probability table",
    )

    args = parser.parse_args()
    # print(cscs_add_nucleotide_mutant_rank_column(*vars(args).values()))
