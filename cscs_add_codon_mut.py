from Bio.Data import CodonTable
import pandas as pd
from itertools import permutations
import json
import argparse
import os

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


def codon_from_aa_pos(wt_pos_df, source_position):
    return wt_pos_df[wt_pos_df["pos"] == source_position]["Codon"].values[0]


def pos_indep_prob(wt_pos_df, pos_indep_prob_dict, aa_pos, target: str):
    """wildtype codon to mutant amino_acid independent probability"""
    source_codon = codon_from_aa_pos(wt_pos_df, aa_pos)
    target_codons = filter(lambda x: x != source_codon, to_codons(target))

    def sum_iter():
        for target_codon in target_codons:
            p = 1
            for s, t in zip(source_codon, target_codon):
                p *= pos_indep_prob_dict[s + t]
            yield p

    return sum(sum_iter())


def pos_dep_weight(wt_pos_df, pos_dep_prob_dict, pos_dep_sum, source_position, mut_aa):
    def pos_dep_codon_weight(source_position, mut_codon):
        """Compute the position-dependent mutation weight"""

        def codon_weight(current: str, remain_permutation_number):
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

            p = pos_dep_prob_dict[keyword] / pos_dep_sum
            new_current = current[:head] + mut_codon[head - 1] + current[head + 1 :]
            result = p * codon_weight(new_current, remain_permutation_number[1:])
            return 0 if result == 1 else result

        pre = codon_from_aa_pos(wt_pos_df, source_position - 1)[2]
        codon = codon_from_aa_pos(wt_pos_df, source_position)
        end = codon_from_aa_pos(wt_pos_df, source_position + 1)[0]
        return sum(
            codon_weight(pre + codon + end, _)
            for _ in permutations(
                _ for _ in range(1, 4) if codon[_ - 1] != mut_codon[_ - 1]
            )
        )

    return sum(
        map(lambda x: pos_dep_codon_weight(source_position, x), to_codons(mut_aa))
    )


def clean_cscs(cscs: pd.DataFrame):
    cscs_copy = cscs.copy()
    cscs_copy["pos"] += 1  # Increase position to be started from 1, not 0.

    # get rid of row that matches the condition "pos == 1" (starting position)
    cscs_copy = cscs_copy[cscs_copy["pos"] != 1]
    # get rid of unusable alphabet
    cscs_copy = cscs_copy[
        ~cscs_copy["wt"].str.contains("X|B|Z|J|U", case=False)
        & ~cscs_copy["mut"].str.contains("X|B|Z|J|U", case=False)
    ]
    return cscs_copy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rank escape mutants by combining ranks which are based on grammar, \
            semantic change and codon mutant score."
    )
    parser.add_argument("cscs", type=str, help="path of cscs table")
    parser.add_argument("wt", type=str, help="path of wild type position table(csv)")
    parser.add_argument(
        "--indep",
        type=str,
        help="path of position independent codon mutation probability(json)",
    )
    parser.add_argument(
        "--dep",
        type=str,
        help="path of position dependent codon mutation probability table(json)",
    )
    args = parser.parse_args()
    cscs = clean_cscs(pd.read_csv(args.cscs, sep="\t"))
    wt_pos = pd.read_csv(args.wt)

    # get codon mutation value and rank respectively
    if args.indep:
        with open(args.indep) as json_fs:
            indep = json.load(json_fs)
            cscs["codon indep"] = cscs.apply(
                lambda row: pos_indep_prob(wt_pos, indep, row["pos"], row["mut"]),
                axis=1,
            )
    if args.dep:
        with open(args.dep) as json_fs:
            dep = json.load(json_fs)
            dep_sum = sum(dep.values())
            cscs["codon dep"] = cscs.apply(
                lambda row: pos_dep_weight(
                    wt_pos, dep, dep_sum, row["pos"], row["mut"]
                ),
                axis=1,
            )

    file_name = "extended_" + os.path.basename(args.cscs)
    cscs.to_csv(file_name, index=False, sep="\t")
    print(file_name + " has been written")
