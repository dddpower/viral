from Bio.Data import CodonTable
import pandas as pd
from itertools import permutations, combinations
import json
import argparse
from enum import Enum, auto
from collections import namedtuple

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


class CSCSColumnAdder:
    def __init__(
        self,
        wt_positions,
        pos_indep_prob,
        pos_dep_prob,
        pos_dep_sum,
        cscs,
    ) -> None:
        self.wt_positions = wt_positions
        self.pos_indep_prob_df = pos_indep_prob
        self.pos_dep_prob_dict = pos_dep_prob
        self.pos_dep_sum = pos_dep_sum
        self.cscs = cscs

    def codon_from_aa_pos(self, source_position):
        wt_positions = self.wt_positions
        return wt_positions[wt_positions["pos"] == source_position]["Codon"].values[0]

    def pos_indep_prob(self, aa_pos, target: str):
        """wildtype codon to mutant amino_acid independent probability"""
        source_codon = self.codon_from_aa_pos(aa_pos)
        target_codons = filter(lambda x: x != source_codon, to_codons(target))

        def sum_iter():
            for target_codon in target_codons:
                p = 1
                for s, t in zip(source_codon, target_codon):
                    p *= self.pos_indep_prob_df[s + t]
                yield p

        return sum(sum_iter())

    def pos_dep_weight(self, source_position, mut_aa):
        # must ignore first 'M'

        def pos_dep_codon_weight(source_position, mut_codon):
            codon_from_aa_pos = self.codon_from_aa_pos
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

                p = self.pos_dep_prob_dict[keyword] / self.pos_dep_sum
                new_current = current[:head] + mut_codon[head - 1] + current[head + 1 :]
                result = p * codon_weight(new_current, remain_permutation_number[1:])
                return 0 if result == 1 else result

            pre = codon_from_aa_pos(source_position - 1)[2]
            codon = codon_from_aa_pos(source_position)
            end = codon_from_aa_pos(source_position + 1)[0]
            return sum(
                codon_weight(pre + codon + end, _)
                for _ in permutations(
                    _ for _ in range(1, 4) if codon[_ - 1] != mut_codon[_ - 1]
                )
            )

        return sum(
            map(lambda x: pos_dep_codon_weight(source_position, x), to_codons(mut_aa))
        )

    def make_extra_column(self):
        cscs = self.cscs
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
            lambda row: self.pos_dep_weight(row["pos"], row["mut"]), axis=1
        )
        cscs["pos_indep_prob"] = cscs.apply(
            lambda row: self.pos_indep_prob(row["pos"], row["mut"]), axis=1
        )

        # make rank columns
        cscs["grammar rank"] = cscs["prob"].rank(method="min", ascending=False)
        cscs["semantic rank"] = cscs["change"].rank(method="min", ascending=False)
        cscs["pos_indep_rank"] = cscs["pos_indep_prob"].rank(
            method="min", ascending=False
        )
        cscs["pos_dep_rank"] = cscs["pos_dep_weight"].rank(
            method="min", ascending=False
        )

        base_list = ["grammar rank", "semantic rank"]
        for x, y in combinations(base_list + ["pos_indep_rank"], 2):
            cscs[x + "+" + y] = cscs[x] + cscs[y]
        for x, y in combinations(base_list + ["pos_dep_rank"], 2):
            cscs[x + "+" + y] = cscs[x] + cscs[y]

        # name1 = "grammar rank"
        # name2 = "semantic rank"
        # name3 = "pos_indep_rank"
        # name4 = "pos_dep_rank"
        # sum_rank0 = "+".join([name1, name2])
        # sum_rank1 = "+".join([name1, name2, name3])
        # sum_rank2 = "+".join([name1, name2, name4])
        # cscs[sum_rank0] = cscs[name1] + cscs[name2]
        # cscs[sum_rank1] = cscs[name1] + cscs[name2] + cscs[name3]
        # cscs[sum_rank2] = cscs[name1] + cscs[name2] + cscs[name4]
        # cscs[sum_rank0 + " rank"] = cscs[sum_rank0].rank(method="min", ascending=True)
        # cscs[sum_rank1 + " rank"] = cscs[sum_rank1].rank(method="min", ascending=True)
        # cscs[sum_rank2 + " rank"] = cscs[sum_rank2].rank(method="min", ascending=True)
        return cscs


def pos_dep_probs_helper(virus_name, pos_dep_probs):
    return (
        pd.Series(
            pos_dep_probs[virus_name].values, index=pos_dep_probs["Substitution type"]
        ).to_dict(),
        pos_dep_probs[virus_name].sum(),
    )


def calc_combined_rank(col_added_df):
    df = col_added_df

    def mean_std_helper(name):
        return f"mean and std of {name} is {df[name].mean()}, {df[name].std()}"

    for x, y in combinations(
        ["grammar rank", "semantic rank", "pos_indep_rank", "pos_dep_rank"], 2
    ):
        df["+".join((x, y)) + " rank"] = (df[x] + df[y]).rank(
            method="min", ascending=True
        )

    for x, y, z in combinations(
        ["grammar rank", "semantic rank", "pos_indep_rank", "pos_dep_rank"], 3
    ):
        df["+".join((x, y, z)) + " rank"] = (df[x] + df[y] + df[z]).rank(
            method="min", ascending=True
        )
    df = df[df["is_escape"] == True]
    return map(
        mean_std_helper, filter(lambda string: "rank rank" in string, df.columns)
    )


# Load position-dependent probability table, transform it to dictionary
input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
pos_dep_probs = pd.read_excel(input_file_name, sheet_name=3)
cov_pos_dep_probs, cov_pos_dep_sum = pos_dep_probs_helper("SARS-CoV-2", pos_dep_probs)
# print(cov_pos_dep_probs)
# print(pos_dep_probs["SARS-CoV-2"])
flu_pos_dep_probs, flu_pos_dep_sum = pos_dep_probs_helper("Influenza A", pos_dep_probs)
# print(flu_pos_dep_probs)
hiv_pos_dep_probs, hiv_pos_dep_sum = pos_dep_probs_helper("HIV", pos_dep_probs)


with open("cov_static_probability.json", "r") as json_fs:
    cov_pos_indep_probs = json.load(json_fs)

with open("flu_static_probability.json", "r") as json_fs:
    flu_pos_indep_probs = json.load(json_fs)

Virus = namedtuple(
    "Virus", "wt_positions pos_indep_probs, pos_dep_probs pos_dep_sum cscs"
)
covid = Virus(
    wt_positions=pd.read_csv("cov_wildtype_codon.csv"),
    pos_indep_probs=cov_pos_indep_probs,
    pos_dep_probs=cov_pos_dep_probs,
    pos_dep_sum=cov_pos_dep_sum,
    cscs=pd.read_csv(
        "results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", sep="\t"
    ),
)
influenza = Virus(
    wt_positions=pd.read_csv("flu_wildtype_codon.csv"),
    pos_indep_probs=flu_pos_indep_probs,
    pos_dep_probs=flu_pos_dep_probs,
    pos_dep_sum=flu_pos_dep_sum,
    cscs=pd.read_csv(
        "results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", sep="\t"
    ),
)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="provide file paths for column generation"
    )
    cov_df = CSCSColumnAdder(*covid).make_extra_column()
    cov_df.to_csv("0509_covid.csv")
    with open("cov_rank.txt", "w") as f:
        f.write("\n".join(calc_combined_rank(cov_df)))

    flu_df = CSCSColumnAdder(*influenza).make_extra_column()
    flu_df.to_csv("0509_flu.csv")
    with open("flu_rank.txt", "w") as f:
        f.write("\n".join(calc_combined_rank(flu_df)))
