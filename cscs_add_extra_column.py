from Bio.Data import CodonTable
import pandas as pd
from itertools import permutations, combinations
import json
from collections import namedtuple
import argparse

VirusData = namedtuple(
    "VirusData",
    "name wt_positions pos_indep_prob_df, pos_dep_prob_dict pos_dep_sum cscs",
    defaults=["", None, None, None, None, None],
)
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


"""Add position dependent/independent codon mutation chance,
value columns to the original cscs result table. For this purpose, this class stores
informations - cscs result table, position dependent/independent mutation table respectively,
wildtype codon position table for each virus.
"""


def codon_from_aa_pos(wt_positions, source_position):
    return wt_positions[wt_positions["pos"] == source_position]["Codon"].values[0]


def pos_indep_prob(pos_indep_prob_df, aa_pos, target: str):
    """wildtype codon to mutant amino_acid independent probability"""
    source_codon = codon_from_aa_pos(aa_pos)
    target_codons = filter(lambda x: x != source_codon, to_codons(target))

    def sum_iter():
        for target_codon in target_codons:
            p = 1
            for s, t in zip(source_codon, target_codon):
                p *= pos_indep_prob_df[s + t]
            yield p

    return sum(sum_iter())


def pos_dep_weight(pos_dep_prob_dict, pos_dep_sum, source_position, mut_aa):
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


def codon_col(cscs: pd.DataFrame, func):
    cscs["pos"] += 1
    """remove unusable elements"""
    # 1. get rid of row that matches condition "pos == 1"
    cscs = cscs[cscs["pos"] != 1]
    # 2. get rid of unusable alphabet
    cscs = cscs[
        ~cscs["wt"].str.contains("X|B|Z|J|U", case=False)
        & ~cscs["mut"].str.contains("X|B|Z|J|U", case=False)
    ]
    return cscs.apply(lambda row: func(row["pos"], row["mut"]), axis=1)

    # if virus.pos_dep_prob_dict:
    #     cscs["pos dep"] = cscs.apply(
    #         lambda row: pos_dep_weight(row["pos"], row["mut"]), axis=1
    #     )

    # if virus.pos_indep_prob_df:
    #     cscs["pos indep"] = cscs.apply(
    #         lambda row: pos_indep_prob(row["pos"], row["mut"]), axis=1
    #     )

    # make rank columns


def rank(cscs: pd.DataFrame, col_name):
    return cscs[col_name].rank(method="min", ascending=False)

    cscs["grammar rank"] = cscs["prob"].rank(method="min", ascending=False)
    cscs["semantic rank"] = cscs["change"].rank(method="min", ascending=False)

    if virus.pos_indep_prob_df:
        cscs["pos_indep_rank"] = cscs["pos_indep_prob"].rank(
            method="min", ascending=False
        )

    if virus.pos_dep_prob_dict:
        cscs["pos_dep_rank"] = cscs["pos_dep_weight"].rank(
            method="min", ascending=False
        )

    base_list = ["grammar rank", "semantic rank"]

    if virus.pos_indep_prob_df:
        for x, y in combinations(base_list + ["pos_indep_rank"], 2):
            cscs[x + "+" + y] = cscs[x] + cscs[y]

    if virus.pos_dep_prob_dict:
        for x, y in combinations(base_list + ["pos_dep_rank"], 2):
            cscs[x + "+" + y] = cscs[x] + cscs[y]

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

    target = [_ for _ in col_added_df.columns if "rank" in _]
    # target = ["grammar rank", "semantic rank", "pos_indep_rank", "pos_dep_rank"]
    for x, y in combinations(target, 2):
        df["+".join((x, y)) + " rank"] = (df[x] + df[y]).rank(
            method="min", ascending=True
        )

    assert len(target) >= 3
    for x, y, z in combinations(target, 3):
        df["+".join((x, y, z)) + " rank"] = (df[x] + df[y] + df[z]).rank(
            method="min", ascending=True
        )
    df = df[df["is_escape"]]
    return map(mean_std_helper, filter(lambda string: "rank" in string, df.columns))


# Load position-dependent probability table, transform it to dictionary
input_file_name = "12276_2021_658_MOESM2_ESM.xlsx"
pos_dep_probs = pd.read_excel(input_file_name, sheet_name=3)
cov_pos_dep_probs, cov_pos_dep_sum = pos_dep_probs_helper("SARS-CoV-2", pos_dep_probs)
flu_pos_dep_probs, flu_pos_dep_sum = pos_dep_probs_helper("Influenza A", pos_dep_probs)
hiv_pos_dep_probs, hiv_pos_dep_sum = pos_dep_probs_helper("HIV", pos_dep_probs)
with open("cov_static_probability.json", "r") as json_fs:
    cov_pos_indep_probs = json.load(json_fs)

with open("flu_static_probability.json", "r") as json_fs:
    flu_pos_indep_probs = json.load(json_fs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rank escape mutants by combining rank based on grammar, \
            semantic change and codon mutant score."
    )
    parser.add_argument("cscs", type=str, help="path of cscs table")
    parser.add_argument("wt", type=str, help="path of wild type postion table(csv)")
    parser.add_argument(
        "--indep",
        type=str,
        help="path of position independent codon mutation probability(json)",
    )
    parser.add_argument(
        "--dep",
        type=str,
        help="path of position dependent codon mutation probability table(csv)",
    )
    args = parser.parse_args()
    cscs = pd.read_csv(args.cscs, sep="\t")
    wt_pos = pd.read_csv(args.wt)

    # make rank of "prob", "change"
    cscs["grammar rank"] = rank(cscs, "prob")
    cscs["semantic change rank"] = rank(cscs, "change")
    # get codon mutation value and rank respectively
    if args.indep:
        with open(args.indep, "r") as json_fs:
            indep = json.load(json_fs)
        cscs["indep"] = codon_col(cscs, pos_indep_prob)
        cscs["indep rank"] = rank(cscs, "indep")
    if args.dep:
        with open(args.dep, "r") as json_fs:
            dep = json.load(json_fs)
        cscs["dep"] = codon_col(cscs, pos_dep_weight)
        cscs["dep rank"] = rank(cscs, "dep")

    # add codon mutation colum

    # covid = Virus(
    #     wt_positions=pd.read_csv("cov_wildtype_codon.csv"),
    #     pos_indep_probs=cov_pos_indep_probs,
    #     pos_dep_probs=cov_pos_dep_probs,
    #     pos_dep_sum=cov_pos_dep_sum,
    #     cscs=pd.read_csv(
    #         "results/cov/semantics/analyze_semantics_cov_bilstm_512.txt", sep="\t"
    #     ),
    # )
    # def analyze_and_create(prefix, virus)
    # # covid
    # cov_df = CSCSColumnAdder(*covid).make_extra_column()
    # cov_df.to_csv(today + "covid.csv")
    # cov_df[cov_df["is_escape"]].to_csv(today + "covid_escape.csv")
    # with open(today + "cov_rank.txt", "w") as f:
    #     f.write("\n".join(calc_combined_rank(cov_df)))

    # influenza = Virus(
    #     wt_positions=pd.read_csv("flu_wildtype_codon.csv"),
    #     pos_indep_probs=flu_pos_indep_probs,
    #     pos_dep_probs=flu_pos_dep_probs,
    #     pos_dep_sum=flu_pos_dep_sum,
    #     cscs=pd.read_csv(
    #         "results/flu/semantics/analyze_semantics_flu_h1_bilstm_512.txt", sep="\t"
    #     ),
    # )

    # # influenza
    # flu_df = CSCSColumnAdder(*influenza).make_extra_column()
    # flu_df.to_csv(today + "flu.csv")
    # flu_df[flu_df["is_escape"]].to_csv(today + "flu_escape.csv")
    # with open(today + "flu_rank.txt", "w") as f:
    #     f.write("\n".join(calc_combined_rank(flu_df)))

    # hiv = Virus(
    #     wt_positions=pd.read_csv("hiv_wildtype_codon.csv"),
    #     pos_indep_probs=None,
    #     pos_dep_probs=hiv_pos_dep_probs,
    #     pos_dep_sum=hiv_pos_dep_sum,
    #     cscs=pd.read_csv(
    #         "results/hiv/semantics/analyze_semantics_hiv_bilstm_512.txt", sep="\t"
    #     ),
    # )
    # # hiv
    # hiv_df = CSCSColumnAdder(*hiv).make_extra_column()
    # hiv_df.to_csv(today + "hiv.csv")
    # hiv_df[hiv_df["is_escape"]].to_csv(today + "hiv_escape.csv")
    # with open(today + "hiv_rank.txt", "w") as f:
    #     f.write("\n".join(calc_combined_rank(hiv_df)))
