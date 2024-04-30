import argparse
from Bio import SeqIO
import matplotlib.pyplot as plt
from collections import Counter


# alphabet counting for each position for every sequence
if __name__ == "__main__":
    FIXED_LENGTH = 1273
    not_fixed_cases = []
    cov_path = "data/cov/cov_all.fa"
    flu_path = "data/influenza/"
    cov_records = SeqIO.parse(open("data/cov/cov_all.fa"), "fasta")
    flu_records = SeqIO.parse(open("data/influenza/ird_influenzaA_HA_allspecies.fa"), "fasta")
    hiv_records = SeqIO.parse(open("data/hiv/HIV-1_env_samelen.fa"), "fasta")

    def process(val):
        return Counter(list(map(lambda x: len(x.seq), val)))

    cov_count_result = process(cov_records)
    flu_count_result = process(flu_records)
    hiv_count_result = process(hiv_records)
    