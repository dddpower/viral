# get the statistics of given protein sequences
# cli

from Bio import SeqIO
from collections import Counter
import sys
import matplotlib.pyplot as plt

filename = sys.argv[1]
XposY = sys.argv[2]


def decode(XposY):
    wildtype_aa = XposY[0]
    position = int(XposY[1:-1]) + 1
    mutation_aa = XposY[-1]
    return wildtype_aa, position, mutation_aa


def isPositionY(seq, position, y):
    return seq[position] == y


def whatIsPositionY(seq, position):
    return seq[position]


def show_statistics(AAs):
    def draw_histogram(ratios):
        # 히스토그램 그리기
        alphabets, ratios = zip(*ratios)
        plt.bar(alphabets, ratios)
        plt.xlabel("Alphabet")
        plt.ylabel("Ratio")
        plt.title("Alphabet Histogram")
        plt.show()

    counter = Counter(AAs)
    total_number = sum(counter.values())
    alphabet_ratios = []
    for char, count in sorted(counter.items()):
        ratio = count / total_number
        alphabet_ratios.append((char, ratio))
        print(f"{char}: {ratio:.2f}")

    draw_histogram(alphabet_ratios)


if __name__ == "__main__":
    parses = SeqIO.parse(filename, "fasta")
    seqs = [record.seq for record in parses]
    wildtype_aa, position, mutatino_aa = decode(XposY)
    nth_positional_chars = map(lambda x: whatIsPositionY(x, position), seqs)
    show_statistics(nth_positional_chars)
