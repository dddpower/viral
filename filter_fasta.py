import sys
from Bio import SeqIO
path = sys.argv[1]
fasta = SeqIO.parse(open(path, encoding='latin-1'), 'fasta')
i = 0
for fa in fasta:
    nation = fa.id.split("|")[1].split("/")[1]
    if "Korea" in nation or "China" in nation:
        i += 1

print(i)
