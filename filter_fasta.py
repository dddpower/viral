from Bio import SeqIO
fasta_seqs = SeqIO.parse(open("unique.fasta"), 'fasta')
i = 0
for fa in fasta_seqs:
    if i == 10:
        break
    i += 1
    # type(fa.id) is str
    id_list = fa.id.split("|")
    fa.id = id_list[5]
with open("test.fasta", "w") as out_file:
    out_file.write(fasta_seqs)
i1 = 0
i2 = 0
i3 = 0
# with open("Korea.fasta", "w") as out_file:
#     try:
#         for fasta in fasta_seqs:
#             i1 += 1
#             name, seq = fasta.id, str(fasta.seq)
#             if "Korea" in name:
#                 i3 += 1
#                 if seq[-1] == '*':
#                     seq = seq[:-1]
#                 out_file.write('>' + name + '\n')
#                 out_file.write(seq + '\n')
#     except UnicodeDecodeError as ue:
#         print("in ", i1, "th iteration, exception occurred")
#         i2 += 1
#         pass
# 
# out_file.close()
print(i1)
print(i2)
print(i3)
print("finished")
