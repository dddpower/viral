import argparse
import pandas as pd
import os
from itertools import combinations

if __name__ == "__main__":
    """Rank the given scores(basically 'prob', 'change', 'codon_dep', 'codon_indep') & rank the 
    sum of each ranks.
    """
    parser = argparse.ArgumentParser(description= "rank the given cscs table")
    parser.add_argument("path", type=str, help="path of cscs table")
    parser.add_argument('--dep', action='store_true', help='position dependent codon mut score is used')
    parser.add_argument('--indep', action='store_true', help='position independent codon mut score is used')
    args = parser.parse_args()
    path = args.path
    file = os.path.basename(path)
    cscs = pd.read_csv(path, sep='\t')
    gram_rank = "grammar rank"
    sem_rank = "semantic change rank"
    print(args)
    mut_type = "dep" if args.dep is True else 'indep'
    codon_rank =mut_type + " rank"
    rank_kargs = {"method":'min', "ascending":True}
    cscs[[gram_rank, sem_rank, codon_rank]] = cscs[['prob', 'change', 'codon ' + mut_type]].rank(**rank_kargs)
    for combs in combinations((gram_rank, sem_rank, codon_rank), 2):
        cscs['+'.join(combs) + " rank"] = (cscs[combs[0]] + cscs[combs[1]]).rank(**rank_kargs)
    cscs["+".join((gram_rank, sem_rank, codon_rank)) + " rank"] = (cscs[gram_rank] + cscs[sem_rank] + cscs[codon_rank]).rank(**rank_kargs)

    new_file = "ranked_" + file
    cscs.to_csv(new_file, sep='\t', index=False)
    print(new_file + " has been written")
