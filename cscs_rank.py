import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description= "rank grammar & semantic change score of cscs table")
parser.add_argument("--cscs-paths", nargs="+", type=str, help="path of cscs table")


if __name__ == "__main__":
    args = parser.parse_args()
    paths = args.cscs_paths
    for path in paths:
        file = os.path.basename(path)
        cscs = pd.read_csv(path, sep='\t')
        gram_rank = "grammar rank"
        sem_rank = "semantic change rank"
        cscs[[gram_rank, sem_rank]] = cscs[['prob', 'change']].rank(method='min', ascending=False)
        cscs['+'.join([gram_rank, sem_rank]) + " rank"] = (cscs[gram_rank] + cscs[sem_rank]).rank(method='min', ascending=False)
        new_file = "ranked_" + file
        cscs.to_csv(new_file, sep='\t', index=False)
        print(new_file + " created")