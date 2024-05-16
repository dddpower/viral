import argparse
import pandas as pd
import functools


def nofify_file(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        print(" ".join((kwargs.get("file_name"), "has been written.")))
        return result

    return wrapper


if __name__ == "__main__":
    """Shrink the original position "12276_2021_658_MOESM2_ESM.xlsx" to covid, flu, hiv
    since the file size is too big."""
    parser = argparse.ArgumentParser(
        description="shrink the size of the given position dependent mutation table"
    )
    parser.add_argument(
        "path", type=str, help="path of position dependent mutation table(excel)"
    )
    args = parser.parse_args()
    pos_dep_probs = pd.read_excel(args.path, sheet_name=3)

    @nofify_file
    def df_to_json_helper(virus_col, file_name):
        pd.Series(
            pos_dep_probs[virus_col].values, index=pos_dep_probs["Substitution type"]
        ).to_json(file_name)

    prefix = "dep_table_"
    postifx = ".json"
    for args in [("SARS-CoV-2", "covid"), ("Influenza A", "flu"), ("HIV", "hiv")]:
        df_to_json_helper(args[0], file_name=prefix + args[1] + postifx)
