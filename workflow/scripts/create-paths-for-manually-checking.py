sys.stderr = open(snakemake.log[0], "w")

from os import sep
import pandas as pd


def save_df(df: pd.DataFrame, out_path: str):
    df.to_csv(out_path, sep="\t", index=False, header=False)


def no_redaction(summary_df: pd.DataFrame, out_path: str):
    save_df(summary_df[summary_df["# personal data"] == 0][["processed img"]], out_path)


def high_degree_of_redaction(summary_df: pd.DataFrame, out_path: str):
    save_df(
        summary_df[summary_df["# personal data"] >= 10][["processed img"]], out_path
    )


def partly_found_address(summary_df: pd.DataFrame, out_path: str):
    if "city" in summary_df.columns and "address" in summary_df.columns:
        df = summary_df[summary_df["city"] != summary_df["address"]][["processed img"]]
    else:
        df = pd.DataFrame(columns=["processed img"])
    
    save_df(df, out_path)


def partly_found_name(summary_df: pd.DataFrame, out_path: str):
    if "name_family" in summary_df.columns and "name_first_0" in summary_df.columns:
        df = summary_df[summary_df["name_family"] != summary_df["name_first_0"]][["processed img"]]
    else:
        df = pd.DataFrame(columns=["processed img"])
    
    save_df(df, out_path)


if __name__ == "__main__":
    summary_df = pd.read_csv(snakemake.input[0], sep="\t")
    summary_df.fillna(999999999.0, inplace=True)
    no_redaction(summary_df, snakemake.output.no_redaction)
    high_degree_of_redaction(summary_df, snakemake.output.high_degree_of_redaction)
    partly_found_address(summary_df, snakemake.output.partly_found_address)
    partly_found_name(summary_df, snakemake.output.partly_found_name)
