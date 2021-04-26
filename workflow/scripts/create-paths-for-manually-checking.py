sys.stderr = open(snakemake.log[0], "w")

from os import sep
import pandas as pd


def save_df(df: pd.DataFrame, out_path: str):
    df.to_csv(out_path, sep="\t", index=False, header=False)


def no_personal_data_found(summary_df: pd.DataFrame, out_path: str):
    save_df(summary_df[summary_df["# personal data"] == 0][["processed img"]], out_path)


def lots_personal_data_found(summary_df: pd.DataFrame, out_path: str):
    save_df(
        summary_df[summary_df["# personal data"] >= 10][["processed img"]], out_path
    )


def address_not_entirely_found(summary_df: pd.DataFrame, out_path: str):
    save_df(
        summary_df[summary_df["city"] != summary_df["city"]][["processed img"]],
        out_path,
    )


def name_not_entirely_found(summary_df: pd.DataFrame, out_path: str):
    save_df(
        summary_df[summary_df["name_family"] != summary_df["name_first_0"]][
            ["processed img"]
        ],
        out_path,
    )


if __name__ == "__main__":
    summary_df = pd.read_csv(snakemake.input[0], sep="\t")
    no_personal_data_found(summary_df, snakemake.output.no_personal_data_found)
    lots_personal_data_found(summary_df, snakemake.output.lots_personal_data_found)
    address_not_entirely_found(summary_df, snakemake.output.address_not_entirely_found)
    name_not_entirely_found(summary_df, snakemake.output.name_not_entirely_found)
