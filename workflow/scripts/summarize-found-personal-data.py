sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict
import pandas as pd


def summarize_found_personal_data(data_path_list: list[str], img_path_list: list[str], sm_output: str):
    summary_list=[]
    for data_path, img_path in zip(data_path_list, img_path_list):
        page_summary = defaultdict()
        found_data_df = pd.read_csv(data_path, sep="\t")

        page_summary["processed img"] = img_path
        page_summary["# personal data"] = found_data_df.shape[0]

        summary_list.append(page_summary)

    pd.DataFrame(summary_list).to_csv(sm_output, index=False, sep="\t")


if __name__ == "__main__":
    summarize_found_personal_data(snakemake.input.data, snakemake.input.pages, snakemake.output[0])
