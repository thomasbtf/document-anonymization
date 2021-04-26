sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict
import pandas as pd


def summarize_found_personal_data(sm_input: list[str], sm_output: str):
    summary_list=[]
    for path in sm_input:
        page_summary = defaultdict()
        found_data_df = pd.read_csv(path, sep="\t")

        page_summary["# personal data"] = found_data_df.shape[0]

        summary_list.append(page_summary)

    pd.DataFrame(summary_list).to_csv(sm_output, index=False, sep="\t")


if __name__ == "__main__":
    summarize_found_personal_data(snakemake.input, snakemake.output[0])
