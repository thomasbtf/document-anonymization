sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict
import pandas as pd


def summarize_found_personal_data(
    data_path_list: list[str], img_path_list: list[str], sm_output: str, max_dist: int
):
    """Summarizes the found personal data. Saves the summary as a tsv-file.

    Args:
        data_path_list (list[str]): Paths of identify personal data.
        img_path_list (list[str]): Paths to redacted images.
        sm_output (str): Path to write summary to.
        max_dist (int): Mixmal Levenshtein distance.
    """
    summary_list = []
    for data_path, img_path in zip(data_path_list, img_path_list):

        found_data_df = pd.read_csv(data_path, sep="\t")
        page_summary = defaultdict()

        page_summary["processed img"] = img_path
        page_summary["# personal data"] = found_data_df.shape[0]

        tesseract_output = {"left", "top", "width", "height", "conf", "text"}
        personal_data_columns = set(found_data_df.columns) - tesseract_output
        for column in personal_data_columns:
            no_found_data = found_data_df[found_data_df[column] <= max_dist][
                column
            ].shape[0]
            if no_found_data > 0:
                page_summary[column] = no_found_data

        summary_list.append(page_summary)

    pd.DataFrame(summary_list).to_csv(sm_output, index=False, sep="\t")


if __name__ == "__main__":
    summarize_found_personal_data(
        data_path_list=snakemake.input.data,
        img_path_list=snakemake.input.pages,
        sm_output=snakemake.output[0],
        max_dist=snakemake.config["max-distance"],
    )
