sys.stderr = open(snakemake.log[0], "w")
# parameter = snakemake.params.get("parameter", "")

from collections import defaultdict
from os.path import basename, splitext

import pandas as pd


def summarize_manuel_checks(
    paths_to_manuell_check_files: list[str], path_to_total_summary: str, out_path: str
):
    summary_dict = defaultdict()

    summary_dict["total pages processed"] = pd.read_csv(
        path_to_total_summary, sep="\t"
    ).shape[0]

    for path in paths_to_manuell_check_files:
        header = splitext(basename(path))[0].replace("_", " ").replace("-", " ")
        count = pd.read_csv(path, sep="\t", names=[header]).shape[0]
        summary_dict[header] = count

    manuel_check_summary_df = pd.DataFrame(
        summary_dict.items(), columns=["Check", "Count"]
    )
    manuel_check_summary_df.to_csv(out_path, sep="\t", index=False)


if __name__ == "__main__":
    summarize_manuel_checks(
        snakemake.input.manuel_checks,
        snakemake.input.total_imgs_processed,
        snakemake.output[0],
    )
