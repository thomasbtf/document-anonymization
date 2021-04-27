sys.stderr = open(snakemake.log[0], "w")
#parameter = snakemake.params.get("parameter", "")

import altair as alt
import pandas as pd


def plot_manuel_check_summary(path_to_manuel_check_summary: str, out_path: str):
    source = pd.read_csv(path_to_manuel_check_summary, sep="\t")

    bars = alt.Chart(source).mark_bar().encode(
        x="Count:Q",
        y="Check:O",
    )

    text = bars.mark_text(
        align='left',
        baseline='middle',
        dx=3  # Nudges text to right so it doesn't appear on top of the bar
    ).encode(
        text='Count:Q'
    )

    (bars + text).save(out_path)


if __name__ == "__main__":
    plot_manuel_check_summary(snakemake.input[0], snakemake.output[0])
