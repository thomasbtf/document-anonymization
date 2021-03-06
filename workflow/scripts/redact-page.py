import typing

import cv2
import pandas as pd


def process_page(image_path: str, out_path: str, data_to_redact: str, version: str):
    """Analyzes the passed image and removes personal information.

    Args:
        image_path (str): path to the image
        out_path (str): path where the redacted image should be written to
        version (str): version number of the workflow
    """

    df = pd.read_csv(data_to_redact, sep="\t")
    img = cv2.imread(image_path)

    img = add_watermark(img, version)
    img = redact(df, img)

    if not ".jpg" in out_path[-3:]:
        "".join([out_path, ".jpg"])

    cv2.imwrite(out_path, img)


def add_watermark(img: typing.Any, version: str) -> typing.Any:
    x, y = 50, 50
    watermark_text = "anonymized by DocNo {}".format(version)
    cv2.putText(img, watermark_text, (x, y), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 0), 2)
    return img


def redact(personal_data_df: pd.DataFrame, img: typing.Any) -> typing.Any:
    """Redacts personal data.

    Args:
        personal_data_df (pd.DataFrame): personal data with location on image.
        img (typing.Any): image with personal data on it.

    Returns:
        typing.Any: redacted image.
    """

    for i in personal_data_df.index:
        (x, y, w, h) = (
            int(personal_data_df.loc[i].left),
            int(personal_data_df.loc[i].top),
            int(personal_data_df.loc[i].width),
            int(personal_data_df.loc[i].height),
        )
        img = cv2.rectangle(img, (x, y), (x + w, y + h), (0, 0, 0), -1)
    return img


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    version = snakemake.params.get("version", "")
    process_page(
        image_path=snakemake.input.orginal_page,
        out_path=snakemake.output[0],
        data_to_redact=snakemake.input.data_to_redact,
        version=version,
    )
