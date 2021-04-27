import typing
import json

import cv2
import Levenshtein
import pandas as pd
import pytesseract
from pytesseract import Output


def parse_page(
    image_path: str,
    out_path_all_text: str,
    out_path_personal_data: str,
    personal_data_path: dict,
    min_conf: float = 0.6,
    max_dist: int = 2,
):
    """Analyzes the passed image and identifies personal information on it.

    Args:
        image_path (str): path to the image
        out_path_all_text (str): path where all text should be written to
        out_path_personal_data (str): path where personal data should be written to
        personal_data (dict): path to personal data that should be made unrecognizable
        min_conf (float, optional): minimal OCR confidence score. Defaults to 0.6.
        max_dist (int, optional): maximum Levenshtein distance of the found text on the image to the personal data. Defaults to 2.
    """

    img = cv2.imread(image_path)

    with open(personal_data_path) as json_file:
        personal_data = json.load(json_file)

    df = detect_text(img, min_conf)
    df.to_csv(out_path_all_text, index=False, sep="\t")

    df = select_personal_data(df, personal_data, max_dist)
    df.to_csv(out_path_personal_data, index=False, sep="\t")


def detect_text(img: typing.Any, min_conf: float) -> pd.DataFrame:
    """Recognizes text on the image.

    Args:
        img (typing.Any): image with text to be recognized.
        min_conf (float): minimum OCR Confidence Scores.

    Returns:
        pd.DataFrame: all found text on image with text field data, filtered by min_conf.
    """

    # ocr
    detected_text_df = pytesseract.image_to_data(img, output_type=Output.DATAFRAME)

    # filter ocr table
    detected_text_df = detected_text_df[detected_text_df.conf >= min_conf]
    detected_text_df.drop(
        columns=["level", "page_num", "block_num", "par_num", "line_num", "word_num"],
        inplace=True,
    )

    detected_text_df.text = detected_text_df.text.astype(str)
    detected_text_df.text = detected_text_df.text.str.lower()

    return detected_text_df


def select_personal_data(
    detected_text_df: pd.DataFrame, personal_data: dict, max_dist: int
) -> pd.DataFrame:
    """Identifies personal data from the detected text.

    Args:
        detected_text_df (pd.DataFrame): detected text on the image.
        personal_data (dict): personal data to be masked out
        max_dist (int): maximum Levenshtein distance of the found text on the image to the personal data.

    Returns:
        pd.DataFrame: person data with location on image, filtered by max_dist.
    """

    final_df = pd.DataFrame(columns=["left", "top", "width", "height", "conf", "text"])

    max_spaces = max([value.count(" ") for value in personal_data.values()])
    for no_spaces in range(max_spaces + 1):
        tmp_df = detected_text_df.copy()
        tmp_df.rename(columns={"text": "text_0", "width": "width_0"}, inplace=True)

        # subset of personal data dict, according to # spaces
        tmp_dict = {
            key: value
            for key, value in personal_data.items()
            if value.count(" ") == no_spaces
        }

        # shift text to get longer phrases
        if no_spaces > 0:
            shift_colums = []
            for shift in range(1, no_spaces + 1):
                # shift text column and aggreagte
                shift_colums.append("text_" + str(shift))
                highest_text_column_name = "text_" + str(shift)
                tmp_df[highest_text_column_name] = (
                    tmp_df["text_" + str(shift - 1)]
                    + " "
                    + tmp_df.text_0.shift(-shift).fillna("")
                )

                # shift width column and aggreagte
                shift_colums.append("width_" + str(shift))
                highest_width_column_name = "width_" + str(shift)
                tmp_df[highest_width_column_name] = (
                    tmp_df["width_" + str(shift - 1)]
                    + tmp_df.width_0.shift(-shift).fillna(0)
                    + tmp_df.width_0 / tmp_df.text_0.str.len()
                )

            tmp_df["text_0"] = tmp_df[highest_text_column_name]
            tmp_df["width_0"] = tmp_df[highest_width_column_name].astype(int)
            tmp_df.drop(columns=shift_colums, inplace=True)

        tmp_df.rename(columns={"text_0": "text", "width_0": "width"}, inplace=True)

        # calc edit distances for each key in subsampled dict
        for key, value in tmp_dict.items():
            tmp_df[key] = tmp_df["text"].apply(
                lambda text: Levenshtein.distance(value, text)
            )

        # select entries, where distance is below or equal a threshhold
        query_list = []
        for col in tmp_dict.keys():
            query_list.append(f"{col}<={max_dist}")

        # if query list is not empty, thus personal data was found,
        # then append filtered df
        if query_list:
            tmp_df = tmp_df.query(" | ".join(query_list))
            final_df = final_df.append(tmp_df, ignore_index=True)

    return final_df


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    parse_page(
        image_path=snakemake.input.preprocessed_page,
        out_path_all_text=snakemake.output.all_text,
        out_path_personal_data=snakemake.output.text_to_redact,
        personal_data_path=snakemake.input.personal_data,
        min_conf=snakemake.config["min-confidence"],
        max_dist=snakemake.config["max-distance"],
    )
