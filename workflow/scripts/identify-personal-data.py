from os import replace
from re import split
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
    out_path_non_personal_data: str,
    personal_data_path: dict,
    replacements_path:str,
    min_conf: float = 0.6,
    max_dist: int = 2,
):
    """Analyzes the passed image and identifies personal information on it.

    Args:
        image_path (str): path to the image
        out_path_all_text (str): path where all text should be written to
        out_path_personal_data (str): path where personal data should be written to
        out_path_non_personal_data (str): path where non personal data should be written to
        personal_data (dict): path to personal data that should be made unrecognizable
        replacements: Path to replacement json
        min_conf (float, optional): minimal OCR confidence score. Defaults to 0.6.
        max_dist (int, optional): maximum Levenshtein distance of the found text on the image to the personal data. Defaults to 2.
    """

    img = cv2.imread(image_path)

    with open(personal_data_path) as json_file:
        personal_data = json.load(json_file)

    all_text = detect_text(img, min_conf)
    all_text.reset_index(inplace=True)
    all_text.to_csv(out_path_all_text, index=False, sep="\t")

    personal_text = select_personal_data(all_text, personal_data, max_dist)

    replace_and_save_personal_text(personal_text=personal_text, all_text=all_text, out_path_non_personal_data=out_path_non_personal_data, replacements_path=replacements_path, max_dist=max_dist)

    personal_text.drop(columns=["index"], inplace=True)
    personal_text.to_csv(out_path_personal_data, index=False, sep="\t")


def replace_and_save_personal_text(personal_text:pd.DataFrame, all_text:pd.DataFrame, out_path_non_personal_data:str, replacements_path:str, max_dist: int = 2):
    """Replaces and saves text on page.

    Args:
        personal_text (pd.DataFrame): DataFrame with all detected text
        all_text (pd.DataFrame): DataFrame with identifed personal
        out_path_non_personal_data (str): Path to write replaced text to.
    """
    personal_text = personal_text.copy()
    # remove personal data
    indices_to_remove = [
        ele.split(",") for ele in personal_text["index"].astype(str).values
    ]
    indices_to_remove = [
        int(float(item)) for sublist in indices_to_remove for item in sublist
    ]
    non_personal_text = all_text[~all_text["index"].isin(indices_to_remove)]

    # extract reason
    non_distance_columns=["index", "left", "top", "width", "height", "conf", "text"]
    distance_columns = list(set(personal_text.columns) - set(non_distance_columns))


    personal_text["reason"] = ""
    for col in distance_columns:
        personal_text[col] = personal_text[col].mask(personal_text[col] > float(max_dist))
        personal_text[col] = personal_text[col].mask(personal_text[col] <= float(max_dist), col)
        personal_text["reason"] = personal_text["reason"] +  personal_text[col].fillna("")

    personal_text = personal_text[["index", "reason"]]

    # insert replacements
    with open(replacements_path) as json_file:
        replacements = json.load(json_file)

    replaced_text = personal_text.copy().rename(columns={"reason" : "text"})
    for key in replacements.keys():
        replaced_text["text"][personal_text["reason"].str.contains(key)] = replacements[key]

    # if no replacement was found, replace identified personal data with "PrivateDataPrivateData"
    replaced_text["text"][replaced_text["text"] == personal_text["reason"]] == "PrivateDataPrivateData"

    # append replacements to whol text
    replaced_text['index'] = replaced_text['index'].astype(str)
    replaced_text['index'] = [x.split(',') for x in replaced_text['index']]
    replaced_text = replaced_text.explode("index")
    replaced_text['index'] = replaced_text['index'].astype(float).astype(int)
    non_personal_text = non_personal_text.append(replaced_text, ignore_index=True)
    non_personal_text.sort_values(by=["index"], inplace=True)

    with open(out_path_non_personal_data, "w") as out_txt:
        out_txt.write(" ".join(non_personal_text["text"].values))


def detect_text(img: typing.Any, min_conf: float) -> pd.DataFrame:
    """Recognizes text on the image.

    Args:
        img (typing.Any): image with text to be recognized.
        min_conf (float): minimum OCR Confidence Scores.

    Returns:
        pd.DataFrame: all found text on image with text field data, filtered by min_conf.
    """

    # ocr
    detected_text_df = pytesseract.image_to_data(
        img, lang="deu", output_type=Output.DATAFRAME
    )

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
    final_df = pd.DataFrame()

    max_spaces = max([value.count(" ") for value in personal_data.values()])
    for no_spaces in range(max_spaces + 1):
        tmp_df = detected_text_df.copy()
        tmp_df.rename(
            columns={"text": "text_0", "width": "width_0", "index": "index_0"},
            inplace=True,
        )

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
                # shift index column and aggreagte
                shift_colums.append("index_" + str(shift))
                highest_index_column_name = "index_" + str(shift)
                tmp_df[highest_index_column_name] = (
                    tmp_df["index_" + str(shift - 1)].astype(str)
                    + ","
                    + tmp_df.index_0.shift(-shift).fillna("").astype(str)
                )

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

            tmp_df["index_0"] = tmp_df[highest_index_column_name]
            tmp_df["text_0"] = tmp_df[highest_text_column_name]
            tmp_df["width_0"] = tmp_df[highest_width_column_name].astype(int)
            tmp_df.drop(columns=shift_colums, inplace=True)

        tmp_df.rename(
            columns={"text_0": "text", "width_0": "width", "index_0": "index"},
            inplace=True,
        )

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
        out_path_non_personal_data=snakemake.output.non_personal_data,
        min_conf=snakemake.config["min-confidence"],
        max_dist=snakemake.config["max-distance"],
        replacements_path=snakemake.params["replacements"]
    )
