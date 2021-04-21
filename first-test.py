import datetime

import cv2
import Levenshtein
import pandas as pd
import pytesseract
from pytesseract import Output

personal_data = {
    "last_name": "Mustermann",
    "middel_name": "Erika",
    "first_name": "Maxim",
    "birthday" : "01.01.1920",
    "street" : "Münchener Heidestraße 17",
    "city": "50667 Köln",
    "years": "99 Jahre",
    "case_number": "999999999"
}

image_file = "in_img"

out_path = "img.jpg"


def parse_page(image_path: str, out_path: str, personal_data: dict, min_conf: float):
    img = cv2.imread(image_path)

    df = detect_text(img, min_conf)

    personal_data = enrich_personal_data(personal_data)
    df = select_personal_data(df, personal_data, 2)

    img = redact(df, img)

    cv2.imwrite(out_path, img)

def enrich_personal_data(personal_data):
    personal_data["names_combined"] = ",".join([personal_data["last_name"], personal_data["middel_name"]])
    personal_data["slash_birthday"] = datetime.datetime.strptime(personal_data["birthday"], "%d.%m.%Y").strftime("%d/%m/%Y")
    personal_data = {key: value.lower() for key, value in personal_data.items()}
    personal_data = {key: value.strip() for key, value in personal_data.items()}
    return personal_data 


def detect_text(img, min_conf):
    # ocr
    detected_text_df = pytesseract.image_to_data(img, output_type=Output.DATAFRAME)

    # filter ocr table
    detected_text_df = detected_text_df[detected_text_df.conf >= min_conf]
    detected_text_df.drop(columns = ["level", "page_num", "block_num", "par_num", "line_num", "word_num"], inplace = True)
    detected_text_df.text = detected_text_df.text.str.lower()

    return detected_text_df


def select_personal_data(detected_text_df, personal_data, min_distance):

    final_df = pd.DataFrame()
    max_spaces = max([value.count(" ") for value in personal_data.values()])

    for no_spaces in range(max_spaces + 1):
        tmp_df = detected_text_df.copy()
        tmp_df.rename(columns={"text": "text_0", "width" : "width_0"}, inplace = True)

        # subset of personal data dict, according to # spaces
        tmp_dict = {key: value.replace(" ", "") for key, value in personal_data.items() if value.count(" ") == no_spaces}

        # shif text to get longer phrases
        if no_spaces > 0:
            for shift in range(1, no_spaces + 1):
                # shift text column and aggreagte
                highest_text_column_name = "text_" + str(shift)
                tmp_df[highest_text_column_name] = tmp_df["text_" + str(shift - 1)] + tmp_df.text_0.shift(-shift).fillna("")

                # shift width column and aggreagte
                highest_width_column_name = "width_" + str(shift)
                tmp_df[highest_width_column_name] = tmp_df["width_" + str(shift - 1)] + tmp_df.width_0.shift(-shift).fillna(0) + tmp_df.width_0 / tmp_df.text_0.str.len()
            
            tmp_df["width_0"] = tmp_df[highest_width_column_name].astype(int)
        else:
            highest_text_column_name = "text_0"
        
        # calc edit distances for each key in subsampled dict
        for key, value in tmp_dict.items():
            tmp_df[key] = tmp_df[highest_text_column_name].apply(lambda text: Levenshtein.distance(value, text))

        tmp_df.rename(columns={"text_0": "text", "width_0":"width"}, inplace = True)

        # select entries, where distance is below or equal a threshhold
        query_list = []
        for col in tmp_dict.keys():
            query_list.append(f'{col}<={min_distance}') 
        tmp_df = tmp_df.query(' | '.join(query_list))
        
        final_df = final_df.append(tmp_df, ignore_index=True)
    
    return final_df[["left", "top", "width", "height", "conf", "text"]]


def redact(personal_data_df, img):
    for i in personal_data_df.index:
        (x, y, w, h) = personal_data_df.loc[i].left, personal_data_df.loc[i].top, personal_data_df.loc[i].width, personal_data_df.loc[i].height, 
        img = cv2.rectangle(img, (x, y), (x + w, y + h), (0, 0, 0), -1)
    return img


if __name__ == "__main__":
    parse_page(image_file, out_path, personal_data, 0.6)
