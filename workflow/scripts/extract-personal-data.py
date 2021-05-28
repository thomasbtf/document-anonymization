sys.stderr = open(snakemake.log[0], "w")

import json
from collections import defaultdict
import sys
import itertools


def parse_meta_data(json_path: str) -> defaultdict:
    """Parses the FHIR metadata and extracts personal data.
    The extracted data is redacted in the further course of the workflow.

    Args:
        json_path (str): path to FHIR metadata

    Returns:
        defaultdict: personal data of the patient
    """

    with open(json_path) as json_file:
        data = json.load(json_file)

    # select the patient reosurce from the bundel data export
    for ele in data.get("entry", {}):
        # iterate of entries
        for key, value in ele.get("resource", {}).items():
            if key == "resourceType" and value == "Patient":
                data = ele.get("resource")
                break

    # TODO design this part more flexible, maybe via the snakemake config file
    # ---------------------------------------
    personal_data = defaultdict()
    first_name_count = 0
    for i, first_name in enumerate(data.get("name")[0].get("given")):
        first_name_count += 1
        personal_data["name_first_{}".format(i)] = first_name
    personal_data["name_family"] = data.get("name")[0].get("family")
    personal_data["birthDate"] = data.get("birthDate")
    personal_data["address"] = data.get("address")[0].get("line")[0]
    personal_data["city"] = " ".join([data.get("address")[0].get("postalCode"), data.get("address")[0].get("city")])
    personal_data["case_number"] = json_path.split("/")[-1].split(".")[0]
    for com in data.get("telecom", {}):
        com_type = com.get("system", {})
        personal_data[com_type] = com.get("value", {})
    # personal_data["gender"] = data.get("gender")
    # personal_data["country"] = data.get("address")[0].get("country")
    # ---------------------------------------
    return personal_data, first_name_count

def variate_personal_data(personal_data: dict, first_name_count: int) -> defaultdict:
    # permutate names
    names_simple = set((personal_data["name_first_0"], personal_data["name_family"]))
    names_all = set()
    for i in range(first_name_count):
        names_all.add(personal_data["name_first_{}".format(i)])
    names_all.add(personal_data["name_family"])

    name_perms = list(itertools.permutations(list(names_simple)))
    if names_simple != names_all:
        names_all_perm = list(itertools.permutations(list(names_all)))
        name_perms.extend(names_all_perm)
    
    for i, perm in enumerate(name_perms):
        personal_data[f"name_perm_{i}"] = ",".join(perm)
    
    # variate phone number
    provider_local_codes = [
                            "01511", "01512", "01514", "01515", "01516", "01517", "01520", "01522", "01523", 
                            "01525", "015566", "01570", "01573", "01575", "01577", "01578", "01590", "0160",
                             "0162", "0163", "0170", "0171", "0172", "0173", "0174", "0175", "0176", "0177",
                            "0178", "0179", 
                            ]

    # it would be much better to generate this list only once centrally instead for every patient sample again
    with open("resources/Vorwahlen_Festnetz_Bundesnetzagentur.csv", "r") as local_codes:
        for line in local_codes:
            if line.startswith("Ortsnetzkennzahl"):
                pass
            else:
                provider_local_codes.append("0" + line.split(";")[0])
    
    nums = ["1","2","3","4","5","6","7","8","9","0"]
    tmp_phone = personal_data.get("phone", "")
    for letter in tmp_phone:
        if letter not in nums:
            tmp_phone = personal_data["phone"].replace(letter, "")
    personal_data["phone_perm0"] = tmp_phone

    for code in provider_local_codes:
        if  tmp_phone.startswith(code):
            pre_code = code
            break
        else:
            pre_code = tmp_phone[:4]
    
    seperators = ["/", "\\", "-", " ", "_", ".", ":"]

    for i, sep in enumerate(seperators):
       personal_data[f"phone_perm{i+1}"] = tmp_phone[:len(pre_code)+1] + sep + tmp_phone[len(pre_code)+1:]

    #variate birthdate
    yr, m, dy = personal_data["birthDate"].split("-")
    for i, sep in enumerate(seperators):
        personal_data[f"birthDate_perm{i}"] = f"{dy}{sep}{m}{sep}{yr}"
        personal_data[f"birthDate_perm{i}"] = f"*{dy}{sep}{m}{sep}{yr}"
        personal_data[f"birthDate_perm{i}{i}"] = f"{yr}{sep}{m}{sep}{dy}"
    
    #variate country

    return personal_data


def save_personal_data(personal_data: dict, out_path: str):
    """Save the final dic with the personal data as json.

    Args:
        personal_data (dict): dict with the personal data, that is to be removed
        out_path (str): path to save the json to
    """

    with open(out_path, "w") as fp:
        json.dump(personal_data, fp, indent=2)


if __name__ == "__main__":
    personal_data = parse_meta_data(snakemake.input[0])
    var_data = variate_personal_data(personal_data[0], personal_data[1])
    # TODO enrich the personal data. Other examples below
    # if personal_data.get("birthDate"):
    #     personal_data = format_birthday(personal_data)

    # if personal_data.get("gender"):
    #     personal_data = format_gender(personal_data)

    # if personal_data.get("country"):
    #     personal_data = format_country(personal_data)

    #personal_data = {key: value.lower().strip() for key, value in personal_data.items()}
    var_data = {key: value.lower().strip() for key, value in var_data.items()}
    save_personal_data(var_data, snakemake.output[0])

# def enrich_personal_data(personal_data: dict) -> dict:
#     personal_data["names_combined"] = ",".join(
#         [personal_data["last_name"], personal_data["middel_name"]]
#     )
#     personal_data["slash_birthday"] = datetime.datetime.strptime(
#         personal_data["birthday"], "%d.%m.%Y"
#     ).strftime("%d/%m/%Y")
#     personal_data = {key: value.lower() for key, value in personal_data.items()}
#     personal_data = {key: value.strip() for key, value in personal_data.items()}
#     return personal_data
