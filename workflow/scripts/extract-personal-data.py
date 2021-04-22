import json
from collections import defaultdict


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

    # TODO design this part more flexible, maybe via the snakemake config file
    # ---------------------------------------
    personal_data = defaultdict()
    personal_data["name_family"] = data.get("name")[0].get("family")
    personal_data["birthDate"] = data.get("birthDate")
    personal_data["gender"] = data.get("gender")
    personal_data["address"] = data.get("address")[0].get("line")[0]
    personal_data["city"] = " ".join(
        [data.get("address")[0].get("postalCode"), data.get("address")[0].get("city")]
    )
    personal_data["country"] = data.get("address")[0].get("country")
    personal_data["phone"] = data.get("telecom")[0].get("value")

    for i, first_name in enumerate(data.get("name")[0].get("given")):
        personal_data["name_first_{}".format(i)] = first_name

    personal_data["case_number"] = json_path.split("/")[-1].split(".")[0]
    # ---------------------------------------

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
    # TODO replace placeholders
    json_path = "SNAKEMKAE INPUT"
    out_path = "SNAKEMAKE OUTPUT"

    personal_data = parse_meta_data(json_path)

    # TODO enrich the personal data. Other examples below
    # if personal_data.get("birthDate"):
    #     personal_data = format_birthday(personal_data)

    # if personal_data.get("gender"):
    #     personal_data = format_gender(personal_data)

    # if personal_data.get("country"):
    #     personal_data = format_country(personal_data)

    save_personal_data(personal_data, out_path)


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
