rule preprocess_page:
    input:
        get_uncompressed_image,
    output:
        temp("results/{id}/tmp/preprocessed-docs/{img}"),
    log:
        "logs/{id}/preprocess-page/{img}.log",
    conda:
        "../envs/opencv.yaml"
    script:
        "../scripts/preprocess-page.py"


rule identify_personal_data:
    input:
        preprocessed_page="results/{id}/tmp/preprocessed-docs/{img}",
        personal_data="results/{id}/tmp/personal-data.json",
    output:
        text_to_redact=temp("results/{id}/tmp/data-to-redact/{img}.tsv"),
        all_text=temp("results/{id}/tmp/detected-text/{img}.tsv"),
        non_personal_data="results/{id}/text/{img}.txt",
    params:
        replacements = "resources/replacements.json"
    log:
        "logs/{id}/identify-personal-data/{img}.log",
    conda:
        "../envs/pytesseract.yaml"
    script:
        "../scripts/identify-personal-data.py"


rule redact_page:
    input:
        orginal_page=get_uncompressed_image,
        data_to_redact="results/{id}/tmp/data-to-redact/{img}.tsv",
    output:
        "results/{id}/processed-docs/{img}",
    params:
        version=get_version(),
    log:
        "logs/{id}/redact-page/{img}.log",
    conda:
        "../envs/opencv.yaml"
    script:
        "../scripts/redact-page.py"
