rule preprocess_page:
    input:
        "results/{id}/uncompressed-docs/{img}",
    output:
        "results/{id}/preprocessed-docs/{img}",
    log:
        "logs/{id}/preprocess-page/{img}.log",
    conda:
        "../envs/opencv.yaml"
    script:
        "../scripts/preprocess-page.py"


rule identify_personal_data:
    input:
        preprocessed_page="results/{id}/preprocessed-docs/{img}",
        personal_data="results/{id}/personal-data.json",
    output:
        text_to_redact="results/{id}/data-to-redact/{img}.tsv",
        all_text="results/{id}/detected-text/{img}.tsv",
    log:
        "logs/{id}/identify-personal-data/{img}.log",
    conda:
        "../envs/pytesseract.yaml"
    script:
        "../scripts/identify-personal-data.py"


rule redact_page:
    input:
        orginal_page="results/{id}/uncompressed-docs/{img}",
        data_to_redact="results/{id}/data-to-redact/{img}.tsv",
    output:
        "results/{id}/processed-docs/{img}",
    log:
        "logs/{id}/redact-page/{img}.log",
    conda:
        "../envs/opencv.yaml"
    script:
        "../scripts/redact-page.py"
