rule preprocess_page:
    input:
        "results/uncompressed-docs/{id}/{path_to_img}",
    output:
        deskewed_image="results/preprocessed-docs/{id}/{path_to_img}",
        processed_image="results/preprocessed-docs/{id}/{path_to_img}/placeholder",
    log:
        "logs/preprocess-page/{id}/{path_to_img}.log",
    conda:
        "../envs/opencv.yaml"
    script:
        "../scripts/preprocess-page.py"


rule redact_page:
    input:
        orginal_img="results/uncompressed-docs/{id}/{path_to_img}",
        personal_data="results/personal-data/{id}.json",
    output:
        "results/processed-docs/{id}/{path_to_img}",
    log:
        "logs/redact_page/{id}/{path_to_img}.log",
    conda:
        "../envs/pytesseract.yaml"
    script:
        "../scripts/redact-page.py"
