rule redact_page:
    input:
        orginal_img = ""
        personal_data="results/personal-data/{id}.json",
    output:
        "results/processed-docs/{id}/{path_to_img}",
    log:
        "logs/redact_page/test.log",
    conda:
        "../envs/pytesseract.yaml"
    script:
        "../scripts/redact-page.py"
