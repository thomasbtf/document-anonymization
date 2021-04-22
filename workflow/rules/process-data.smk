rule redact_page:
    input:
        original_img="results/uncompressed-docs/1/1021076867.20210420.docs/2011-08-05.25500363.WithoutTitle/0",
        personal_data="results/personal-data/1.json",
    output:
        "results/processed-docs/1/1021076867.20210420.docs/2011-08-05.25500363.WithoutTitle/0",
    log:
        "logs/redact_page/test.log",
    conda:
        "../envs/pytesseract.yaml"
    script:
        "../scripts/redact-page.py"
