rule summarize_found_personal_data:
    input:
        data=get_personal_data,
        pages=get_processed_pages,
    output:
        "results/{id}/personal-data-summary.tsv",
    log:
        "logs/{id}/summarize_found_personal_data.log",
    script:
        "../scripts/summarize-found-personal-data.py"
