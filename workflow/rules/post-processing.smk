rule summarize_found_personal_data:
    input:
        get_personal_data
    output:
        "results/{id}/personal-data-summary.tsv",
    log:
        "logs/{id}/summarize_found_personal_data.log"
    script:
        "../scripts/summarize-found-personal-data.py"