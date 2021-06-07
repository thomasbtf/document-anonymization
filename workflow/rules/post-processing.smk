rule summarize_found_personal_data:
    input:
        data=get_personal_data,
        pages=get_processed_pages,
    output:
        temp("results/{id}/tmp/personal-data-summary.tsv"),
    log:
        "logs/{id}/summarize_found_personal_data.log",
    script:
        "../scripts/summarize-found-personal-data.py"


checkpoint create_paths_for_manually_checking:
    input:
        "results/{id}/tmp/personal-data-summary.tsv",
    output:
        no_redaction="results/{id}/tabels/no-redaction.tsv",
        high_degree_of_redaction="results/{id}/tabels/high-degree-of-redaction.tsv",
        partly_found_address="results/{id}/tabels/partly-found-address.tsv",
        partly_found_name="results/{id}/tabels/partly-found-name.tsv",
    log:
        "logs/{id}/create_paths_for_manually_checking.log",
    script:
        "../scripts/create-paths-for-manually-checking.py"


rule cp_no_redaction:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/no_redaction/{img}",
            caption="../report/no_redaction.rst",
            category="2. No redaction",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_no_redaction/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule cp_high_degree_of_redaction:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/high_degree_of_redaction/{img}",
            caption="../report/high_degree_of_redaction.rst",
            category="3. High degree of redaction",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_high_degree_of_redaction/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule cp_partly_found_address:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/partly_found_address/{img}",
            caption="../report/partly_found_address.rst",
            category="4. Partly found address",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_partly_found_address/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule cp_partly_found_name:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/partly_found_name/{img}",
            caption="../report/partly_found_name.rst",
            category="5. Partly found names",
        ),
    log:
        "logs/{id}/cp_partly_found_name/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule copy_questionable_imgs:
    input:
        lambda wildcards: get_questionable_imgs(wildcards, case="no_redaction"),
        lambda wildcards: get_questionable_imgs(
            wildcards, case="high_degree_of_redaction"
        ),
        lambda wildcards: get_questionable_imgs(wildcards, case="partly_found_address"),
        lambda wildcards: get_questionable_imgs(wildcards, case="partly_found_name"),
    output:
        temp(touch("results/{id}/tmp/moved")),
        touch("results/{id}/tmp/filler"),
    log:
        "logs/{id}/move_questionable_imgs.log",


rule remove_questionable_imgs:
    input:
        paths=lambda wildcards: get_questionable_imgs(wildcards, case="all"),
        moved="results/{id}/tmp/moved",
    output:
        temp(touch("results/{id}/tmp/deleted")),
    log:
        "logs/{id}/remove_questionable_imgs",
    shell:
        "(rm {input.paths}) 2> {log}"


rule summarize_manuel_checks:
    input:
        manuel_checks=rules.create_paths_for_manually_checking.output,
        total_imgs_processed="results/{id}/tmp/personal-data-summary.tsv",
    output:
        "results/{id}/tabels/manuel-check-summary.tsv",
    log:
        "logs/{id}/summarize_manuel_checks.log",
    script:
        "../scripts/summarize-manuel-checks.py"


rule plot_manuel_check_summary:
    input:
        "results/{id}/tabels/manuel-check-summary.tsv",
    output:
        report(
            "results/{id}/plots/summary-for-{id}.svg",
            caption="../report/manuel_check_summary.rst",
            category="1. Overview",
        ),
    log:
        "logs/{id}/plot_manuel_check_summary.log",
    conda:
        "../envs/altair.yaml"
    script:
        "../scripts/plot-manuel-check-summary.py"
