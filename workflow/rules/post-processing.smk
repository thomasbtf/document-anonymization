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


checkpoint create_paths_for_manually_checking:
    input:
        "results/{id}/personal-data-summary.tsv",
    output:
        no_personal_data_found="results/{id}/no_personal_data_found.tsv",
        lots_personal_data_found="results/{id}/lots_personal_data_found.tsv",
        address_not_entirely_found="results/{id}/address_not_entirely_found.tsv",
        name_not_entirely_found="results/{id}/name_not_entirely_found.tsv",
    log:
        "logs/{id}/create_paths_for_manually_checking.log",
    script:
        "../scripts/create-paths-for-manually-checking.py"


rule cp_no_personal_data_found:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/no_personal_data_found/{img}",
            caption="../report/test.rst",
            category="No Data Found",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_no_personal_data_found/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule cp_lots_personal_data_found:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/lots_personal_data_found/{img}",
            caption="../report/test.rst",
            category="Lot of Data Found",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_lots_personal_data_found/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule cp_address_not_entirely_found:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/address_not_entirely_found/{img}",
            caption="../report/test.rst",
            category="Address Not Entirely Found",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_address_not_entirely_found/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule cp_name_not_entirely_found:
    input:
        "results/{id}/processed-docs/{img}",
    output:
        report(
            "results/{id}/to-check/name_not_entirely_found/{img}",
            caption="../report/test.rst",
            category="Name Not Entirely Found",
            subcategory="{id}",
        ),
    log:
        "logs/{id}/cp_name_not_entirely_found/{img}.log",
    shell:
        "(cp '{input}' '{output}') 2> '{log}'"


rule move_questionable_imgs:
    input:
        lambda wildcards: get_questionable_imgs(
            wildcards, case="no_personal_data_found"
        ),
        lambda wildcards: get_questionable_imgs(
            wildcards, case="lots_personal_data_found"
        ),
        lambda wildcards: get_questionable_imgs(
            wildcards, case="address_not_entirely_found"
        ),
        lambda wildcards: get_questionable_imgs(
            wildcards, case="name_not_entirely_found"
        ),
    output:
        touch("results/{id}/moved"),
    log:
        "logs/{id}/move_questionable_imgs.log",
