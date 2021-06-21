rule extract_personal_data:
    input:
        get_fhir_metadata,
        get_additional_metadata,
    output:
        temp("results/{id}/tmp/personal-data.json"),
    log:
        "logs/{id}/extract_personal_data.log",
    script:
        "../scripts/extract-personal-data.py"

# if get_additional_metadata

rule extract_lz4_docs:
    input:
        get_compressed_docs,
    output:
        directory("results/{id}/tmp/uncompressed-lz4-docs"),
    log:
        "logs/{id}/extract_lz4_docs.log",
    shell:
        "(mkdir -p {output} && lz4 -dc --no-sparse {input} | tar -xf - -C {output}) 2> {log}"


rule extract_zipped_doc:
    input:
        get_path_of_filename,
    output:
        temp(directory("results/{id}/tmp/single-uncompressed-zip-docs/{filename}")),
    log:
        "logs/{id}/extract_zip/{filename}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        '(unzip "{input}" -d "{output}") > "{log}" 2>&1'


rule extract_all_zipped_docs:
    input:
        lambda wildcards: expand(
            "results/{{id}}/tmp/single-uncompressed-zip-docs/{filename}",
            filename=get_zip_files_in_dir(wildcards),
        ),
    output:
        directory("results/{id}/tmp/uncompressed-zip-docs/"),
    params:
        in_dir=lambda w, input: os.path.dirname(input[0]),
    log:
        "logs/{id}/extract_zip/all.log",
    shell:
        "(mkdir -p {output} && cp -r  {params.in_dir}/* {output}) 2> {log}"


rule scan_decomp_dir:
    input:
        get_uncompressed_docs_dir,
    output:
        temp("results/{id}/tmp/file_paths.csv"),
    log:
        "logs/{id}/scan_decomp_dir.log",
    script:
        "../scripts/scan_decomp_data.py"


# TODO combine this rule with the rule scan_decomp_dir
checkpoint fix_file_ext:
    input:
        files="results/{id}/tmp/file_paths.csv",
        dirs=get_uncompressed_docs_dir,
    output:
        temp("results/{id}/tmp/fixed_paths.csv"),
    log:
        "logs/{id}/fix_file_ext.log",
    conda:
        "../envs/filetype.yaml"
    script:
        "../scripts/fix_filenames.py"
