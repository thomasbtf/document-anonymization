rule extract_personal_data:
    input:
        get_fhir_metadata,
    output:
        "results/{id}/personal-data.json",
    log:
        "logs/{id}/extract_personal_data.log",
    script:
        "../scripts/extract-personal-data.py"


rule extract_lz4_docs:
    input:
        get_compressed_docs,
    output:
        directory("results/{id}/uncompressed-lz4-docs"),
    log:
        "logs/{id}/extract_lz4_docs.log",
    shell:
        "(mkdir -p {output} && lz4 -dc --no-sparse {input} | tar -xf - -C {output}) 2> {log}"


rule extract_zipped_doc:
    input:
        get_path_of_filename,
    output:
        temp(directory("results/{id}/single-uncompressed-zip-docs/{filename}")),
    log:
        "logs/{id}/extract_zip/{filename}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        '(unzip "{input}" -d "{output}") > "{log}" 2>&1'


rule extract_all_zipped_docs:
    input:
        lambda wildcards: expand(
            "results/{{id}}/single-uncompressed-zip-docs/{filename}",
            filename=get_zip_files_in_dir(wildcards),
        ),
    output:
        directory("results/{id}/uncompressed-zip-docs/"),
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
        "results/{id}/file_paths.csv",
    log:
        "logs/{id}/scan_decomp_dir.log",
    script:
        "../scripts/scan_decomp_data.py"


checkpoint fix_file_ext:
    input:
        "results/{id}/file_paths.csv",
    output:
        "results/{id}/fixed_paths.csv",
    log:
        "logs/{id}/fix_file_ext.log",
    script:
        "../scripts/fix_filenames.py"
