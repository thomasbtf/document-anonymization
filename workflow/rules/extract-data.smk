rule extract_docs:
    input:
        get_compressed_docs,
    output:
        directory("results/{id}/uncompressed-docs"),
    log:
        "logs/{id}/extract_docs.log",
    shell:
        "(mkdir -p {output} && lz4 -dc --no-sparse {input} | tar -xf - -C {output}) 2> {log}"


rule scan_decomp_dir:
    input:
        "results/{id}/uncompressed-docs",
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


rule extract_personal_data:
    input:
        get_fhir_metadata,
    output:
        "results/{id}/personal-data.json",
    log:
        "logs/{id}/extract_personal_data.log",
    script:
        "../scripts/extract-personal-data.py"
