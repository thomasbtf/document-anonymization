rule extract_docs:
    input:
        get_compressed_docs,
    output:
        directory("results/uncompressed-docs/{id}"),
    log:
        "logs/extract_docs/{id}.log",
    shell:
        "(mkdir -p {output} && lz4 -dc --no-sparse {input} | tar -xf - -C {output}) 2> {log}"


rule extract_personal_data:
    input:
        get_fhir_metadata,
    output:
        "results/personal-data/{id}.json",
    log:
        "logs/extract_personal_data/{id}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract-personal-data.py"
