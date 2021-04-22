rule extract_docs:
    input:
        get_compressed_docs
    output:
        directory("results/uncompressed-docs/{id}"),
    log:
        "logs/extract_docs/{id}.log"
    shell:
        "(mkdir -p {output} && lz4 -dc --no-sparse {input} | tar -xf - -C {output}) 2> {log}"