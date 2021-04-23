rule scan_decomp_dir:
    input:
        directory("results/uncompressed-docs/{id}")
    output:
        "results/uncompressed-docs/{id}_file_paths.csv"
    log:
        "logs/scan_decomp_dir/{id}.log"
    script:
        "../scripts/scan_decomp_data.py"