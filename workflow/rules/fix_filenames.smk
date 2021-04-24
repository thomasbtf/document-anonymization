checkpoint fix_file_ext:
    input:
        "results/{id}/file_paths.csv",
    output:
        "results/{id}/fixed_paths.csv",
    log:
        "logs/{id}/fix_file_ext.log",
    script:
        "../scripts/fix_filenames.py"
