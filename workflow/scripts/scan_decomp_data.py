import os
import sys
import typing


def scan_folder(subfolder_path: str, writeable_file_object: typing.TextIO):
    ignore = [".snakemake_timestamp", ".DS_Store"]
    
    for entry in os.scandir(subfolder_path):
        if entry.is_dir(follow_symlinks=False):
            scan_folder(entry.path, writeable_file_object)
        elif entry.is_file() and not any(
            ignore_element in entry.path for ignore_element in ignore
        ):
            writeable_file_object.write(f"{entry.path}\n")
        else:
            pass


def recursive_folder_scan(decomp_data_dir: str, results_csv_paths: str):
    with open(results_csv_paths, "w") as paths_csv:
        scan_folder(decomp_data_dir, paths_csv)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    recursive_folder_scan(snakemake.input[0], snakemake.output[0])
