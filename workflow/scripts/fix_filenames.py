import os
import sys


def add_ext(paths_file: str, fixed_paths: str):
    with open(paths_file, 'r') as path_list:
        with open(fixed_paths, 'w') as new_paths:
            cwd = os.getcwd()
            accepted_ext = ["jpg", "jpeg", "tiff", "tif", "bmp"]
            for path in path_list:
                path = path.strip()
                filedir, filename = os.path.split(path)
                ext = filename.split(".")[-1]
                if ext not in accepted_ext:
                    new_paths.write(path.strip() + ".jpg\n")
                    os.chdir(os.path.join(cwd, filedir))
                    os.rename(filename, filename + ".jpg")
                    os.chdir(cwd)
                else:
                    new_paths.write(path + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    add_ext(snakemake.input[0], snakemake.output[0])
    