import os
import sys
import imghdr


def add_ext(paths_file: str, fixed_paths: str):
    with open(paths_file, 'r') as path_list:
        with open(fixed_paths, 'w') as new_paths:
            cwd = os.getcwd()
            accepted_ext = ["jpg", "jpeg", "tiff", "tif", "bmp"] # may add "pdf" etc.
            ext_pairs = [{"jpg", "jpeg"}, {"tiff", "tif"}]
            for path in path_list:
                path = path.strip()
                filedir, filename = os.path.split(path)
                ext = filename.split(".")[-1]
                ftype = imghdr.what(path)
                if ext != ftype:
                    if ext not in accepted_ext:
                        new_paths.write(path + f".{ftype}\n")
                        os.chdir(os.path.join(cwd, filedir))
                        os.rename(filename, filename + f".{ftype}")
                        os.chdir(cwd)
                    elif set((ext, ftype)) in ext_pairs:
                        print(ext, ftype)
                        new_paths.write(path + "\n")
                    elif ftype == "None" and ext in accepted_ext:
                        print(f"file {path} is not an image file.")
                        # This elif clause allows to specify what shall happen to the file in question.
                        # This will become relevant if we are going to allow non-img file-types like pdf.
                        # In this case the file needs to be channeled into another branch of the workflow.
                    else:
                        print(f"file {path} is in an incompatible file format.")
                        # file won´t be written to file-list for further processing
                else:
                    new_paths.write(path + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    add_ext(snakemake.input[0], snakemake.output[0])