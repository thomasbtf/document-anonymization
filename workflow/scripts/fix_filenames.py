import os
import sys
import filetype
import pdf2image


def add_ext(paths_file: str, fixed_paths: str):
    with open(paths_file, "r") as path_list:
        with open(fixed_paths, "w") as new_paths:
            cwd = os.getcwd()
            accepted_ext = ["jpg", "jpeg", "tiff", "tif", "bmp"]
            ext_pairs = [{"jpg", "jpeg"}, {"tiff", "tif"}]

            for path in path_list:
                path = path.strip()
                filedir, filename = os.path.split(path)
                ext = filename.split(".")[-1]
                ftype = filetype.guess(path).extension

                if ext != ftype:
                    #  file won´t be written to file-list for further processing
                    if ftype is None:
                        print(f"file {path} is in an incompatible file format.")

                    # convert pdf
                    elif ftype == "pdf":
                        print("pdf")
                        pages = pdf2image.convert_from_path(path)
                        for i, page in enumerate(pages):
                            new_paths.write(path + f"_{i}.tif\n")
                        os.chdir(os.path.join(cwd, filedir))
                        for i, page in enumerate(pages):
                            page.save(f"{filename}_{i}.tif", "TIFF")
                        os.chdir(cwd)
                        os.remove(path)

                    # change the img file type
                    elif ext not in accepted_ext and ftype in accepted_ext:
                        print("change extension")
                        new_paths.write(path + f".{ftype}\n")
                        os.chdir(os.path.join(cwd, filedir))
                        os.rename(filename, filename + f".{ftype}")
                        os.chdir(cwd)

                    # This elif clause allows to leave files with alternative but adequate extension untouched
                    elif set((ext, ftype)) in ext_pairs:
                        print("set((ext, ftype)) in ext_pairs")
                        new_paths.write(path + "\n")

                    # This elif clause allows to specify what shall happen to the file in question.
                    # This will become relevant if we are going to allow non-img file-types like pdf.
                    # In this case the file needs to be channeled into another branch of the workflow.
                    elif ftype == "None" and ext in accepted_ext:
                        print(f"file {path} is not an image file.")

                # file extension equals equals the detected extension
                else:
                    print("else")
                    new_paths.write(path + "\n")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    add_ext(snakemake.input.files, snakemake.output[0])
