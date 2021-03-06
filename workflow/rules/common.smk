from snakemake.utils import validate
import pandas as pd
from os import walk


# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"


##### load config and sample sheets #####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv("config/pep/documents.csv").set_index("sample_name", drop=False)
samples.index.names = ["sample_name"]
validate(samples, schema="../schemas/documents.schema.yaml")


##### helper functions #####


def get_compressed_docs(wildcards):
    return pep.sample_table.loc[wildcards.id][["compressed_docs"]]


def get_fhir_metadata(wildcards):
    return pep.sample_table.loc[wildcards.id][["fhir_metadata"]]

def get_additional_metadata(wildcards):
    if "additional_metadata" in pep.sample_table.columns:
        if type(pep.sample_table.loc[wildcards.id]["additional_metadata"]) == str:
            return pep.sample_table.loc[wildcards.id][["additional_metadata"]]
        else:
            return []
    else:
        return []

def get_all_ids():
    return pep.sample_table["sample_name"].to_list()

def get_image_paths_for_id(wildcards):
    with checkpoints.fix_file_ext.get(id=wildcards.id).output[0].open() as f:
        paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    # TODO
    # This is ugly. We should rewrite the file handling.
    paths = [
        path.removeprefix(
            "results/{id}/tmp/uncompressed-lz4-docs/".format(id=wildcards.id)
        )
        for path in paths
        if ".snakemake" not in path and ".DS_Store" not in path
    ]

    paths = [
        path.removeprefix(
            "results/{id}/tmp/uncompressed-zip-docs/".format(id=wildcards.id)
        )
        for path in paths
    ]
    return paths


def get_processed_pages(wildcards):
    paths = get_image_paths_for_id(wildcards)
    pattern = "results/{{id}}/processed-docs/{img}"
    return expand(pattern, img=paths)


def get_personal_data(wildcards):
    paths = get_image_paths_for_id(wildcards)
    pattern = "results/{{id}}/tmp/data-to-redact/{img}.tsv"
    return expand(pattern, img=paths)


def load_tsv(f):
    try:
        return pd.read_csv(f, sep="\n", header=None, squeeze=True)
    except pd.errors.EmptyDataError:
        return []


def get_questionable_imgs(wildcards, case):
    if case == "no_redaction":
        pattern = "results/{id}/to-check/no_redaction/{{img}}".format(id=wildcards.id)
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.no_redaction.open() as f:
            paths = load_tsv(f)

    elif case == "high_degree_of_redaction":
        pattern = "results/{id}/to-check/high_degree_of_redaction/{{img}}".format(
            id=wildcards.id
        )
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.high_degree_of_redaction.open() as f:
            paths = load_tsv(f)

    elif case == "partly_found_address":
        pattern = "results/{id}/to-check/partly_found_address/{{img}}".format(
            id=wildcards.id
        )
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.partly_found_address.open() as f:
            paths = load_tsv(f)

    elif case == "partly_found_name":
        pattern = "results/{id}/to-check/partly_found_name/{{img}}".format(
            id=wildcards.id
        )
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.partly_found_name.open() as f:
            paths = load_tsv(f)
    elif case == "all":
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.no_redaction.open() as f:
            no_redaction = load_tsv(f)
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.high_degree_of_redaction.open() as f:
            high_degree_of_redaction = load_tsv(f)
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.partly_found_address.open() as f:
            partly_found_address = load_tsv(f)
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.partly_found_name.open() as f:
            partly_found_name = load_tsv(f)

        all_paths = pd.Series("results/{id}/tmp/filler".format(id=wildcards.id))
        all_paths = (
            all_paths.append(no_redaction, ignore_index=True)
            .append(high_degree_of_redaction, ignore_index=True)
            .append(partly_found_address, ignore_index=True)
            .append(partly_found_name, ignore_index=True)
        )
        all_paths = all_paths.unique()

        return all_paths

    paths = [
        path.replace("results/{id}/processed-docs/".format(id=wildcards.id), "")
        for path in paths
    ]

    return expand(pattern, img=paths)


def get_resource(name):
    return str((Path(workflow.snakefile).parent.parent.parent / "resources") / name)


def get_version():
    with open("resources/version.txt", "r") as v:
        return v.readline().replace("\n", "")


def get_uncompressed_docs_dir(wildcards):
    docs = get_compressed_docs(wildcards).values[0]
    if docs.endswith(".docs"):
        return ("results/{id}/tmp/uncompressed-zip-docs",)
    elif docs.endswith(".tar.lz4"):
        return ("results/{id}/tmp/uncompressed-lz4-docs",)


def get_zip_files_in_dir(wildcards):
    docs = get_compressed_docs(wildcards).values[0]
    filenames = os.listdir(docs)
    return [
        filename.removesuffix(".zip")
        for filename in filenames
        if filename.endswith(".zip")
    ]


def get_path_of_filename(wildcards):
    docs = get_compressed_docs(wildcards).values[0]
    return os.path.join(docs, wildcards.filename + ".zip")


def get_uncompressed_image(wildcards):
    dir = get_uncompressed_docs_dir(wildcards)[0]
    return os.path.join(dir, wildcards.img)
