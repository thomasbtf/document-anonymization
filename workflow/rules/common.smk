from snakemake.utils import validate
import pandas as pd


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


def get_all_ids():
    return pep.sample_table["sample_name"].to_list()


def get_image_paths_for_id(wildcards):
    with checkpoints.fix_file_ext.get(id=wildcards.id).output[0].open() as f:
        paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    paths = [
        path.replace("results/{id}/uncompressed-docs/".format(id=wildcards.id), "")
        for path in paths
        if ".snakemake" not in path
    ]

    return paths


def get_processed_pages(wildcards):
    paths = get_image_paths_for_id(wildcards)
    pattern = "results/{{id}}/processed-docs/{img}"
    return expand(pattern, img=paths)


def get_personal_data(wildcards):
    paths = get_image_paths_for_id(wildcards)
    pattern = "results/{{id}}/data-to-redact/{img}.tsv"
    return expand(pattern, img=paths)


def get_questionable_imgs(wildcards, case):
    if case == "no_redaction":
        pattern = "results/{id}/to-check/no_redaction/{{img}}".format(id=wildcards.id)
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.no_redaction.open() as f:
            paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    elif case == "high_degree_of_redaction":
        pattern = "results/{id}/to-check/high_degree_of_redaction/{{img}}".format(
            id=wildcards.id
        )
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.high_degree_of_redaction.open() as f:
            paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    elif case == "partly_found_address":
        pattern = "results/{id}/to-check/partly_found_address/{{img}}".format(
            id=wildcards.id
        )
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.partly_found_address.open() as f:
            paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    elif case == "partly_found_name":
        pattern = "results/{id}/to-check/partly_found_name/{{img}}".format(
            id=wildcards.id
        )
        with checkpoints.create_paths_for_manually_checking.get(
            id=wildcards.id
        ).output.partly_found_name.open() as f:
            paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    paths = [
        path.replace("results/{id}/processed-docs/".format(id=wildcards.id), "")
        for path in paths
    ]

    return expand(pattern, img=paths)
