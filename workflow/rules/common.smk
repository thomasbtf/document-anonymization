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


def get_processed_pages(wildcards):
    with checkpoints.scan_decomp_dir.get(id=wildcards.id).output[0].open() as f:
        paths = pd.read_csv(f, sep="\n", header=None, squeeze=True)

    paths = [
        path.replace("results/{id}/uncompressed-docs/".format(id=wildcards.id), "")
        for path in paths
        if ".snakemake" not in path
    ]

    pattern = "results/{{id}}/processed-docs/{img}"
    return expand(pattern, img=paths)
