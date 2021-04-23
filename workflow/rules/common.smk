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


def get_compressed_docs(wildcards):
    return pep.sample_table.loc[wildcards.id][["compressed_docs"]]


def get_fhir_metadata(wildcards):
    return pep.sample_table.loc[wildcards.id][["fhir_metadata"]]
