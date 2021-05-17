FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="f17b7e4073b7eede2245a158bc1adbcc168049b0d4c3bef010ebf235c00be057"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/altair.yaml
#   prefix: /conda-envs/967ce1dbbaa66e76e6e805900689443e
#   channels:
#     - conda-forge
#     - anaconda
#   dependencies:  
#     - altair =4.1
#     - altair_saver =0.5
RUN mkdir -p /conda-envs/967ce1dbbaa66e76e6e805900689443e
COPY workflow/envs/altair.yaml /conda-envs/967ce1dbbaa66e76e6e805900689443e/environment.yaml

# Conda environment:
#   source: workflow/envs/opencv.yaml
#   prefix: /conda-envs/36cefac058e5ea5e6c75c7aa5ceb6041
#   channels:
#     - conda-forge
#   dependencies:
#     - opencv =4.5.1
#     - numpy =1.20.2
RUN mkdir -p /conda-envs/36cefac058e5ea5e6c75c7aa5ceb6041
COPY workflow/envs/opencv.yaml /conda-envs/36cefac058e5ea5e6c75c7aa5ceb6041/environment.yaml

# Conda environment:
#   source: workflow/envs/pytesseract.yaml
#   prefix: /conda-envs/000b5fd8f54e617a93595385111b2095
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#   dependencies:
#     - python =3.9.2
#     - pytesseract =0.3.7
#     - tesseract =4.1.1
#     - opencv =4.5.1
#     - python-levenshtein =0.12.2
#     - pandas =1.2.4
RUN mkdir -p /conda-envs/000b5fd8f54e617a93595385111b2095
COPY workflow/envs/pytesseract.yaml /conda-envs/000b5fd8f54e617a93595385111b2095/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/967ce1dbbaa66e76e6e805900689443e --file /conda-envs/967ce1dbbaa66e76e6e805900689443e/environment.yaml && \
    mamba env create --prefix /conda-envs/36cefac058e5ea5e6c75c7aa5ceb6041 --file /conda-envs/36cefac058e5ea5e6c75c7aa5ceb6041/environment.yaml && \
    mamba env create --prefix /conda-envs/000b5fd8f54e617a93595385111b2095 --file /conda-envs/000b5fd8f54e617a93595385111b2095/environment.yaml && \
    mamba clean --all -y
