FROM condaforge/mambaforge:latest

RUN mkdir /conda-envs && \
    mamba install -c conda-forge jupyterlab nb_conda_kernels && \
    mamba create -p /conda-envs/bioinfo_tools -c conda-forge -c bioconda samtools bcftools delly alfred ipykernel && \
    mamba create -p /conda-envs/r_vizu -c conda-forge -c bioconda -c r r-rcircos ipykernel && \
    mamba clean --all -y

