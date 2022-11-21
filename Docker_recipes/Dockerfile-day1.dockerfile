FROM condaforge/mambaforge:latest

RUN mkdir /conda-envs && \
    mamba create -p /conda-envs/bioinfo_tools -c conda-forge -c bioconda samtools bcftools delly alfred && \
    mamba create -p /conda-envs/r_vizu -c conda-forge -c bioconda -c r r-rcircos r-irkernel jupyterlab && \
    mamba clean --all -y

