FROM condaforge/mambaforge:latest

RUN mkdir /conda-envs

RUN mamba create -n bioinfo_tools -p /conda-envs/bioinfo_tools -c conda-forge -c bioconda samtools bcftools delly alfred && \
    mamba create -n r_vizu -p /conda-envs/r_vizu -c conda-forge -c bioconda -c r r-circos && \
    mamba clean --all -y
