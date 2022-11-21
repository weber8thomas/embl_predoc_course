# install dependencies
mamba install -c conda-forge jupyterlab nb_conda_kernels

mamba create -n bioinfo_tools -c conda-forge -c bioconda samtools bcftools delly alfred &&
    mamba create -n r_vizu -c conda-forge -c bioconda -c r r-rcircos r-irkernel &&
    mamba clean --all -y
