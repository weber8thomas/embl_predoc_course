bam_input = ["BM510_chr10_chr13_chr22.bam", "RPEWtp30_chr10_chr13_chr22.bam"]
ref_genome = "hg38.fa"
sv_types = ["TRA", "INV", "DUP", "INS", "DEL"]
samples_file = "spl.tsv"


rule all:
    input:
        expand("DELLY_FILTER/somatic_{sv_type}.bcf", sv_type=sv_types, allow_missing=True)

rule delly_call:
    input:
        bam=bam_input,
        ref=ref_genome,
    output:
        "DELLY_CALL/{sv_type}.bcf"
    log:
        "log/delly_call/{sv_type}.log"
    envmodules:
        "delly/0.7.6-foss-2016b"
    shell:
        """
        delly call \
            --type {wildcards.sv_type} \
            --noindels \
            --map-qual 20 \
            --genome {input.ref} \
            --outfile {output} \
            {input.bam} > {log}
        """

rule delly_filter:
    input:
        bcf="DELLY_CALL/{sv_type}.bcf",
        samples=samples_file
    output:
        "DELLY_FILTER/somatic_{sv_type}.bcf"
    log:
        "log/delly_filter/{sv_type}.log"
    envmodules:
        "delly/0.7.6-foss-2016b"
    shell:
        """
        delly filter \
            --type {wildcards.sv_type} \
            --pass \
            --filter somatic \
            --altaf 0.25 \
            --samples {input.samples} \
            --outfile {output} \
            {input.bcf} > {log}
        """

