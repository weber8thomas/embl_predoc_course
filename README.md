# SV computational practical - Predoc course

- [SV computational practical - Predoc course](#sv-computational-practical---predoc-course)
  - [Introduction](#introduction)
  - [0. Prerequesite](#0-prerequesite)
  - [1. Quality checking of the BAM file](#1-quality-checking-of-the-bam-file)
  - [2. Structural Variant Calling](#2-structural-variant-calling)
    - [Example of duplication calling](#example-of-duplication-calling)
    - [Check the result of duplication calling](#check-the-result-of-duplication-calling)
  - [3. Somatic Filtering](#3-somatic-filtering)
    - [Preparation of sample list](#preparation-of-sample-list)
    - [Example of somatic duplication filtering](#example-of-somatic-duplication-filtering)
    - [Check the result of somatic duplication filtering](#check-the-result-of-somatic-duplication-filtering)
  - [4. Visualize SVs in circos plot](#4-visualize-svs-in-circos-plot)
  - [5. Detection of subclonal SVs](#5-detection-of-subclonal-svs)
    - [Call TRA using mixture data](#call-tra-using-mixture-data)
    - [Check the result of TRA calling](#check-the-result-of-tra-calling)
    - [Filtering somatic TRA using mixture data](#filtering-somatic-tra-using-mixture-data)
    - [Check the result of somatic TRA calling](#check-the-result-of-somatic-tra-calling)

## Introduction

In this practical we will learn how to identify somatic structural variations using whole genome sequencing of bulk samples (WGS) and Delly (Rausch et al. 2012). We will start from the bam file which includes three example chromosomes (chr10, chr13, chr22) and then make somatic SV calling.

This practical is largely based on the tutorial from Tobias Rausch for SV calling

https://tobiasrausch.com/courses/vc/sv/

We will firstly start warming up with some basic linux command lines we will frequently use during the practical session.

## 0. Prerequesite

- Go to https://jupyterhub.embl.de/ & enter your EMBL credentials

- Select the last option (Data Science VSCode IDE)

- In the terminal: paste the following

`git clone https://github.com/weber8thomas/embl_predoc_course.git && cd embl_predoc_course`

- Link the data into your folder

`ln -s /g/korbel/jeong/predoc_course/WGS/*bam* .`

`ln -s /g/korbel/jeong/predoc_course/WGS/spl.tsv .`

`ln -s /g/korbel/shared/datasets/refgenomes/human/GRCh38_full_analysis_set_plus_decoy_hla.fa hg38.fa`

- Activate conda environment including bcftools, samtools, delly, alfred:

`conda activate /g/korbel2/weber/miniconda3/envs/predoc_course_bioinfo_tools`

- Get a look at the following BAM file

`samtools view BM510_WT_chr10_chr13_chr22_mixture_RG.bam | less -S`

**What is the structure of the BAM file?**

## 1. Quality checking of the BAM file

Paired-end methods can be affected by a skewed insert size distribution, read-depth methods by non-uniform coverage and split-read methods suffer from high sequencing error rates that cause mis-mappings. Prior to any structural variant discovery you should therefore evaluate the quality of the data such as the percentage of mapped reads, singletons, duplicates, properly paired reads and the insert size and coverage distributions. Picard, SAMtools, FastQC and Alfred compute some of these alignment statistics as shown below.

Regarding the QC interpretation, there are some general things to watch out for such as mapping percentages below 70%, larger than 20% duplicates or multiple peaks in the insert size distribution. Be aware that many alignment statistics vary largely by protocol and hence, it's usually best to compare multiple different sequencing runs using the same protocol (DNA-seq, RNA-seq, ChIP-seq, paired-end, single-end or mate-pair) against each other, which then highlights the outliers.

```bash
alfred qc \
    -r hg38.fa  \
    -o qc_RPEWtp30_chr10_chr13_chr22.tsv.gz \
    -j qc_RPEWtp30_chr10_chr13_chr22.json.gz \
    RPEWtp30_chr10_chr13_chr22.bam
```

```bash
zcat qc_RPEWtp30_chr10_chr13_chr22.tsv.gz | grep ^ME > qc_RPEWtp30_chr10_chr13_chr22.txt
```

```bash
awk '
{
    for (i=1; i<=NF; i++) {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' qc_RPEWtp30_chr10_chr13_chr22.txt | less -S
```

1. **What is the median coverage of the data set?**
2. **What is the mapping percentage of the data set? and what will be the cutoff?**
3. **What is the duplicate rate? and what will be the cutoff?**

## 2. Structural Variant Calling

Delly is a state of the art tool used to identify structural variants from normal or tumor samples. The example provided below is to detect duplication (DUP).
You can detect several SV types (DEL, DUP, INV, TRA) using the similar command line with -t option.

### Example of duplication calling

```bash
delly call \
    --type DUP \
    --noindels \
    --map-qual 20 \
    --genome hg38.fa \
    --outfile sv_DUP.bcf \
    BM510_chr10_chr13_chr22.bam RPEWtp30_chr10_chr13_chr22.bam
```

### Check the result of duplication calling

```bash
bcftools view sv_DUP_chr10_chr13_chr22.bcf | grep -v "^##" | less -S
```

## 3. Somatic Filtering

This is to identify structural variants specific to the tumor samples compared to normal. This example is to detect duplication (DUP).

### Preparation of sample list

```bash
cat spl.tsv
```

### Example of somatic duplication filtering

```bash
delly filter \
    --type DUP \
    --pass \
    --filter somatic \
    --altaf 0.25 \
    --samples spl.tsv \
    --outfile somatic_DUP.bcf \
    sv_DUP.bcf
```

### Check the result of somatic duplication filtering

```bash
bcftools view somatic_DUP.bcf | grep "^##" | less -S
```

**Do the same thing for the other types of SVs (INV, TRA, INS, DEL).**

**What kind of somatic structural variation did you detect in BM510 compared to WT?**

## 4. Visualize SVs in circos plot

In this section, we will visualize detected somatic SV in circos plot using R.

```bash
conda activate /g/korbel2/weber/miniconda3/envs/predoc_course_r_vizu
```

```bash
Rscript rcircos_1.R
```

```bash
Rscript rcircos_2.R
```

## 5. Detection of subclonal SVs

So far, we have tried to identify clonal SVs in BM510 compared to WT. What about subclonal SVs? To see if we can detect subclonal SVs with low variant allele frequency we have prepared synthetic mixture data which includes 80 percentage of reads from WT and 20 percentage of reads from BM510. **If you repeat Delly SV call and somatic filtering on translocation (TRA), can you detect this subclonal SVs which is carried by 20 percentage of the cells?**

### Call TRA using mixture data

```bash
delly call \
    -n \
    -q 20 \
    -t TRA \
    -g hg38.fa \
    -o sv_TRA_chr10_chr13_chr22_mixture.bcf \
    BM510_WT_chr10_chr13_chr22_mixture_RG.bam RPEWtp30_chr10_chr13_chr22.bam
```

### Check the result of TRA calling

```bash
bcftools view sv_TRA_chr10_chr13_chr22_mixture.bcf | less -S
```

### Filtering somatic TRA using mixture data

```bash
delly filter \
    -t TRA \
    -p \
    -f somatic \
    -o somatic_TRA_mixture.bcf \
    -a 0.05 \
    -s spl.tsv \
    sv_TRA_chr10_chr13_chr22_mixture.bcf
```

### Check the result of somatic TRA calling

```bash
bcftools view somatic_TRA_mixture.bcf | less -S
```

```bash
for type in DUP TRA INS DEL INV
do
    delly call \
        --type "{$type}" \
        --noindels \
        --map-qual 20 \
        --genome hg38.fa \
        --outfile sv_"{$type}".bcf \
        BM510_chr10_chr13_chr22.bam RPEWtp30_chr10_chr13_chr22.bam
done
```

```bash
for type in DUP TRA INS DEL INV
do
    delly filter \
        --type "$type" \
        --pass \
        --filter somatic \
        --altaf 0.25 \
        --samples spl.tsv \
        --outfile somatic_"${type}".bcf \
        sv_"$type".bcf
done
```
