##The BAM files are here, RPEWtp30 is the wildtype:
/g/korbel2/rausch/RPECellLine/results/nRex/BM150/BM150.bam
/g/korbel2/rausch/RPECellLine/results/nRex/BM510/BM510.bam
/g/korbel2/rausch/RPECellLine/results/nRex/RPEWtp30/RPEWtp30.bam



##Preparing RPE WT bam
module load SAMtools/1.3.1-foss-2016b
samtools view -b RPEWtp30.bam chr10 chr13 chr22 > RPEWtp30_chr10_chr13_chr22.bam
samtools index RPEWtp30_chr10_chr13_chr22.bam
#samtools view -H RPEWtp30_chr10_chr22.bam > header_for_alfred.sam
#samtools view RPEWtp30_chr10_chr22.bam | cat header_for_alfred.sam - | samtools view -Sb - > RPEWtp30_chr10_chr22.reheadered.bam
#samtools index RPEWtp30_chr10_chr22.reheadered.bam
#alfred qc -r /g/korbel/shared/datasets/refgenomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o qc_RPEWtp30_chr10_chr22.tsv.gz -j qc_RPEWtp30_chr10_chr22.json.gz RPEWtp30_chr10_chr22.reheadered.bam
#zcat qc_RPEWtp30_chr10_chr22.tsv.gz | grep ^ME > qc_RPEWtp30_chr10_chr22.txt


##Preparing BM510 bam
module load SAMtools/1.3.1-foss-2016b
samtools view -b BM510.bam chr10 chr13 chr22 > BM510_chr10_chr13_chr22.bam
samtools index BM510_chr10_chr22.bam
#samtools view -H BM510_chr10_chr22.bam > header_for_alfred_BM510.sam
#samtools view BM510_chr10_chr22.bam | cat header_for_alfred_BM510.sam - | samtools view -Sb - > BM510_chr10_chr22.reheadered.bam
#samtools index BM510_chr10_chr22.reheadered.bam
#alfred qc -r /g/korbel/shared/datasets/refgenomes/human/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -o qc_BM510_chr10_chr22.tsv.gz -j qc_BM510_chr10_chr22.json.gz BM510_chr10_chr22.reheadered.bam
#zcat qc_BM510_chr10_chr22.tsv.gz | grep ^ME > qc_BM510_chr10_chr22.txt



##Preparing in-silico mixture (Let's make 20% subsampling of BM510)
module load BWA/0.7.15-foss-2016b SAMtools/1.3.1-foss-2016b
module load BEDTools/2.26.0-foss-2016b
module load GATK/3.7-Java-1.8.0_112
module load picard/2.9.0-Java-1.8.0_112
module load R/3.4.0-foss-2016b

java -jar $EBROOTPICARD/picard.jar DownsampleSam I=BM510_chr10_chr13_chr22.bam O=BM510_chr10_chr13_chr22_downsample10.bam P=0.2
java -jar $EBROOTPICARD/picard.jar DownsampleSam I=RPEWtp30_chr10_chr13_chr22.bam O=RPEWtp30_chr10_chr13_chr22_downsample.bam P=0.6




##Merge bamfiles and Unify read group ID
samtools merge BM510_WT_chr10_chr13_chr22_mixture.bam RPEWtp30_chr10_chr13_chr22_downsample.bam BM510_chr10_chr13_chr22_downsample10.bam
samtools index BM510_WT_chr10_chr13_chr22_mixture.bam

java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=BM510_WT_chr10_chr13_chr22_mixture.bam O=BM510_WT_chr10_chr13_chr22_mixture_RG.bam RGID=BM510/BM510 RGLB=BM510/BM510 RGPL=Illumina RGPU=unit1 RGSM=BM510/BM510
samtools index BM510_WT_chr10_chr13_chr22_mixture_RG.bam


##Let's check the % of reads coming from BM510 and RPEWT
-bash-4.2$ samtools view BM510_chr10_chr13_chr22_downsample10.bam | wc -l
#12279511
-bash-4.2$ samtools view RPEWtp30_chr10_chr13_chr22_downsample.bam | wc -l
#50302469

#> 12279511/(12279511+50302469)
#[1] 0.1962148

