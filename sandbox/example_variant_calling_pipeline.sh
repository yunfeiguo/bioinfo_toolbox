#$ -cwd
#$ -V
#$ -M guoyunfei1989@gmail.com
#$ -m ea
#$ -l h_vmem=1G
#$ -pe smp 12

#generate SA coordinates for alignment (note: -I denotes illumina 1.3+ format)
bwa aln -t 12 -q 10 -I -f test.sai /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta s4.1.fastq
bwa aln -t 12 -q 10 -I -f test.fq2.sai /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta s4.2.fastq
#align the FASTQ files according to SA coordinates
bwa sampe -r '@RG\tID:READGROUP\tSM:SAMPLE' -f test.sam /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta test.sai test.fq2.sai s4.1.fastq s4.2.fastq

#convert SAM file to BAM file to improve performance later on
samtools view -b -S -o test.bam test.sam

#sort the cleaned BAM files (Picard does not recognize SamTools-sorted BAM files)
java -Xmx4g -jar /home/kaiwang/usr/picard/picard-tools/SortSam.jar INPUT=test.bam OUTPUT=test.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT

#remove duplicates in the cleaned sorted BAM files
java -Xmx4g -jar /home/kaiwang/usr/picard/picard-tools/MarkDuplicates.jar INPUT=test.sort.bam OUTPUT=test.rmdup.bam METRICS_FILE=test.rmdup.metrics REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT

#generate BAM index for rmdupped BAM files
samtools index test.rmdup.bam

#generate GATK intervals
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T RealignerTargetCreator  -I test.rmdup.bam -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta -o test.intervals

#realign BAM files based on GATK intervals
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T IndelRealigner  -I test.rmdup.bam -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta -targetIntervals test.intervals --out test.realign.bam

#index realigned BAM files
samtools index test.realign.bam

#count covariates for base quality recalibration
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T CountCovariates  -I test.realign.bam -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --default_platform illumina -knownSites /home/kaiwang/project/seqlib/gatk/dbsnp_132.b37.vcf -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate -recalFile test.recal_data.csv

#plot recalibration performance before recalibration
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/AnalyzeCovariates.jar  -ignoreQ 5 -recalFile test.recal_data.csv -outputDir test.before_recal/ -resources /home/kaiwang/usr/gatk/GenomeAnalysisTK/resources/

#recalibrate base score
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T TableRecalibration  -I test.realign.bam -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --default_platform illumina -o test.recal.bam -recalFile test.recal_data.csv

#index recalibrated BAM files
samtools index test.recal.bam

#count covariates for recalibration performance plot
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T CountCovariates  -I test.recal.bam -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --default_platform illumina -knownSites /home/kaiwang/project/seqlib/gatk/dbsnp_132.b37.vcf -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate -recalFile test.2.recal_data.csv

#plot recalibration performance after recalibration
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/AnalyzeCovariates.jar  -ignoreQ 5 -recalFile test.2.recal_data.csv -outputDir test.after_recal/ -resources /home/kaiwang/usr/gatk/GenomeAnalysisTK/resources/

#generate GATK SNP and indel calls
java -Xmx4g -jar /home//usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T UnifiedGenotyper  -l INFO -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta -I test.recal.bam -o test.raw.vcf --genotype_likelihoods_model BOTH -nt 12

#separate variant calls to a SNP file and an indel file
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T SelectVariants  -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --variant test.raw.vcf -o test.raw.indel.vcf -selectType INDEL
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T SelectVariants  -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --variant test.raw.vcf -o test.raw.snp.vcf -selectType SNP

#apply variant recalibration to assign confidence values to SNP calls
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T VariantRecalibrator  -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta -input test.raw.snp.vcf \
--maxGaussians 4 --percentBadVariants 0.05 \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/kaiwang/project/seqlib/gatk/hapmap_3.3.b37.sites.vcf \
-resource:omni,known=false,training=true,truth=false,prior=12.0 /home/kaiwang/project/seqlib/gatk/1000G_omni2.5.b37.sites.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 /home/kaiwang/project/seqlib/gatk/dbsnp_132.b37.vcf \
-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
-mode SNP -recalFile test.recal_snp -tranchesFile test.tranches.metrics -rscriptFile test.plots.R -nt 12

java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T ApplyRecalibration  -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta -input test.raw.snp.vcf --ts_filter_level 99.0 -tranchesFile test.tranches.metrics -recalFile test.recal_snp -o test.filter.snp.vcf

#filter high-quality indel calls by hard criteria
java -Xmx4g -jar /home//usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T VariantFiltration  -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --variant test.raw.indel.vcf -o test.filter.indel.vcf --clusterWindowSize 10 --filterExpression "QD < 2.0" --filterName QDFilter --filterExpression "ReadPosRankSum < -20.0" --filterName ReadPosFilter --filterExpression "FS > 200.0" --filterName FSFilter

#combine SNP and indel calls
java -Xmx4g -jar /home/kaiwang/usr/gatk/GenomeAnalysisTK/GenomeAnalysisTK.jar -T CombineVariants  -R /home/kaiwang/project/seqlib/g1k_v37/human_g1k_v37.fasta --variant test.filter.snp.vcf --variant test.filter.indel.vcf -o test.filter.vcf

#generate ANNOVAR input files
/home/kaiwang/usr/kgenome/trunk/bin/convert2annovar.pl -format vcf4 test.filter.snp.vcf -coverage 10 -filter PASS > test.avinput
/home/kaiwang/usr/kgenome/trunk/bin/convert2annovar.pl -format vcf4 test.filter.indel.vcf -coverage 10 -filter PASS >> test.avinput

#Two important modifications: 
#change '~' to full path such that I can run commands locally
#name all metrics files with suffix 'metrics'
