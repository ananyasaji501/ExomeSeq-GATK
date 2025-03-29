#!/usr/bin/bash

## Enter Data directory
cd /mnt/d/gatk/run  ##1
REF='/mnt/d/gatk/ref' ##2  ## The reference directory should contain all the reference files including the snp databases.
GATK4='/mnt/d/gatk/gatk-4.3.0.0' ##3
BED='/mnt/d/gatk/ref/truseq-exome-targeted-regions-manifest-v1-2.bed' ##4

## Reading files in folder
for i in $(ls *_R1.fastq.gz | cut -d "_" -f 1)   ##5
do
mkdir ${i}
echo "$i"
mv ${i}_R1.fastq.gz ${i}
mv ${i}_R2.fastq.gz ${i}
done

cd /mnt/d/gatk/run   ##6

ls | while read line
do 
echo $line

cd /mnt/d/gatk/run/$line/    ##7

echo "DATE"
date

## 1. FastqToSam Start
echo "# FastqToSam Start #"
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar FastqToSam -F1 *_R1.fastq.gz -F2 *_R2.fastq.gz -O unmapped.bam -RG H0164.2 -SM ${line} -LB Solexa-272222 -PU H0164ALXX140820.2 -PL ILLUMINA --MAX_RECORDS_IN_RAM 10000000
echo "# FastqToSam End #"

## 2. Adapter Marking Start
echo "# Adapter Marking Start #"  
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar MarkIlluminaAdapters -I unmapped.bam -O unmapped_markilluminaadapters.bam -M unmapped_markilluminaadapters_metrics.txt --MAX_RECORDS_IN_RAM 10000000
echo "# Adapter Marking End #"

## 3. Sam To Fastq Start
echo "# Sam To Fastq Start #" 
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar SamToFastq -I unmapped_markilluminaadapters.bam --FASTQ unmapped_markilluminaadaptors_R1.fastq.gz --SECOND_END_FASTQ unmapped_markilluminaadaptors_R2.fastq.gz -CLIP_ATTR XT -CLIP_ACT 2 -NON_PF true --MAX_RECORDS_IN_RAM 10000000
echo "# Sam To Fastq End #"

## 4. Bwa Alignment
echo "# Bwa Alignment Start #"  
bwa mem -t 2 -M $REF/ucsc.hg19.fasta unmapped_markilluminaadaptors_R1.fastq.gz unmapped_markilluminaadaptors_R2.fastq.gz -o Alignment.sam
echo "# Bwa Alignment End #"

## 5. Sorting & BAM conversion
echo "# Sorting & BAM conversion Start #"
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar SortSam -I Alignment.sam -O Sorted_Alignment.bam -SO coordinate --CREATE_INDEX true --MAX_RECORDS_IN_RAM 10000000
echo "# Sorting & BAM conversion End #"

## 6. Bam file qc using bamqc
echo "# BAMQC Start #"  
bamqc --outdir=sample_dir SortedAligned.bam
echo "# BAMQC End #"

## 7. Merge BAM Start #Ref should have Seq dictionary (i.e Reference.dict file) and .fai index file ($ samtools faidx reference.fa)
echo "# Merge BAM Start #" 
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar MergeBamAlignment -ALIGNED Sorted_Alignment.bam -UNMAPPED unmapped.bam -O Mapped.bam -R $REF/ucsc.hg19.fasta -MC true --CLIP_ADAPTERS false --CLIP_OVERLAPPING_READS true --INCLUDE_SECONDARY_ALIGNMENTS true --MAX_INSERTIONS_OR_DELETIONS -1 --PRIMARY_ALIGNMENT_STRATEGY MostDistant --ATTRIBUTES_TO_RETAIN XS --MAX_RECORDS_IN_RAM 10000000
echo "# Merge BAM End #"

## 9. Validate Sam
echo "# Validate Sam Start #"
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar ValidateSamFile -I Mapped.bam --REFERENCE_SEQUENCE $REF/ucsc.hg19.fasta -M SUMMARY -O SamValidation_metrics.txt --MAX_RECORDS_IN_RAM 10000000
echo "# Validate Sam End #"

## 10. MArkDuplicate
echo "# MarkDuplicate Start #"
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar MarkDuplicates -I Mapped.bam -O Markdup.bam -M Markduplicate_metrics.txt --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 --CREATE_INDEX true --CREATE_MD5_FILE true --MAX_RECORDS_IN_RAM 10000000
echo "# MarkDuplicate End #"

## 12. BASE RECALIBERATION
## Base Quality Score Recalibration (BQSR) 1
echo "# BQSR-1 Start #" 
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar BaseRecalibrator -R $REF/ucsc.hg19.fasta -I Markdup.bam -L $BED -ip 100 --known-sites $REF/dbsnp_138.hg19_sorted.vcf --known-sites $REF/Mills_and_1000G_gold_standard.indels.hg19.sites.sorted.vcf -O recal_data.table
echo "# BQSR-1 End #"

## 13. Apply Recalibration 2
echo "# Apply Recalibration-2 Start #"  
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar ApplyBQSR -bqsr recal_data.table -R $REF/ucsc.hg19.fasta -I Markdup.bam -O recal_reads.bam --create-output-bam-index true
echo "# Apply Recalibration-2 End #"

## 14. Variant Calling
echo "# Haplotypecaller Start #" 
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar HaplotypeCaller  -R $REF/ucsc.hg19.fasta -I recal_reads.bam -L $BED -ip 100 -stand-call-conf 30 -O ${line}_Raw_variants.vcf
echo "# Haplotypecaller End #"
# -ERC GVCF - Use this parameter for joint genotyping samples

## 15. Extract SNP and Idels
echo "# SNP & Indel Extraction Start #"
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar SelectVariants -R $REF/ucsc.hg19.fasta -V ${line}_Raw_variants.vcf -select-type SNP -O Raw_snps.vcf

java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar SelectVariants -R $REF/ucsc.hg19.fasta -V ${line}_Raw_variants.vcf -select-type INDEL -O Raw_indels.vcf
echo "# SNP & Indel Extraction End #"

## 16.1. SNP Filtering
echo "# SNP Filtering Start #"
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar VariantFiltration -R $REF/ucsc.hg19.fasta -V Raw_snps.vcf --filter-expression "QD < 2.0" --filter-name "QualityByDepth2" --filter-expression "FS > 60.0" --filter-name "FisherStrand60" --filter-expression "MQ < 40.0" --filter-name "RMSMappingQuality40" --filter-expression "MQRankSum < -12.5" --filter-name "MappingQualityRankSumTest-12.5" --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumTest-8" -O snps_filtered.vcf
echo "# SNP Filtering End #"

## 16.2. Indel Filtering
echo "# Indel Filtering Start #"  
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar VariantFiltration -R $REF/ucsc.hg19.fasta -V Raw_indels.vcf --filter-expression "QD < 2.0" --filter-name "QualityByDepth2" --filter-expression "FS > 200.0" --filter-name "FisherStrand200" --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O indels_filtered.vcf
echo "# Indel Filtering End #"

## 17. Combine variants
echo "# Combine Variants Start #" 
java -Xmx6G -jar $GATK4/gatk-package-4.3.0.0-local.jar MergeVcfs -I snps_filtered.vcf -I indels_filtered.vcf -O ${line}_Final_Filtered_Variants.vcf
echo "# Combine Variants End #"

# deletion of intermediate files generated
echo "# File Deletion Start #" 
rm unmapped_markilluminaadaptors_R1.fastq.gz
rm unmapped_markilluminaadaptors_R2.fastq.gz
rm Alignment.sam
rm Markdup.bai
rm Markdup.bam
rm recal_data.table
rm Sorted_Alignment.bam
rm unmapped.bam
rm unmapped_markilluminaadapters.bam
rm Raw_snps.vcf
rm Raw_indels.vcf
echo "## File Deletion End ##"

echo "DATE"
date

done

### Annotate the vcf file Generated
## 18. Annovar Annotation
  ## Convert the final filtered VCF file to Annovar input format
  ## and perform functional annotation using Annovar databases.
  ## Uncomment the below lines to enable Annovar annotation.

  # echo "# Annovar Annotation Start #"
  # perl /data/applications/annovar/convert2annovar.pl --format vcf4 ${line}_Final_Filtered_Variants.vcf --includeinfo --outfile ${line}.avinput --withzyg
  # perl /data/applications/annovar/table_annovar.pl ${line}.avinput /data/applications/annovar/humandb/ -buildver hg19 -out ${line} -remove -protocol refGene,dbnsfp31a_interpro,clinvar_20200316,avsnp150,1000g2015aug_all,gnomad211_exome -operation g,f,f,f,f,f -nastring NA -csvout -otherinfo
  # echo "# Annovar Annotation End #"

  ## 19. VEP Annotation
  ## Use VEP (Variant Effect Predictor) to annotate the final VCF file
  ## with gene symbols, SIFT, PolyPhen, and other variant effect predictions.
  ## Uncomment the below lines to enable VEP annotation.

  # echo "# VEP Annotation Start #"
  # vep --input_file ${line}_Final_Filtered_Variants.vcf --output_file ${line}_vep_annotated.vcf --format vcf --vcf --symbol --sift b --polyphen b --everything --custom /data/vep_plugins/custom_database.vcf.gz,CustomDB,vcf,exact,0 --dir_cache /data/vep_cache/
  # echo "# VEP Annotation End #"
