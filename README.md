# Whole Exome Sequencing (WES) Pipeline Using GATK

## Overview
This pipeline automates the processing of Whole Exome Sequencing (WES) data using GATK4. It starts from FASTQ files, performs alignment, variant calling, and filtering, and produces a final VCF file. The final variants can be annotated using either **Annovar** or **VEP**.

## Prerequisites
Ensure the following dependencies are installed:
- **GATK 4.3.0.0**
- **BWA** (for alignment)
- **SAMtools** (for BAM file manipulation)
- **Annovar** (for variant annotation, optional)
- **VEP** (for variant annotation, optional)

## Usage

1. **Prepare the environment**
   ```bash
   chmod +x run_pipeline.sh
   ./run_pipeline.sh
   ```

2. **Pipeline steps:**
   - Converts FASTQ to BAM.
   - Marks adapters and performs alignment using BWA.
   - Sorts and merges BAM files.
   - Performs Base Quality Score Recalibration (BQSR).
   - Calls variants using GATK HaplotypeCaller.
   - Filters variants to generate a final VCF.
   - Annotates variants using **Annovar** or **VEP** (optional).

## Variant Annotation

### Using ANNOVAR
```bash
perl annovar/convert2annovar.pl --format vcf results/sample.vcf --outfile results/sample.avinput
perl annovar/table_annovar.pl results/sample.avinput annovar/humandb/ -buildver hg19 -out results/sample_annotated -remove \
-protocol refGene,1000g2015aug_all,clinvar_20200316 -operation g,f,f -nastring NA -csvout -otherinfo
```

### Using VEP
```bash
vep --input_file results/sample.vcf --output_file results/sample_annotated.vep.vcf --cache --dir_cache vep_cache/ \
--assembly GRCh37 --offline --fasta ref/ucsc.hg19.fasta --everything
```

## Notes
- **Ensure your FASTQ files follow the naming convention:** `Sample*_R1.fastq.gz` and `Sample*_R2.fastq.gz`.
- **Reference files** must be indexed properly (`.dict`, `.fai`, and sorted databases).
- Modify the **BED file** path if using a different target region.

## Output
- `results/sample_Final_Filtered_Variants.vcf` - Final high-confidence variants.

## References
- GATK Best Practices: https://gatk.broadinstitute.org
- Annovar: https://annovar.openbioinformatics.org
- VEP (Variant Effect Predictor): https://www.ensembl.org/info/docs/tools/vep/index.html
---


