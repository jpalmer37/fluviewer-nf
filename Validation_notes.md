# Validation notes

To keep track of validation steps for FluViewer

## Overview

I have put the FluViewer pipeline into a nextflow wrapper, with:
- Fastp upstream to clean up reads
  - In addition, I have supplied a list of primers/adapters to add to the trimming
  - Fastp outputs a collated into a MultiQC report


## Initial trials

- All initial runs were done on the samples in test_data
- 2022-04-22 added gblock primer to primers_adapters.fa to account for presence in contemporary samples in the validation dataset where this had been used.  Test pipeline with this and seemed to work ok.


## Datasets

- Kevin's avian influenza, these were copied into `AIV_egg_culture`
- Samples sequenced by the lab.
- These are copied into:
  -

Command run
```
nextflow run FluViewer_installation/main.nf --fastq_input test_data_files/ --ref ref/FluViewer_db.fa -profile conda --cache ~/.conda/envs --outdir output_test -with-dag flowchart.html
```

**Note that:**
  - the flowchart is just for ensuring the pipeline looks correct and can be dropped.
  - `-profile` and `--cache` are essential for making sure that the conda environment works properly
  - `--fastq_input` is the location of the fastq files for this.
  - `-ref`
  - There are other options available (e.g. --adapters for the location of a file to add to the fastp trimming) - see the nextflow.config file for details.


## Running the pipeline

1. Kevin's Datasets

2. 2020_samples from the lab

3. Contempoary samples (also from the lab)

4. Avian sample run 4 times

5. 2019_samples from the lab

Commands run:


2020_samples from the lab.

Align mode:

```
nextflow run FluViewer_installation/main.nf --fastq_input validation/2020_samples/ --ref ref/FluViewer_db.fa -profile conda --cache ~/.conda/envs --outdir 2020_samples_align
```
Assemble mode:

```
nextflow run FluViewer_installation/main.nf --fastq_input validation/2020_samples/ --ref ref/FluViewer_db.fa -profile conda --cache ~/.conda/envs --outdir 2020_samples_assemble --mode "assemble"
```

Contemporary samples (also from the lab)

Align mode:


```
nextflow run FluViewer_installation/main.nf --fastq_input validation/contemporary_samples/ --ref ref/FluViewer_db.fa -profile conda --cache ~/.conda/envs --outdir contemporary_align
```


Assemble mode:

```
nextflow run FluViewer_installation/main.nf --fastq_input validation/contemporary_samples/ --ref ref/FluViewer_db.fa -profile conda --cache ~/.conda/envs --outdir contemporary_assemble --mode "assemble"

```
