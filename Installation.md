# FluViewer-nf

Nextflow pipeline to run the FluViewer tool

## Analyses

- Read trimming & QC: `fastp`
- Sequence analysis with FluViewer

## To add in

- ? `MultiQC` process
- ? `cutadapt` (or do through `fastp`)
- Aggregate results in a linelist

## Usage

To come

## Output

To come

## Provencance files

Details to come.






# Standard Operating Procedure to Install FluViewer_conda

** Note this is undergoing changes so this is slightly out of date and will be updated shortly **



- This will install Kevin Kuchinski's tool FluViewer and all its necessary dependencies in a conda environment known as FluViewer.

- The dependencies are apparently fragile, so use only the yml file to build the environment.

- Also note, the included Spades version in the yml file does not have a MacOS version available on Conda, so it does not install on a Mac.

- When running FluViewer, the output directory needs to be in the same directory that you are running the script from (no subdirectories)


## Install instructions

1. Log-on to Linux server (e.g. Sabin or Almeida)

2. Enter the following on the command line:
```
git clone https://github.com/JamesZlosnik/FluViewer_installation
```

3. This will create a directory called: FluViewer_installation.

4. Go into this directory
```
cd FluViewer_installation
```

5. Make the conda environment:

```
conda env create -f FluViewer_conda.yml
```

6.  Once installed then you should be good to go once you activate the environment, e.g.
```
conda activate FluViewer
```
