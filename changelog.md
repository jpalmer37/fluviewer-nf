# Changelog

All notable updates to this pipeline will be tracked here.
Follows: https://keepachangelog.com/en/1.0.0/ recommended by nf-core

## [Unreleased]

## Version 0.1.0 - 2022-10-06

### Added

- Cutadapt step immediately after fastp to remove primers
- FastQC step after cutadapt for final assesment of reads
- Run both align and assemble mode
- Provenance outputs, outputted into own directory under 
- Changelog
- Ability to add FluViewer options explicitly

### Changed
- Update ReadMe file including, guidance on running, outputs and diagram.
- Updated `nextflow.config` to always output a timeline and report for each analysis.
- Updated outputs of `align` mode to have `_align` as a suffix and `assemble` to have `assemble` as a suffix.
- Brought cutadapt outputs into multiqc.
- Appended sampleID with either `_align` or `_assemble` to make clear which output is which. 
- Changed output older locations so that folders are: Runname/FluViewer_version/sample/[align or assemble]/outputs.


### Removed
- Qualimap process


