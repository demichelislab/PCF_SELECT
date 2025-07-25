# PCF\_SELECT - Singularity package

[![DOI](https://zenodo.org/badge/DOI/10.1093/narcan/zcac016.svg)](https://doi.org/10.1093/narcan/zcac016)

- [Introduction](#intro)
- [Execute the PCF\_SELECT pipeline](#execute-pcf)
  * [Declare the environment variables](#declare-the-environment-variables)
  * [Run the whole PCF\_SELECT pipeline](#run-the-whole-pcfs-pipeline)
    + [Sample Information File (SIF)](#sample-information-file--sif-)
  * [List apps present in the singularity package](#list-apps-present-in-the-singularity-package)
  * [Run apps other than PCFS](#run-apps-other-than-pcfs)
- [Acknowledgments](#acknow)

## Introduction <a name="intro"></a>

This repository provides the **Singularity packages** to run the PCF\_SELECT pipeline described in:  
**[Allele-informed copy number evaluation of plasma DNA samples from metastatic prostate cancer patients: the PCF_SELECT consortium assay - Orlando et al., NAR Cancer (2022)](https://doi.org/10.1093/narcan/zcac016)**  

Depending on the panel version, a specific Singularity package must be used (available under [`singularity/versions/`](singularity/versions/)).  

## Execute the PCF\_SELECT pipeline <a name="execute-pcf"></a>

### Declare the environment variables <a name="declare-the-environment-variables"></a>
In case your HPC do have issues or special requirements in managing the temporary space (/tmp folder), it's better to declare a set of environment variables that would allow you to choose a custom folder for storing temporary files and folders:
```
export SINGULARITY_BINDPATH=<Folder you want to bind to the singularity image>
export SINGULARITY_TMPDIR=<Temporary folder>
export SINGULARITY_CACHEDIR=<Temporary folder>
export SINGULARITY_LOCALCACHEDIR=<Temporary folder>
export SINGULARITY_CLEANENV=1
```
### Run the whole PCF\_SELECT pipeline <a name="run-the-whole-pcfs-pipeline"></a>
```
singularity run --app PCFS pcfselect.sif -h
usage: pcfs.py [-h] -s SIF -o OUTDIR -t TMPDIR [-n NCORES]

Arguments for running PCF-SELECT

optional arguments:
  -h, --help            show this help message and exit
  -s SIF, --sif SIF     Sample Info Description
  -o OUTDIR, --outDir OUTDIR
                        Directory where to store results
  -t TMPDIR, --tmpDir TMPDIR
                        Temporary directory
  -n NCORES, --nCores NCORES
                        N cores to be allocated for the analysis
```
#### Sample Information File (SIF) <a name="sample-information-file--sif-"></a>
The Sample Information File (SIF) is a TSV (tab-separated) text file that contains the location of Tumour samples to analyse, together with their Normal counterpart.
It's **crucial** to name the columns "Tumour" and "Normal" in order to recognise the right columns to analyse. Samples are intended in [BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf). 

| Patient | Tumour | Normal |
| ------- | ------ | ------ |
| Patient_1 | BAM 1T | BAM 1N |
| Patient_2 |BAM 2T | BAM 2N |
| .... | ....  | ....   |

### List apps present in the singularity package <a name="list-apps-present-in-the-singularity-package"></a>
```
singularity inspect --list-apps pcfselect.sif

ABEMUS
ABEMUS_2_ANNOVAR
ALLELIC_IMBALANCE
ANNOVAR_2_ONCOTATOR
CN_SNVs_CALLS
CORRECT_VAF
GENERATE_FOCAL_TABLES
PCFS
PEAK_CORRECTION
SEGMENTATION_BETA_COMPUTATION
SIF_2_ABEMUS
TC_ESTIMATION
```

### Run apps other than PCFS <a name="run-apps-other-than-pcfs"></a>
In this example we list the arguments to run the app "ABEMUS_2_ANNOVAR"
```
run --app ABEMUS_2_ANNOVAR pcf_select_latest.sif -h
usage: /usr/local/bin/pcfs/scripts/ABEMUS_to_annovar.R [-h] -s ABEMUSFILE -o OUTDIR -n NCORES

optional arguments:
  -h, --help            show this help message and exit
  -s ABEMUSFILE, --abemusFile ABEMUSFILE
                        File ABEMUS containing the SNPs to be converted
  -o OUTDIR, --outDir OUTDIR
                        Output directory
  -n NCORES, --nCores NCORES
                        Number of available cores for the analysis
```

## Acknowledgments <a name="acknow"></a>
| Author | Institution | Contribution |
| --------- | ----------- | ------------ |
| Francesco Orlando, PhD | CIBIO, University of Trento | Main author of PCF\_SELECT |
| Stefano Pirrò, PhD | Cancer Institute, UCL | Porting of PCF\_SELECT on singularity |
| Osvaldas Vainauskas | Cancer Institute, UCL | Contributor |
| Prof. Francesca Demichelis | CIBIO, University of Trento | Principal Investigator |
| Prof. Gerhardt Attard | Cancer Institute, UCL | Principal Investigator |