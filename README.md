# PCF-SELECT - Porting to Singularity

| Author | Institution | Contribution |
| --------- | ----------- | ------------ |
| Francesco Orlando | CIBIO, University of Trento | Main author of PCF-SELECT |
| Stefano Pirrò, PhD | Cancer Institute, UCL | Porting of PCF-SELECT on singularity |
| Prof. Francesca Demichelis | CIBIO, University of Trento | Principal Investigator |
| Prof. Gerhardt Attard | Cancer Institute, UCL | Principal Investigator |

- [Introduction to Singularity](#introduction-to-singularity)
  * [Installing Singularity on the host system](#installing-singularity-on-the-host-system)
- [Build PCF-SELECT singularity from Scratch](#build-pcf-select-singularity-from-scratch)
  * [%environment](#-environment)
  * [%files](#-files)
  * [%post](#-post)
    + [Linux essential packages](#linux-essential-packages)
    + [Python and related packages](#python-and-related-packages)
    + [R and related packages](#r-and-related-packages)
    + [ABEMUS and pacBAM](#abemus-and-pacbam)
  * [%apprun](#-apprun)
- [Execute an analysis](#execute-an-analysis)
  * [Declare the environment variables](#declare-the-environment-variables)
  * [Run the whole PCFS pipeline](#run-the-whole-pcfs-pipeline)
    + [Sample Information File (SIF)](#sample-information-file--sif-)
  * [List apps present in the singularity package](#list-apps-present-in-the-singularity-package)
  * [Run apps other than PCFS](#run-apps-other-than-pcfs)

## Introduction to Singularity <a name="introduction-to-singularity"></a>
Singularity is the most widely used container system for HPC. It is designed to execute applications at bare-metal performance while being secure, portable, and 100% reproducible. Singularity is an open-source project with a friendly community of developers and users. The user base continues to expand, with Singularity now used across industry and academia in many areas of work.\
\
The main features of Singularity are:
* **Trust and Security**: Singularity is the only container system that supports public/private key signing, providing trust and guarantees of immutability.
* **Compatibility**: Singularity is 100% Docker and OCI (Open Containers Initiative) compatible (but easier to use).
* **Encrypted**: Singularity can encrypt containers and integrates with Vault and other secret management platforms to secure applications, models, and data.
* **Absolute Portability**: The single-file SIF container format allows you to reproducibly build, share, and archive your workload from workstations, to HPC, to the edge.
* **Secure**: Singularity runs "rootless" and prohibits privilege escalation within the container; users are the same inside and outside the container.

Please find more information about Singularity at the [following link](https://singularity.hpcng.org)

### Installing Singularity on the host system <a name="installing-singularity-on-the-host-system"></a>
Please refer to [this guide](https://sylabs.io/guides/3.8/user-guide/quick_start.html#quick-installation-steps) on how to download and install Singularity on your system. This is a pre-build software in most HPC distributions.

## Build PCF-SELECT singularity from Scratch <a name="build-pcf-select-singularity-from-scratch"></a>
As for all Singularity packages, PCF-SELECT can be built from scratch by using a dedicated recipe. A Singularity Recipe is the driver of a custom build, and the starting point for designing any custom container. It includes specifics about installation software, environment variables, files to add, and container metadata.\

A Singularity Recipe file is divided into several parts:
* **Header**: The Header describes the core operating system to build within the container. Here you will configure the base operating system features that you need within your container. Examples of this include, what distribution of Linux, what version, what packages must be part of a core install.
* **Sections**: The rest of the definition is comprised of sections, sometimes called scriptlets or blobs of data. Each section is defined by a % character followed by the name of the particular section. All sections are optional. Sections that are executed at build time are executed with the /bin/sh interpreter and can accept bin/sh options. Similarly, sections that produce scripts to be executed at runtime can accept options intended for /bin/sh.

**Please note!** You need ROOT permissions to build a singularity recipe.

The current PCF-SELECT is based on **Ubuntu - rolling** version and so the header will be:\
```
Bootstrap: docker
From: ubuntu:rolling
```

### %environment <a name="-environment"></a>
As of Singularity 2.3, you can add environment variables to your Singularity Recipe in a section called %environment. Keep in mind that these environment variables are sourced at runtime and not at build time.\
Here we declare the path of pacbam and annovar in order to call the command just with their aliases.
```
export PATH=/usr/local/bin/pacbam:$PATH
export PATH=/usr/local/bin/annovar:$PATH
export DEBIAN_FRONTEND=noninteractive
```

### %files <a name="-files"></a>
If you want to copy files from your host system into the container, you should do so using the %files section. Each line is a pair of <source> and <destination>, where the source is a path on your host system, and the destination is a path in the container. \
In PCF-SELECT recipe, we copy the annovar program and (of course) all the scripts/src files necessary to run the analysis correctly.
```
annovar /usr/local/bin/
pcfs /usr/local/bin/
```

### %post <a name="-post"></a>
Commands in the %post section are executed within the container after the base OS has been installed at build time. This is where the meat of your setup will live, including making directories, and installing software and libraries.
#### Linux essential packages <a name="linux-essential-packages"></a>
```
apt-get -y update
apt-get -y upgrade
apt-get -y install libncurses5-dev liblzma-dev libssl-dev libcurl4-openssl-dev
apt-get -y install git libxml2-dev pandoc
apt-get -y install locales wget
locale-gen en_US.UTF-8
```
#### Python and related packages <a name="python-and-related-packages"></a>
```
apt-get -y install python3 python3-pip
pip3 install argparse glob2 pandas
pip3 install regex datetime tenacity
```
#### R and related packages <a name="r-and-related-packages"></a>
```
apt-get -y install r-base
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
Rscript -e "install.packages('colorRamps')"
Rscript -e "install.packages('stringr')"
Rscript -e "install.packages('data.table')"
Rscript -e "install.packages('BiocManager')"
Rscript -e "install.packages('argparse')"
Rscript -e "BiocManager::install('GenomicRanges', update = TRUE, ask = FALSE)"
Rscript -e "BiocManager::install('GenomeInfoDb', update = TRUE, ask = FALSE)"
Rscript -e "BiocManager::install('IRanges', update = TRUE, ask = FALSE)"
Rscript -e "BiocManager::install('S4Vectors', update = TRUE, ask = FALSE)"
Rscript -e "BiocManager::install('BiocGenerics', update = TRUE, ask = FALSE)"
Rscript -e "install.packages('funr')"
Rscript -e "install.packages('TPES')"
Rscript -e "install.packages('CLONETv2')"
Rscript -e "install.packages('devtools')"
Rscript -e "install.packages('rmarkdown')"
Rscript -e "install.packages('pacman')"
Rscript -e "install.packages('dplyr')"
Rscript -e "install.packages('fitdistrplus')"
Rscript -e "install.packages('stringr')"
Rscript -e "install.packages('tidyr')"
Rscript -e "install.packages('doFuture')"
Rscript -e "install.packages('doRNG')"
Rscript -e "install.packages('foreach')"
Rscript -e "install.packages('parallel')"
```
#### ABEMUS and pacBAM <a name="abemus-and-pacbam"></a>
```
Rscript -e "devtools::install_github('cibiobcg/abemus', build_vignettes = T)"
git clone 'https://CibioBCG@bitbucket.org/CibioBCG/pacbam.git'
cd pacbam
make -f Makefile.linux
mv pacbam /usr/local/bin/ && cd
```

### %apprun <a name="-apprun"></a>
Starting in Singularity 2.4, multiple commands can be used in the context of internal modules called **apps** based on the [Standard Container Integration Format](https://sci-f.github.io/).
**PCFS** app refers to the execution of the whole pipeline of analysis, while all other commands refers to each step, taken singularly.
```
%apprun PCFS
    exec python3 /usr/local/bin/pcfs/pcfs.py $@

%apprun SIF_2_ABEMUS
    exec Rscript /usr/local/bin/pcfs/scripts/SIF_2_ABEMUS.R $@

%apprun ABEMUS
    exec Rscript /usr/local/bin/pcfs/scripts/ABEMUS.R $@

%apprun ABEMUS_2_ANNOVAR
    exec Rscript /usr/local/bin/pcfs/scripts/ABEMUS_to_annovar.R $@

%apprun ANNOVAR_2_ONCOTATOR
    exec Rscript /usr/local/bin/pcfs/scripts/annovar_to_oncotator.R $@

%apprun PEAK_CORRECTION
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/peak_correction.R $@

%apprun SEGMENTATION_BETA_COMPUTATION
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/segmentation_betaComputation.R $@

%apprun TC_ESTIMATION
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/TC_estimation.R $@

%apprun GENERATE_FOCAL_TABLES
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/generate_focal_tables.R $@

%apprun ALLELIC_IMBALANCE
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/generate_allelic_imbalance_log2_tables.R $@

%apprun CN_SNVs_CALLS
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/make_CN_SNVs_calls.R $@

%apprun CORRECT_VAF
    exec Rscript /usr/local/bin/pcfs/scripts/pcf-select/correct_VAF.R $@
```

## Execute an analysis <a name="execute-an-analysis"></a>

### Declare the environment variables <a name="declare-the-environment-variables"></a>
In case your HPC do have issues or special requirements in managing the temporary space (/tmp folder), it's better to declare a set of environment variables that would allow you to choose a custom folder for storing temporary files and folders:
```
export SINGULARITY_BINDPATH=<Folder you want to bind to the singularity image>
export SINGULARITY_TMPDIR=<Temporary folder>
export SINGULARITY_CACHEDIR=<Temporary folder>
export SINGULARITY_LOCALCACHEDIR=<Temporary folder>
export SINGULARITY_CLEANENV=1
```
### Run the whole PCFS pipeline <a name="run-the-whole-pcfs-pipeline"></a>
```
singularity run --app PCFS /cluster/scratch9/spirro_scratch9/pcfs/pcf_select_latest.sif -h
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
It's **crucial** to name the columns "Tumour" and "Normal" in order to recognise the right columns to analyse.

| Patient | Tumour | Normal |
| ------- | ------ | ------ |
| Patient_1 | BAM 1T | BAM 1N |
| Patient_2 |BAM 2T | BAM 2N |
| .... | ....  | ....   |

### List apps present in the singularity package <a name="list-apps-present-in-the-singularity-package"></a>
```
singularity inspect --list-apps pcf_select_latest.sif

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
