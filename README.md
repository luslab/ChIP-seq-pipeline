# ChIP-seq pipeline

This repository contains scripts for the analysis of the Transcription Factor ChIP-seq (TF ChIP-seq) data ([GEO: GSE59530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59530)) from [Franco HL, _et al._, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25752574). For the sake of time, we will focus only on the Input and ChIP-seq _FoxA1 Vehicle_ samples without treatments.

The aim is to create a TF ChIP-seq pipeline that includes quality controls, data analysis and visualisation. The pipeline will be developed using [Snakemake](http://snakemake.readthedocs.io/en/stable/), a tool to create reproducible and scalable data analyses.

## Part 1

### Download and process raw data

The Sequence Read Archives (SRA) files required for this analysis can be retrieved from the the [National Center for Biotechnology Information (NCBI)](https://www.ncbi.nlm.nih.gov/).

__Requirements:__
* `curl`
* `snakemake` from _Snakemake_
* `fastq-dump` from [_SRA Toolkit_](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)

Since it is a bit tricky to provide remote files as input in Snakemake, we will download the data externally using `curl`.

```
# Create the output folder
mkdir -p data

# Cycle through the identifiers in 'SraAccList.txt' and download the corresponding file
for sra in `awk '{print $1}' SraAccList.txt`; do
    curl "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${sra:0:6}/$sra/$sra.sra" --output data/$sra
done
```

Have a look at the _Snakemake_ file in the main directory. I have already created some _rules_ to convert the _sra_ data to _fastq_ format. You can run this first step of the pipeline by simply typing `snakemake`.

```
snakemake
```

### Assess reads quality

__Requirements:__
* `fastqc` from [_FastQC_](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### Pre-processing

__Requirements:__
* `cutadapt` from [_Cutadapt!_](https://cutadapt.readthedocs.io/en/stable/)
* `trim_galore` from [_Trim Galore!_](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

### Map to the genome

__Requirements:__
* `bowtie2-build` and `bowtie2` from [_Bowtie2_](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

### Post-processing

Additional quality controls including _Nonredundant Fraction_ (NRF = Unique start positions of uniquely mappable reads / Uniquely mappable reads) and Cross-correlation (Strand cross-correlation is computed as the Pearson correlation between the positive and the negative strand profiles at different strand shift distances, k)

## Part 2

### Peak calling

__Requirements:__
* `macs2` from [_MACS2_](https://github.com/taoliu/MACS)

### Consistency of replicates

Calculate the Irreproducible Discovery Rate (IDR)

### Visualize the results