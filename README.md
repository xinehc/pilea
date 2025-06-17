# Pilea
Pilea: profiling bacterial growth dynamics from metagenomes with sketching

## Quick Start
### Installation
[![Conda Version](https://anaconda.org/bioconda/pilea/badges/version.svg)](https://anaconda.org/bioconda/pilea)
[![Latest Release](https://anaconda.org/bioconda/pilea/badges/latest_release_date.svg)](https://anaconda.org/bioconda/pilea)

Install Pilea in a new conda environment:

```bash
conda create -n pilea -c bioconda -c conda-forge pilea
conda activate pilea
```

### Executing Pilea
Download the [example dataset](https://doi.org/10.5281/zenodo.15681352) from Zenodo. This dataset contains 1 million paired-end short reads from four bacterial species (*Bacillus subtilis*, *Klebsiella pneumoniae*, *Morganella morganii*, and *Pseudomonas\_E putida\_U*):

```bash
wget -qN --show-progress https://zenodo.org/records/15681353/files/example.tar.gz
tar -xvf example.tar.gz && cd example
```

<ins>***Step 1: Indexing***</ins>

> [!TIP]
> If MAGs aren't available, use `pilea fetch` to download a pre-built [GTDB database](https://doi.org/10.5281/zenodo.15596115). This reference database was constructed with `pilea rebuild`.

The taxonomy mapping file (`-a/--taxonomy`) is optional. If provided, it needs to be a tab-separated file containing at least two columns (genome and its associated taxonomy). This file can be the output of GTDB-Tk (`gtdbtk.bac120.summary.tsv`) and should include only bacteria (no archaea or non-prokaryotes). MAGs must have extensions in `.(fa|fna|fasta)`.

```bash
pilea index mags/*.fna -a gtdbtk.bac120.summary.tsv -o db
```

<ins>***Step 2: Profiling***</ins>

> [!TIP]
> If multiple samples are available, running them in a single batch (`*.fasta`) can help avoid repeatedly loading the reference database, which can be time-consuming if the database is large.

Both FASTA (`fa|fasta`) and FASTQ (`fq|fastq`) files are supported. By default, paired-end reads are identified with pattern `_(1|2|R1|R2|fwd|rev)`. Use `--single` for single-end reads. PTR estimates and other metadata for MAGs that pass basic filters will be saved to `output.tsv`.

```bash
pilea profile *.fasta -d db -o .
```

