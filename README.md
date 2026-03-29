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

### Running Pilea
Download the [example dataset](https://doi.org/10.5281/zenodo.15681352) from Zenodo. This dataset contains 1 million paired-end short reads from four bacterial species (*Bacillus subtilis*, *Klebsiella pneumoniae*, *Morganella morganii*, and *Pseudomonas putida*):

```bash
wget -qN --show-progress https://zenodo.org/records/15681353/files/example.tar.gz
tar -xvf example.tar.gz && cd example
```

<ins>***Step 1: Indexing***</ins>

> [!TIP]
> If MAGs aren't available, use `pilea fetch` to download a pre-built [GTDB database](https://doi.org/10.5281/zenodo.18061130). This reference database was constructed with `pilea rebuild`. If MAGs are available, they can be merged with the GTDB database via `-d/--database`. See `pilea index -h` for details.

The taxonomy mapping file (`-a/--taxonomy`) is optional. If provided, it needs to be a tab-separated file containing at least two columns (genome and its associated taxonomy). This file can be the output of GTDB-Tk (`gtdbtk.bac120.summary.tsv`) and should include only bacteria (no archaea or non-prokaryotes). MAGs must have extensions in `.(fa|fna|fasta)`.

```bash
pilea index mags/*.fna -a gtdbtk.bac120.summary.tsv -o db
```

<ins>***Step 2: Profiling***</ins>

> [!TIP]
> If multiple samples are available, running them in a single batch (`*.fasta`) helps avoid repeated database loading, which can be time-consuming if the database is large.

Both FASTA (`fa|fasta`) and FASTQ (`fq|fastq`) files are supported. By default, paired-end reads are identified with pattern `_(1|2|R1|R2|fwd|rev)`. Use `--single` for single-end reads.

```bash
pilea profile *.fasta -d db -o .
```

PTR estimates and other metadata for MAGs that pass basic filters will be saved to `output.tsv`:
```text
...  genome  taxonomy                   coverage  dispersion  fraction  containment  PTR     log2(PTR)
...  B       ...Bacillus subtilis       5.9939    0.9057      1.0000    0.9936       1.5451  0.6277
...  K       ...Klebsiella pneumoniae   12.0417   0.9715      0.9956    0.9867       1.2804  0.3566
...  M       ...Morganella morganii     9.9278    0.9222      0.9716    0.9766       1.3799  0.4646
...  P       ...Pseudomonas_E putida_U  5.6176    0.8939      0.9955    0.9859       1.7704  0.8240
```
Coverage (`--min-cove`) represents the median per-window coverage estimated from sketched $k$-mers (lower than per-base coverage, depending on read length and $k$). Dispersion (`--max-disp`) and fraction (`--min-frac`) indicate the median per-window dispersion and the fraction of covered windows, respectively. Containment (`--min-cont`) is the proportion of sketched $k$-mers used for PTR estimation. See `pilea profile -h` for more details.

## Citation
Chen, X., Xu, X., Lin, Y., Shi, X., Wang, D., & Zhang, T. (2026). Pilea: profiling bacterial growth dynamics from metagenomes with sketching. *Microbiome*. https://doi.org/10.1186/s40168-026-02374-0
