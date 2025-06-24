# Progressive Unique Kmer Counter (pukc)

A Rust-based tool for visualizing the number of solid canonical kmers while streaming reads.

---

## Features

- Supports **FASTA** and **FASTQ**
- Real-time WebSocket output for monitoring (Disabled in code)
- **Growth**: The number of new solid k-mers between read intervals.
- **Acceleration**: The second derivative of k-mer discovery, indicating whether the rate of discovery is increasing, decreasing, or stabilizing
- Early termination based on low and stable accelerations 
- For fasta and fastq files: estimates the number of reads and extrapolates an estimate of the total number of kmers. 
- /!\ This does not work for .gz files (as we cannot easily estimate the size of the file and so the number of reads)


---

## Installation 

```bash

git clone https://github.com/pierrepeterlongo/progressive_unique_kmers_counter.git
cd progressive_unique_kmers_counter
cargo install --path .  
```

## Usage

```bash

Usage: pukc [OPTIONS] --k <K> --input <INPUT>

Options:
  -k, --k <K>
          Length of k-mers
  -i, --input <INPUT>
          Input FASTA file
      --nb-reads <NB_READS>
          Give the number of reads [default: 0]
      --nb-stable-curves-found-max <NB_STABLE_CURVES_FOUND_MAX>
          Maximal number of stable curves to find before stopping [default: 50]
  -v, --verbose
          Show progress and details
  -d, --debug
          Show debug information
  -h, --help
          Print help
  -V, --version
          Print version

```

## Example
```
pukc --k 25 --input input.fq 
```
- In this case the solid distinct canonical 25-mers are counted 
- Then the number of kmers is extrapolated 

## Full example.
1/ Download the fasta.gz or fastq.gz entry from [SRR33792312](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR33792312&display=metadata). This file contains 2.2M reads (2113824 precisely).
2/ Run :
```bash
pukc --k 25 --input /tmp/SRR33792312.fasta.gz  --nb-reads 2113824

...

Estimated total unique k-mers: 7431136
```

- Note1: this execution takes approximately 2s. 
- Note2: the exact number of distinct solide 25-mers in this file is in fact 5239190. Pukc tends to surestimate the results.

## Visualization
- Redirect logs to a file: 
```bash
pukc --k 25 --input /tmp/SRR33792312.fasta.gz  --nb-reads 2113824 > logs.txt
```
- Plot figure: 
```bash
python plot_progression.py logs.txt
```
- Generated figure: 

TODO