# unique_kmers_evolution

A Rust-based tool for visualizing the number of solid canonical kmers while streaming reads.

---

## Features

- Supports **FASTA** and **FASTQ**
- Supports **gzip-compressed** files (`.gz`). **WARNING** In this case the extrapolation of number of reads and kmers is false. TODO: add an option to inform the number of reads. 
- Real-time WebSocket output for monitoring
- **Growth**: The number of new solid k-mers between read intervals.
- **Acceleration**: The second derivative of k-mer discovery, indicating whether the rate of discovery is increasing, decreasing, or stabilizing
- Early termination based on configurable acceleration threshold
- For fasta and fastq files: estimates the number of reads and extrapolates an estimate of the total number of kmers. 
- /!\ This does not work for .gz files (as we cannot easily estimate the size of the file and so the number of reads)


---

## Installation 

```bash

git clone https://github.com/yourusername/unique_kmers_evolution.git
cd unique_kmers_evolution
cargo build --release
cargo install --path .  
```

## Usage

```bash

Usage: unique_kmers_evolution [OPTIONS] --k <K> --input <INPUT>

Options:
  -k, --k <K>
          Length of k-mers
  -i, --input <INPUT>
          Input FASTA file
      --stop-acceleration <STOP_ACCELERATION>
          Stop if acceleration is below this threshold [default: 10]
  -h, --help
          Print help
  -V, --version
          Print version
```

## Example
```
unique_kmers_evolution --k 25 --input input.fq --stop-acceleration 20
```
- In this case the solid 25-mers are counted as long as the acceleration over 50000 reads is higher than 20.
- Then the number of kmers is extrapolated 


**Visualize** the evolution of the results opening file `plot.html` in a browser (reload the page once the program runs)
