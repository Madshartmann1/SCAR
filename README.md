# mutate_seq - High-Performance Sequence Mutation Tool

**Version 2.0** - Production-ready tool for introducing controlled mutations into genomic sequences.

## Overview

`mutate_seq` introduces random mutations into DNA sequences (FASTA or FASTQ format) at specified rates or with custom mutation spectra. This tool is designed for researchers and bioinformaticians who need to simulate sequence variation with precise control over mutation patterns, whether for testing bioinformatics pipelines, training machine learning models, or simulating evolutionary processes.

### Key Features

- **6 Mutation Modes**: From simple flat rates to custom mutation spectra
- **Multi-threading**: Scalable to 92+ cores with streaming mode
- **Gzip Support**: Transparent compression/decompression for input/output
- **Auto-detection**: Automatically detects file format (FASTA/FASTQ) and compression
- **SNP Receipt**: Generates seqtk-compatible SNP files for tracking mutations
- **Memory Efficient**: Streaming mode for large FASTQ files (50GB+)
- **Quality Preservation**: FASTQ quality scores remain unchanged

## Installation

### Requirements
- C++11 compiler (g++ 4.8+)
- zlib library (`zlib1g-dev` on Ubuntu/Debian)
- pthread support

### Build
```bash
make
```

Check dependencies:
```bash
make check-deps
```

## Usage

### Basic Syntax
```bash
./mutate_seq --input <file> --output <prefix> [mutation_mode] [options]
```

### Quick Examples

#### 1. Flat Mutation Rate (0.1%)
```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --seed 42
```

#### 2. Fixed Number of Mutations
```bash
./mutate_seq --input genome.fa --output mutated --num-mutations 1000 --seed 42
```

#### 3. Ts/Tv Ratio Control
```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --ts-tv-ratio 2.0 --seed 42
```

#### 4. Multi-threaded FASTQ Processing
```bash
./mutate_seq --input reads.fastq.gz --output mutated --mutation-rate 0.0001 --threads 8 --seed 42
```

## Mutation Modes

### Mode 1: Flat Rate (`--mutation-rate`)
Applies uniform mutation rate across all substitution types.

```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --seed 42
```

**Use case**: Simple uniform mutations for baseline simulations.

### Mode 2: Flat Number (`--num-mutations`)
Introduces exactly N mutations across the entire dataset.

```bash
./mutate_seq --input genome.fa --output mutated --num-mutations 500 --seed 42
```

**Use case**: Fixed mutation count for controlled experiments.  
**Note**: Requires in-memory mode (not available for streaming).

### Mode 3: Separate Ts/Tv Rates (`--ts-rate` + `--tv-rate`)
Specifies different rates for transitions and transversions.

```bash
./mutate_seq --input genome.fa --output mutated --ts-rate 0.002 --tv-rate 0.001 --seed 42
```

**Use case**: Modeling biological Ts/Tv bias (typically 2:1 in vertebrates).

### Mode 4: Ts/Tv Ratio (`--mutation-rate` + `--ts-tv-ratio`)
Sets overall mutation rate with specified Ts/Tv ratio.

```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --ts-tv-ratio 2.0 --seed 42
```

**Use case**: Biological realism with automatic rate calculation.  
**Calculation**: `ts_rate = rate × ratio / (ratio + 1)`, `tv_rate = rate / (ratio + 1)`

### Mode 5: Custom Mutation Matrix (`--mutation-matrix` + `--mutation-rate`)
Per-substitution-type mutation rates for maximum control.

**Matrix file format** (`mutation_matrix.txt`):
```
# Format: from_base  to_base  relative_rate
A       C       0.5
A       G       2.0
A       T       1.0
C       A       0.8
C       G       1.2
C       T       2.5
G       A       2.2
G       C       0.9
G       T       1.1
T       A       1.0
T       C       2.3
T       G       0.7
```

```bash
./mutate_seq --input genome.fa --output mutated \
    --mutation-matrix mutation_matrix.txt --mutation-rate 0.001 --seed 42
```

**Use case**: Ancient DNA damage patterns, context-dependent mutations, or any scenario requiring fine-grained control over specific substitution types.

### Mode 6: Custom Mutation Spectrum (`--mutation-spectrum` + `--mutation-rate`)
Proportional distribution of mutations (values normalized to sum=1).

**Spectrum file format** (`mutation_spectrum.txt`):
```
# Format: from_base  to_base  proportion
# Values will be normalized to sum to 1
# Example: Transition-biased spectrum
A       G       2.0
G       A       2.0
C       T       2.0
T       C       2.0
A       C       1.0
A       T       1.0
C       A       1.0
C       G       1.0
G       C       1.0
G       T       1.0
T       A       1.0
T       G       1.0
```

```bash
./mutate_seq --input genome.fa --output mutated \
    --mutation-spectrum mutation_spectrum.txt --mutation-rate 0.001 --seed 42
```

**Use case**: Modeling mutational signatures (e.g., APOBEC, UV damage, smoking signatures) or any proportional mutation distribution.

## Pipeline Integration Example

Example workflow for simulating sequences with varying divergence levels for read simulation:

```bash
# Step 1: Generate reference variants with different mutation levels
REFERENCE="reference_genome.fa"
SEED=12345

# Baseline: no mutations
cp ${REFERENCE} sample_baseline.fa

# High divergence: 0.1% mutation rate
./mutate_seq --input ${REFERENCE} --output sample_high \
    --mutation-rate 0.001 --threads 8 --seed ${SEED}

# Low divergence: 0.01% mutation rate
./mutate_seq --input ${REFERENCE} --output sample_low \
    --mutation-rate 0.0001 --threads 8 --seed $((SEED + 1))

# Step 2: Simulate sequencing reads (example using ART)
for variant in baseline high low; do
    art_illumina -ss HS25 -i sample_${variant}.fa -l 150 -f 30 \
        -o sample_${variant}_reads -rs ${SEED}
done

# Step 3: SNP tracking files are available for analysis
# Files generated: sample_high.snp and sample_low.snp
```

## Command-Line Options

### Required
- `--input FILE`: Input FASTA/FASTQ file (plain or gzipped)
- `--output PREFIX`: Output file prefix (format auto-matches input)

### Mutation Mode (choose ONE)
- `--mutation-rate FLOAT`: Uniform mutation rate (0.0-1.0)
- `--num-mutations INT`: Fixed number of mutations
- `--ts-rate FLOAT --tv-rate FLOAT`: Separate transition/transversion rates
- `--mutation-rate FLOAT --ts-tv-ratio FLOAT`: Overall rate + Ts/Tv ratio
- `--mutation-matrix FILE --mutation-rate FLOAT`: Custom per-type rates
- `--mutation-spectrum FILE --mutation-rate FLOAT`: Custom proportional spectrum

### Optional
- `--seed INT`: Random seed for reproducibility (default: 42)
- `--threads INT`: Number of threads (default: 4, min: 1)
- `--chunk-size INT`: Sequences per chunk in streaming mode (default: 10000)
- `--format [fasta|fastq]`: Force format (auto-detected if omitted)
- `--force-streaming`: Force streaming mode even for FASTA

## Output Files

### Mutated Sequences
Format matches input:
- Input: `genome.fa` → Output: `prefix.fa`
- Input: `reads.fastq.gz` → Output: `prefix.fastq.gz`

### SNP File (`prefix.snp`)
Seqtk-compatible format (tab-delimited):
```
seq_1       45      A       G
seq_1       127     C       T
seq_2       89      G       A
```

Columns: `sequence_name`, `position` (1-based), `original_base`, `mutated_base`

## Performance

### Benchmarks (10k reads × 150bp, mutation rate 0.001)

| Threads | Time    | Throughput    | Mode      |
|---------|---------|---------------|-----------|
| 1       | 10.0s   | 1000 seq/sec  | In-memory |
| 2       | 10.3s   | 1000 seq/sec  | Streaming |
| 4       | 10.3s   | 1000 seq/sec  | Streaming |
| 8       | 10.3s   | 1000 seq/sec  | Streaming |

**Note**: Small files are I/O bound (gzip). Threading benefits appear with larger datasets (100k+ reads).

### Processing Modes

#### In-Memory Mode
- **Triggered by**: FASTA input, `--num-mutations`, or `--threads 1`
- **Characteristics**: Loads all sequences into RAM, enables detailed statistics
- **Best for**: Reference genomes (<10GB), fixed mutation count mode
- **Memory usage**: ~2-3× file size

#### Streaming Mode
- **Triggered by**: FASTQ input with 2+ threads
- **Characteristics**: Processes chunks, merges outputs, minimal memory
- **Best for**: Large FASTQ files (50GB+), read datasets
- **Memory usage**: ~50-100MB regardless of file size
- **Threads**: 1 producer + (N-1) workers (min 2 threads required)

## Validation

### Test Suite
```bash
# Generate test data
python3 generate_test_data.py

# Test all modes
make -f Makefile_v2 test
```

### Expected Output
- **Mutation rate accuracy**: Within ~2% of specified rate
- **Ts/Tv ratio**: Within 10% of target (small genomes have higher variance)
- **Output format**: Matches input format and compression
- **SNP file**: All mutations logged in seqtk format

## Troubleshooting

### Compilation Errors

**Error**: `fatal error: zlib.h: No such file or directory`  
**Fix**: Install zlib development package
```bash
# Ubuntu/Debian
sudo apt-get install zlib1g-dev

# CentOS/RHEL
sudo yum install zlib-devel

# macOS
brew install zlib
```

### Runtime Issues

**Issue**: Multiple mutation modes specified  
**Fix**: Use only ONE mutation mode. Example of incorrect usage:
```bash
# INCORRECT: --num-mutations conflicts with --mutation-rate
./mutate_seq --input file.fa --output out --num-mutations 100 --mutation-rate 0.01
```

**Issue**: No mutations in output (0 in SNP file)  
**Cause**: Mutation rate too low for small genome, or stochastic variation with random seed  
**Fix**: Increase mutation rate or use `--num-mutations` for guaranteed count

**Issue**: Unexpected mutation counts  
**Note**: For small genomes or low mutation rates, stochastic variation is expected. Use `--num-mutations` for deterministic counts.

### Performance Issues

**Issue**: Slow processing with streaming mode  
**Check**: 
1. Is file gzipped? (Decompression is CPU-intensive)
2. Are you using enough threads? (Try 4-8 for large files)
3. Is disk I/O the bottleneck? (Monitor with `iostat -x 1`)

**Tip**: For repeated processing, decompress once and reuse:
```bash
gunzip -c reads.fastq.gz > reads.fastq
# Process multiple times without decompression overhead
./mutate_seq --input reads.fastq --output test1 --mutation-rate 0.001 --threads 8
./mutate_seq --input reads.fastq --output test2 --mutation-rate 0.0001 --threads 8
```

## Technical Details

### Threading Architecture
- **Producer thread**: Reads input file, distributes chunks to work queue
- **Worker threads**: Process mutations, write to thread-specific temporary files
- **Main thread**: Merges worker outputs into final file
- **Progress thread**: Reports status every 10 seconds

### Random Number Generation
- C++ `<random>` library with Mersenne Twister (mt19937_64)
- Seeded for reproducibility (`--seed` option)
- Each worker thread has independent RNG (seeded from main seed)

### Quality Score Handling
- FASTQ quality scores preserved byte-for-byte
- Mutations only affect sequence line (line 2 of each record)
- Quality line (line 4) remains unchanged

## Changelog

### Version 2.0 (Current)
- Added streaming mode with multi-threading
- Added gzip support (input and output)
- Added custom mutation spectrum mode
- Added progress monitoring
- Fixed single-thread mode for FASTQ
- Fixed mutation mode selection logic

### Version 1.0
- Initial implementation with in-memory processing
- 5 mutation modes (flat rate, flat number, separate Ts/Tv, Ts/Tv ratio, custom matrix)
- FASTA support only

## Citation

If you use this tool in published research, please cite this repository:

```
mutate_seq: High-Performance Sequence Mutation Tool
https://github.com/Madshartmann1/mutation_tool
```

## Contact

For questions or support, please:
- Open an issue on GitHub
- Email: madhar@dtu.dk

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This project is open source. You are free to use, modify, and distribute this software for research and academic purposes.
