# mutate_seq - High-Performance Sequence Mutation Tool

**Version 0.5* - Tool for introducing controlled mutations into genomic sequences.

## Overview

`mutate_seq` Introduces controlled mutations, ancient DNA damage, and fragmentation into FASTA/FASTQ sequences.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Usage](#usage)
- [Mutation Modes](#mutation-modes)
- [Ancient DNA Damage Simulation](#ancient-dna-damage-simulation)
- [DNA Fragmentation](#dna-fragmentation)
- [Output Files](#output-files)
- [Processing Modes](#processing-modes)
- [Troubleshooting](#troubleshooting)
- [Technical Details](#technical-details)
- [Citation](#citation)
- [Contact](#contact)
- [Contributing](#contributing)
- [License](#license)

---


## Installation
### Clone Repository

```bash
git clone https://github.com/Madshartmann1/mutation_tool.git
cd mutation_tool
```
### Build
```bash
make
```

Check dependencies:
```bash
make check-deps
```

### Requirements
- C++11 compiler (g++ 4.8+)
- zlib library (`zlib1g-dev` on Ubuntu/Debian)
- pthread support

## Usage

### Basic Syntax
```bash
./mutate_seq --input <file> --output <prefix> [mutation_mode] [options]
```
<!---
### Quick Examples

#### 1. Flat Mutation Rate (0.1%)
```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --seed 27
```

#### 2. Fixed Number of Mutations
```bash
./mutate_seq --input genome.fa --output mutated --num-mutations 1000 --seed 27
```

#### 3. Ts/Tv Ratio Control
```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --ts-tv-ratio 2.0 --seed 27
```

#### 4. Multi-threaded FASTQ Processing
```bash
./mutate_seq --input reads.fastq.gz --output mutated --mutation-rate 0.0001 --threads 8 --seed 27
```
--->


## Mutation Modes

### Mode 1: Flat Rate (`--mutation-rate`)
Applies uniform mutation rate across all substitution types.

```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --seed 27
```

### Mode 2: Flat Number (`--num-mutations`)
Introduces exactly N mutations across the entire dataset.

```bash
./mutate_seq --input genome.fa --output mutated --num-mutations 500 --seed 27
```

**Note**: Requires in-memory mode (not available for streaming).

### Mode 3: Separate Ts/Tv Rates (`--ts-rate` + `--tv-rate`)
Specifies different rates for transitions and transversions.

```bash
./mutate_seq --input genome.fa --output mutated --ts-rate 0.002 --tv-rate 0.001 --seed 27
```

### Mode 4: Ts/Tv Ratio (`--mutation-rate` + `--ts-tv-ratio`)
Sets overall mutation rate with specified Ts/Tv ratio.

```bash
./mutate_seq --input genome.fa --output mutated --mutation-rate 0.001 --ts-tv-ratio 2.0 --seed 27
```

**Calculation**: `ts_rate = rate × ratio / (ratio + 1)`, `tv_rate = rate / (ratio + 1)`

### Mode 5: Custom Mutation Matrix (`--mutation-matrix` + `--mutation-rate`)

**Full control over every mutation rate directly.** A certain bases mutation rate is the sum of its parts. 

**Example**: If A→G:0.005 + A→C:0.02 + A→T:0.01, then A bases mutate at 0.035% total. If C→A:0.008 + C→G:0.012 + C→A:0.022, then C bases mutate at 0.062% total.


**Matrix file format** (`mutation_matrix.txt`):
```
# Format: from_base  to_base  relative_rate
A       C       0.005
A       G       0.02
A       T       0.01
C       A       0.008
C       G       0.012
C       T       0.025
G       A       0.022
G       C       0.009
G       T       0.011
T       A       0.01
T       C       0.023
T       G       0.007
```

```bash
./mutate_seq --input genome.fa --output mutated \
    --mutation-matrix mutation_matrix.txt --mutation-rate 1 --seed 27
```


**PS: The Mutation rate is currently igonred for the matrix**

### Mode 6: Custom Mutation Spectrum (`--mutation-spectrum` + `--mutation-rate` or `--num-mutations`)

**All bases mutate at the same rate; spectrum controls the distribution of mutation types.** Values are normalized so all 12 substitution types sum to 1.0.

**Example**: With `--mutation-rate 0.006`, all bases mutate at 0.6%. If a bases then mutates, the spectrum determines that 50% of the mutations starting from A is A→G etc. See spectrum example below.

**Spectrum file format** (`mutation_spectrum.txt`):
```
# Format: from_base  to_base  proportion
# Values will be normalized to sum to 1
# Example: Transition-biased spectrum
A       C       1.0
A       G       2.0
A       T       1.0
C       A       1.0
C       G       1.0
C       T       2.0
G       A       2.0
G       C       1.0
G       T       1.0
T       A       1.0
T       C       2.0
T       G       1.0
```

```bash
./mutate_seq --input genome.fa --output mutated \
    --mutation-spectrum mutation_spectrum.txt --mutation-rate 0.001 --seed 27
```

## Ancient DNA Damage Simulation

`mutate_seq` can simulate position-dependent ancient DNA (aDNA) damage patterns, example:
- **C→T deamination** enriched at 5' fragment ends
- **G→A deamination** enriched at 3' fragment ends  
- **Exponential decay** of damage signal away from fragment ends

Ancient damage can be **combined with any mutation mode** (or used standalone).


**Just damage (`--ancient-damage`)**
```bash
./mutate_seq --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --threads 8
```

**Damage + mutation rate + Ts/Tv ratio (`--ancient-damage`+`--mutation-rate`+`--ts-tv-ratio`)**
```bash
./mutate_seq --input reads.fastq.gz --output evolved_damaged \
    --mutation-rate 0.001 --ts-tv-ratio 2.0 \
    --ancient-damage damage_profiles/single_stranded_damage.txt \
    --threads 8
```


**Damage + Background mutations (`--ancient-damage`+`--background-rate`)**

Background rate only effects poistions left untouched by the damage profile. 
```bash
./mutate_seq --input reads.fastq.gz --output complex \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --background-rate 0.0001 \
    --threads 8
```

**See also**: [damage_profiles/](damage_profiles/) directory for detailed damage profile documentation.

---

## DNA Fragmentation

`mutate_seq` can fragment long sequences into realistic DNA fragment length distributions.

**Processing order**: Fragmentation happens **first**, then mutations, then damage (if enabled).

### Fragmentation Modes

#### 1. Static Length (Fixed Size)

All fragments exactly N bp:
```bash
./mutate_seq --input genome.fa --output fragments \
    --fragment-length 100 \
    --mutation-rate 0.001
```

#### 2. Empirical / Costum Distribution

Sample from observed fragment lengths (e.g., from ancient DNA datasets):
```bash
./mutate_seq --input genome.fa --output ancient_frags \
    --fragment-distribution empirical \
    --fragment-distribution-file fragment_distributions/ancient_dist_chagyrskaya8.txt \
    --mutation-rate 0.001
```

**Included distribution**: `ancient_dist_chagyrskaya8.txt` contains 1,000,001 fragment lengths from Chagyrskaya 8 Neanderthal sample (mean ~60bp, range 35-130bp).

**Custom distributions**: Provide a text file with one fragment length per line:
```
50
75
60
48
```

#### 3. Exponental Distribution

```bash
./mutate_seq --input genome.fa --output frags \
    --fragment-distribution exponential \
    --mean-length 80 \
    --mutation-rate 0.001
```

#### 4. Normal Distribution

```bash
./mutate_seq --input genome.fa --output frags \
    --fragment-distribution normal \
    --mean-length 150 \
    --sd-length 30 \
    --mutation-rate 0.001
```

#### 5. Lognormal Distribution

```bash
./mutate_seq --input genome.fa --output frags \
    --fragment-distribution lognormal \
    --mean-length 100 \
    --sd-length 50 \
    --mutation-rate 0.001
```

### Minimum Fragment Length Filtering

By default, fragments shorter than **20bp** are discarded. Change this threshold:

```bash
./mutate_seq --input genome.fa --output frags \
    --fragment-distribution empirical \
    --fragment-distribution-file ancient_dist_chagyrskaya8.txt \
    --min-fragment-length 30 \
    --mutation-rate 0.001
```

### Fragmentation Behavior

**Input**: One sequence (e.g., chromosome or long read)  
**Output**: Multiple fragments with sequential naming:
```
>chr1_frag1
ACGTACGT...
>chr1_frag2  
TGCAGCTA...
>chr1_frag3
CGATCGAT...
```

**Fragment size selection**:
1. Sample a fragment length from the distribution
2. If length > remaining sequence → check if remainder ≥ min_fragment_length
   - If yes: keep short fragment
   - If no: discard remainder
3. If entire read < first sampled length → discard entire read 
4. Extract new read and do step 1


### Complete Ancient DNA Simulation Example

Combine fragmentation + mutations + damage for realistic ancient DNA:

```bash
./mutate_seq --input modern_reference.fa --output ancient_sample \
    --fragment-distribution empirical \
    --fragment-distribution-file fragment_distributions/ancient_dist_chagyrskaya8.txt \
    --min-fragment-length 30 \
    --mutation-rate 0.0015 --ts-tv-ratio 2.0 \
    --ancient-damage damage_profiles/single_stranded_damage.txt \
    --background-rate 0.0001 \
    --threads 8 --seed 27
```

**This simulates**:
1. **Fragmentation**: Realistic ancient fragment lengths (30-130bp)
2. **Evolutionary divergence**: 0.15% substitutions with Ts/Tv bias
3. **Ancient damage**: C→T/G→A at fragment ends
4. **Background mutations**: Additional 0.01% random substitutions


<!---
---

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

### Mutation Mode
- `--mutation-rate FLOAT`: Uniform mutation rate (0.0-1.0)
- `--num-mutations INT`: Fixed number of mutations
- `--ts-rate FLOAT --tv-rate FLOAT`: Separate transition/transversion rates
- `--mutation-rate FLOAT --ts-tv-ratio FLOAT`: Overall rate + Ts/Tv ratio
- `--mutation-matrix FILE --mutation-rate FLOAT`: Custom per-type rates
- `--mutation-spectrum FILE --mutation-rate FLOAT`: Custom proportional spectrum
Ancient DNA Damage
- `--ancient-damage FILE`: Position-dependent damage profile
- `--background-rate FLOAT`: Additional uniform mutation rate

### DNA Fragmentation  
- `--fragment-length INT`: Static fragment length (all N bp)
- `--fragment-distribution TYPE`: Distribution type (empirical, exponential, normal, lognormal)
- `--fragment-distribution-file FILE`: Empirical distribution (one length per line)
- `--mean-length INT`: Mean for parametric distributions
- `--sd-length INT`: Standard deviation (normal/lognormal)
- `--min-fragment-length INT`: Minimum fragment to keep (default: 20)

### 
### Optional
- `--seed INT`: Random seed for reproducibility (default: 27)
- `--threads INT`: Number of threads (default: 4, min: 1)
- `--chunk-size INT`: Sequences per chunk in streaming mode (default: 10000)
- `--format [fasta|fastq]`: Force format (auto-detected if omitted)

---
--->


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



### Processing Modes

#### In-Memory Mode
- **Triggered by**: FASTA input, `--num-mutations`, or `--threads 1`
- **Characteristics**: Loads ALL sequences into RAM, enables detailed statistics
- **Best for**: Reference genomes (<10GB), fixed mutation count mode
- **Memory usage**: ~2-3× file size

#### Streaming Mode
- **Triggered by**: FASTQ input with 2+ threads
- **Characteristics**: Processes chunks, merges outputs, minimal memory
- **Best for**: Large FASTQ files (50GB+), read datasets
- **Memory usage**: ~50-100MB regardless of file size
- **Threads**: 1 producer + (N-1) workers (min 2 threads required)


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

This project is currently in reasearch stages.
