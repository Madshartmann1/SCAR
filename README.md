# SCAR — Sequences + Controlled mutations + Ancient damage + fRagmentation

**Version 1.0* 

`SCAR` Introduces controlled mutations, ancient DNA damage, and fragmentation into FASTA/FASTQ sequences.

## Table of Contents

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
git clone https://github.com/Madshartmann1/SCAR.git
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
./scar --input <file> --output <prefix> [mutation_mode] [options]
```

For mapper-ready paired-end mutation, keep R1 and R2 separate:
```bash
./scar --input-r1 R1.fastq.gz --input-r2 R2.fastq.gz \
    --output-r1 mutated_R1 --output-r2 mutated_R2 \
    --mutation-rate 0.001 --threads 8 --seed 27
```

### Expected Input

SCAR expects FASTA or FASTQ input. For read-level FASTQ input, reads are assumed to have already gone through adapter removal and quality trimming before running SCAR.

For paired-end mutation, provide two synchronized FASTQ files with records in the
same order. SCAR reads one R1 record and one R2 record at a time, checks that the
normalized read names match, and aborts on pair-name or record-count mismatch.
Headers such as `@ERR13719744.1` in both mates are valid. `/1` and `/2` suffixes
are ignored for pair matching.

Paired-end mode uses deterministic per-read random seeds from the global seed,
pair index, normalized pair name, and mate label. With the same input and `--seed`,
paired-end output is stable across thread counts and chunk sizes.
<!---
### Quick Examples

#### 1. Flat Mutation Rate (0.1%)
```bash
./scar --input genome.fa --output mutated --mutation-rate 0.001 --seed 27
```

#### 2. Fixed Number of Mutations
```bash
./scar --input genome.fa --output mutated --num-mutations 1000 --seed 27
```

#### 3. Ts/Tv Ratio Control
```bash
./scar --input genome.fa --output mutated --mutation-rate 0.001 --ts-tv-ratio 2.0 --seed 27
```

#### 4. Multi-threaded FASTQ Processing
```bash
./scar --input reads.fastq.gz --output mutated --mutation-rate 0.0001 --threads 8 --seed 27
```
--->


## Mutation Modes

### Mode 1: Flat Rate (`--mutation-rate`)
Applies uniform mutation rate across all substitution types.

```bash
./scar --input genome.fa --output mutated --mutation-rate 0.001 --seed 27
```

Paired-end FASTQ example:
```bash
./scar --input-r1 reads_R1.fastq.gz --input-r2 reads_R2.fastq.gz \
    --output-r1 mutated_R1 --output-r2 mutated_R2 \
    --mutation-rate 0.001 --threads 8 --seed 27
```

### Mode 2: Flat Number (`--num-mutations`)
Introduces exactly N mutations across the entire dataset.

```bash
./scar --input genome.fa --output mutated --num-mutations 500 --seed 27
```

**Note**: Requires in-memory mode (not available for streaming).

**Paired-end note**: `--num-mutations` is not supported in paired-end mode. Use
a rate-based mutation mode for synchronized R1/R2 output.

### Mode 3: Separate Ts/Tv Rates (`--ts-rate` + `--tv-rate`)
Specifies different rates for transitions and transversions.

```bash
./scar --input genome.fa --output mutated --ts-rate 0.002 --tv-rate 0.001 --seed 27
```

### Mode 4: Ts/Tv Ratio (`--mutation-rate` + `--ts-tv-ratio`)
Sets overall mutation rate with specified Ts/Tv ratio.

```bash
./scar --input genome.fa --output mutated \
 --mutation-rate 0.001 --ts-tv-ratio 2.0 --seed 27
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
./scar --input genome.fa --output mutated \
    --mutation-matrix mutation_matrix.txt --mutation-rate 1 --seed 27
```


**PS: The mutation rate is currently ignored for the matrix**

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
./scar --input genome.fa --output mutated \
    --mutation-spectrum mutation_spectrum.txt --mutation-rate 0.001 --seed 27
```

## Ancient DNA Damage Simulation

`scar` can simulate position-dependent ancient DNA (aDNA) damage patterns, example:
- **C→T deamination** enriched at 5' fragment ends
- **G→A deamination** enriched at 3' fragment ends  
- **Exponential decay** of damage signal away from fragment ends

Ancient damage can be **combined with any mutation mode** (or used standalone).


**Just damage (`--ancient-damage`)**
```bash
./scar --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --threads 8
```

`--ancient-damage` also accepts a mapDamage output directory containing
`misincorporation.txt`; SCAR auto-detects the format and converts it in memory.

```bash
./scar --input reads.fastq.gz --output damaged \
    --ancient-damage mapdamage_results/ \
    --threads 8
```

**Damage + mutation rate + Ts/Tv ratio (`--ancient-damage`+`--mutation-rate`+`--ts-tv-ratio`)**
```bash
./scar --input reads.fastq.gz --output evolved_damaged \
    --mutation-rate 0.001 --ts-tv-ratio 2.0 \
    --ancient-damage damage_profiles/single_stranded_damage.txt \
    --threads 8
```


**Damage + Background mutations (`--ancient-damage`+`--background-rate`)**

Background rate only affects positions left untouched by the damage profile. 
```bash
./scar --input reads.fastq.gz --output complex \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --background-rate 0.0001 \
    --threads 8
```

**See also**: [damage_profiles/](damage_profiles/) directory for detailed damage profile documentation.

---

## DNA Fragmentation

`scar` can fragment long sequences into realistic DNA fragment length distributions.

**Processing order**: Fragmentation happens **first**, then mutations, then damage (if enabled).

**Paired-end limitation**: Fragmentation is single-end logic. Paired-end mode
(`--input-r1`/`--input-r2`) rejects fragmentation because independent fragmenting
would break mapper-ready R1/R2 synchronization. For paired-end datasets that need
fragmentation, use the merge-for-fragmentation helper workflow; that produces a
single-end merged FASTQ and is not mapper-ready paired-end output. For mapper-ready
R1/R2 output, use paired-end mutation/damage without fragmentation.

**FASTA Safeguard**: Fragmentation and/or ancient damage with FASTA input is blocked by default (FASTA lacks base quality scores needed for realistic reads). Use a read simulator (ART, wgsim, etc.) to generate FASTQ first, or override with `--allow-fasta-fragmentation` for testing purposes.

### Fragmentation Modes

#### 1. Static Length (Fixed Size)

All fragments exactly N bp:
```bash
./scar --input reads.fastq.gz --output fragments \
    --fragment-length 100 \
    --mutation-rate 0.001
```

#### 2. Empirical / Custom Distribution

Sample fragment lengths from observed data. Empirical files can be simple length lists,
weighted length tables, or mapDamage length distribution output.

Example using a generic empirical distribution file:
```bash
./scar --input input_reads.fastq.gz --output empirical_fragments \
    --fragment-distribution empirical \
    --fragment-distribution-file fragment_lengths.txt \
    --mutation-rate 0.001
```

This samples one row at random for each fragment. If the file has one length per
line, each listed length has equal weight. Repeated values therefore naturally
increase the probability of that length.

**Included distribution**: `ancient_dist_chagyrskaya8.txt` contains 1,000,001 fragment lengths from Chagyrskaya 8 Neanderthal sample (mean ~60bp, range 35-130bp).

**Format 1: one fragment length per line**
```
50
75
60
48
```

Use this format when you have raw observed fragment lengths. Every row is treated
as one observation.

**Format 2: weighted length table**
```
50 120
75 450
100 80
125 20
```

Equivalent probability-style weights also work:
```
50 0.179
75 0.672
100 0.119
125 0.030
```

The first column is fragment length and the second column is the relative
sampling weight. These values are normalized internally, so they can be absolute
occurrence counts from an empirical dataset, percentages that sum to 100, or
probabilities that sum to 1. In this example, 75 bp fragments are sampled more
often than 125 bp fragments.

**Format 3: mapDamage length distribution**
```
# table produced by mapDamage
Std	Length	Occurrences
+	35	42
+	36	58
+	37	91
-	35	38
-	36	61
```

For mapDamage-style files, `scar` reads the length from the second column and the
occurrence count from the third column. Comment lines and text headers are skipped
automatically.

Example using a mapDamage `lgdistribution.txt` file:
```bash
./scar --input simulated_reads.fastq.gz --output mapdamage_fragments \
    --fragment-distribution empirical \
    --fragment-distribution-file mapdamage_results/lgdistribution.txt \
    --mutation-rate 0.001 \
    --seed 27
```

This uses the observed mapDamage length distribution directly, so lengths with
higher occurrence counts are sampled more frequently.

#### 3. Exponential Distribution

```bash
./scar --input reads.fastq.gz --output frags \
    --fragment-distribution exponential \
    --mean-length 80 \
    --mutation-rate 0.001
```

#### 4. Normal Distribution

```bash
./scar --input reads.fastq.gz --output frags \
    --fragment-distribution normal \
    --mean-length 150 \
    --sd-length 30 \
    --mutation-rate 0.001
```

#### 5. Lognormal Distribution

```bash
./scar --input reads.fastq.gz --output frags \
    --fragment-distribution lognormal \
    --mean-length 100 \
    --sd-length 50 \
    --mutation-rate 0.001
```

### Fragment Length Filtering

**Minimum fragment length**: By default, fragments shorter than **20bp** are discarded. Change this threshold:

```bash
./scar --input reads.fastq.gz --output frags \
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
2. If the first sampled length is longer than the whole read, discard that read
3. If a later sampled length is longer than the remaining sequence, keep the remainder only if it is at least `--min-fragment-length`
   - If yes: keep short fragment
   - If no: discard remainder
4. Extract each fragment with sequential `_fragN` names and repeat until the read is consumed


### Complete Ancient DNA Simulation Example

Combine fragmentation + mutations + damage for realistic ancient DNA:

```bash
# First generate FASTQ reads from reference using ART or similar read simulator
# Then apply fragmentation, mutations, and damage:
./scar --input simulated_reads.fastq.gz --output ancient_sample \
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




## Output Files

### Mutated Sequences
Format matches input:
- Input: `genome.fa` → Output: `prefix.fa`
- Input: `reads.fastq.gz` → Output: `prefix.fastq.gz`
- Input: `--input-r1 reads_R1.fastq.gz --input-r2 reads_R2.fastq.gz` →
  Outputs: `output_R1.fastq.gz` and `output_R2.fastq.gz`


### SNP File (`prefix.snp[.gz]`)
Seqtk-compatible format (tab-delimited):
```
seq_1       45      A       G
seq_1       127     C       T
seq_2       89      G       A
```

Columns: `sequence_name`, `position` (1-based), `original_base`, `mutated_base`

In paired-end mode, SCAR writes one SNP receipt per mate:
- `--output-r1 sample_R1` → `sample_R1.fastq[.gz]` and `sample_R1.snp[.gz]`
- `--output-r2 sample_R2` → `sample_R2.fastq[.gz]` and `sample_R2.snp[.gz]`

**Gzip compression**: Output compression matches input by default. Use `--gz` to force compression of uncompressed input:
- `genome.fa` + `--gz` → `prefix.fa.gz` (compressed output from plain input)
- `genome.fa` (no --gz) → `prefix.fa` (plain output from plain input)
- `reads.fastq.gz` → `prefix.fastq.gz` and `prefix.snp.gz` (compressed output from compressed input)
- Paired-end output compression follows each mate input by default; use `--gz`
  to force compressed output.


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

#### Paired-End Synchronized Mode
- **Triggered by**: `--input-r1`, `--input-r2`, `--output-r1`, and `--output-r2`
- **Characteristics**: Validates synchronized R1/R2 chunks, mutates with deterministic per-read RNG, writes mate outputs in original order
- **Best for**: Mapper-ready paired-end mutation/damage without fragmentation
- **Rejects**: Fragmentation and `--num-mutations`
- **Threads**: `--threads 1` or `2` uses sequential PE mode; `--threads >= 3` uses 1 producer, `N-2` workers, and 1 ordered writer
- **Memory bound**: `--max-pending-chunks` limits queued PE chunks. Default is `min(--threads, 32)`. If `--threads > 32` and no explicit cap is set, SCAR warns to `stderr` and uses 32.


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
./scar --input file.fa --output out --num-mutations 100 --mutation-rate 0.01
```

**Issue**: No mutations in output (0 in SNP file)  
**Cause**: Mutation rate too low for small genome, or stochastic variation with random seed  
**Fix**: Increase mutation rate or use `--num-mutations` for guaranteed count

**Issue**: Unexpected mutation counts  
**Note**: For small genomes or low mutation rates, stochastic variation is expected. Use `--num-mutations` for deterministic counts.

**Issue**: Paired-end run fails with pair-name or record-count mismatch  
**Cause**: R1 and R2 are not synchronized, have different read counts, or have different normalized read names at the same record position.  
**Fix**: Re-pair or sort/filter R1 and R2 before running SCAR. Paired-end mode does not reorder reads; it preserves and validates the existing order.

**Issue**: Fragmentation requested with paired-end input  
**Cause**: Fragmentation would break mapper-ready R1/R2 synchronization.  
**Fix**: Use paired-end mode for mutation/damage-only output, or use the merge-for-fragmentation workflow when single-end merged output is intended.

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
./scar --input reads.fastq --output test1 --mutation-rate 0.001 --threads 8
./scar --input reads.fastq --output test2 --mutation-rate 0.0001 --threads 8
```

## Technical Details

### Threading Architecture
- **Producer thread**: Reads input file, distributes chunks to work queue
- **Worker threads**: Process mutations, write to thread-specific temporary files
- **Main thread**: Merges worker outputs into final file
- **Progress thread**: Reports status every 10 seconds
- **Paired-end mode**: Uses sequential processing for `--threads < 3`; otherwise uses a producer/worker/writer pipeline with ordered chunk output.

### Random Number Generation
- C++ `<random>` library with Mersenne Twister (`std::mt19937`)
- Seeded for reproducibility (`--seed` option)
- Paired-end mode derives deterministic per-read RNG seeds from `--seed`, pair index, normalized pair name, and mate label, so results are stable across PE thread counts.

### Quality Score Handling
- FASTQ quality scores preserved byte-for-byte
- Mutations only affect sequence line (line 2 of each record)
- Quality line (line 4) remains unchanged


## Citation

If you use this tool in published research, please cite this repository:

```
SCAR: Sequences + Controlled mutations + Ancient damage + fRagmentation
https://github.com/Madshartmann1/SCAR
```

## Contact

For questions or support, please:
- Open an issue on GitHub
- Email: madhar@dtu.dk

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

This project is currently in research stages.
