# mutate_seq_v2 Quick Reference Card

## Common Use Cases

### 1. Your Pipeline: Bear/Elephant/Rhino Variants (N/H/L)

```bash
# For each species:
SPECIES="bear_black"
REF="${SPECIES}.fa"
SEED=12345

# N (Normal - no mutations)
cp ${REF} ${SPECIES}_N.fa

# H (High - 0.1% divergence)
./mutate_seq_v2 --input ${REF} --output ${SPECIES}_H \
    --mutation-rate 0.001 --threads 12 --seed ${SEED}

# L (Low - 0.01% divergence)
./mutate_seq_v2 --input ${REF} --output ${SPECIES}_L \
    --mutation-rate 0.0001 --threads 12 --seed $((SEED + 1))
```

### 2. Ancient DNA Damage Simulation

```bash
# Create C→T biased spectrum file (ancient.txt):
cat > ancient_damage.txt << 'EOF'
C       T       10.0
G       A       10.0
A       C       1.0
A       G       1.0
A       T       1.0
C       A       1.0
C       G       1.0
G       C       1.0
G       T       1.0
T       A       1.0
T       C       1.0
T       G       1.0
EOF

./mutate_seq_v2 --input ancient_sample.fa --output damaged \
    --mutation-spectrum ancient_damage.txt --mutation-rate 0.02 --seed 42
```

### 3. Batch Processing 200 Individuals

```bash
# Create batch script
cat > batch_mutate.sh << 'EOF'
#!/bin/bash
RATE=$1
THREADS=12

while read SAMPLE; do
    echo "Processing ${SAMPLE}..."
    ./mutate_seq_v2 --input ${SAMPLE}.fa --output ${SAMPLE}_mut \
        --mutation-rate ${RATE} --threads ${THREADS} --seed ${RANDOM}
done < sample_list.txt
EOF

chmod +x batch_mutate.sh
./batch_mutate.sh 0.001
```

## Quick Mode Selection Guide

| Goal | Mode | Command Pattern |
|------|------|-----------------|
| Simple uniform mutations | Flat rate | `--mutation-rate 0.001` |
| Exactly N mutations | Flat number | `--num-mutations 1000` |
| Control Ts/Tv bias | Ts/Tv ratio | `--mutation-rate 0.001 --ts-tv-ratio 2.0` |
| Ancient DNA damage | Spectrum | `--mutation-spectrum damage.txt --mutation-rate 0.02` |
| Per-substitution control | Matrix | `--mutation-matrix rates.txt --mutation-rate 0.001` |

## Performance Tuning

| File Type | Size | Threads | Expected Time |
|-----------|------|---------|---------------|
| Reference genome | 3 GB | 1 | ~30 sec |
| FASTQ reads | 50 GB | 8 | ~5-10 min |
| Gzipped FASTQ | 15 GB | 12 | ~8-15 min |

**Rule of thumb**: Use 1 thread per 5-10 GB of gzipped data.

## File Format Matrix

| Input Format | Output Format | Compression | Mode |
|--------------|---------------|-------------|------|
| genome.fa | genome.fa | None | In-memory |
| genome.fa.gz | genome.fa.gz | gzip | In-memory |
| reads.fastq | reads.fastq | None | Streaming* |
| reads.fastq.gz | reads.fastq.gz | gzip | Streaming* |

*Streaming requires ≥2 threads. Single thread uses in-memory.

## Mutation Rate Guide

| Scenario | Typical Rate | Command |
|----------|--------------|---------|
| Intraspecies variation | 0.001 (0.1%) | `--mutation-rate 0.001` |
| Population divergence | 0.0001 (0.01%) | `--mutation-rate 0.0001` |
| Ancient DNA damage | 0.02-0.05 (2-5%) | `--mutation-rate 0.03` |
| Interspecies divergence | 0.05-0.1 (5-10%) | `--mutation-rate 0.075` |

## Troubleshooting Checklist

### No mutations in output?
- [ ] Check mutation rate (is it too low for genome size?)
- [ ] Try `--num-mutations 100` to guarantee mutations
- [ ] Verify SNP file is not empty: `wc -l output.snp`

### Memory issues?
- [ ] Use streaming mode: FASTQ with `--threads 4`
- [ ] Or use `--force-streaming` for FASTA

### Slow processing?
- [ ] Increase threads: `--threads 8` or `--threads 12`
- [ ] Check if gzip is bottleneck (decompress first)
- [ ] Verify disk I/O isn't saturated

### "Multiple mutation modes" error?
- [ ] Use ONLY ONE mode (don't mix `--num-mutations` with `--mutation-rate`)
- [ ] Exception: `--mutation-rate` + `--ts-tv-ratio` is allowed

## Output File Examples

### SNP File Format (seqtk compatible)
```
seq_1       45      A       G
seq_1       127     C       T
seq_2       89      G       A
read_100    23      T       C
```

### Integration with seqtk
```bash
# Apply SNP file to original sequence
seqtk mutfa original.fa mutations.snp > mutated.fa

# Or use our tool directly (faster, includes multi-threading)
./mutate_seq_v2 --input original.fa --output mutated --mutation-rate 0.001
```

## One-Liners for Common Tasks

```bash
# Quick test with small file
./mutate_seq_v2 --input test.fa --output out --mutation-rate 0.01 --seed 42

# Production run with all features
./mutate_seq_v2 --input data.fastq.gz --output results \
    --mutation-rate 0.001 --ts-tv-ratio 2.0 --threads 12 --seed 12345

# Check output
wc -l results.snp                    # Count mutations
head -5 results.snp                  # Preview SNPs
zcat results.fastq.gz | head -8      # View mutated reads
```

## Notes for Your 200-Sample Pipeline

1. **Reproducibility**: Use sequential seeds (`--seed ${SAMPLE_ID}`)
2. **Parallelization**: Process multiple samples in parallel (not threads per sample)
3. **Storage**: Gzipped outputs save ~70% disk space
4. **Validation**: Check SNP file counts match expected rates
5. **Tracking**: SNP files are your "ground truth" for validation

```bash
# Example parallel batch (4 samples at once)
parallel -j 4 './mutate_seq_v2 --input {}.fa --output {}_H --mutation-rate 0.001 --threads 4 --seed {}' ::: $(seq 1 200)
```

## Version Info

- **Tool**: mutate_seq_v2
- **Language**: C++11
- **Dependencies**: zlib, pthread
- **Compilation**: `make -f Makefile_v2`
- **Help**: `./mutate_seq_v2 --help`
