# Ancient DNA Damage Profiles

This directory contains damage profiles converted from mapDamage output for simulating ancient DNA damage patterns.

## Available Profiles

### 1. Double-Stranded DNA (Camelid - Alpaca)
**File**: `double_stranded_damage.txt`  
**Source**: ds_Camelid14-El-Olivar (Ancient alpaca from El-Olivar archaeological site)  
**Characteristics**:
- Moderate C→T deamination at 5' end (~4.1% at position 1)
- Moderate G→A deamination at 3' end (~4.2% at position 1)
- Symmetric damage pattern typical of double-stranded ancient DNA
- Well-preserved ancient mammal DNA from biodiversity extinction study
- **Biologically relevant**: Real data from extinct/endangered mammal studies

**Usage**:
```bash
./mutate_seq_v2 --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --threads 8 --seed 42
```

### 2. Single-Stranded DNA (Ccsp015)
**File**: `single_stranded_damage.txt`  
**Source**: ss_Ccsp015 ancient sample  
**Characteristics**:
- Strong C→T deamination at 5' end (~15.9% at position 1)
- Moderate G→A deamination at 3' end (~3.7% at position 1)  
- Asymmetric damage pattern typical of single-stranded ancient DNA
- Higher degradation suitable for ancient samples with poor preservation

**Usage**:
```bash
./mutate_seq_v2 --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/single_stranded_damage.txt \
    --threads 8 --seed 42
```

## Optional Background Mutation Rate

You can add a background mutation rate to simulate evolutionary divergence alongside damage:

```bash
./mutate_seq_v2 --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --background-rate 0.0001 \
    --threads 8 --seed 42
```

The background rate applies uniformly across all positions, while damage is position-dependent (enriched at read ends).

## Profile Format

The damage profile format is simple tab-delimited:

```
# Comment lines start with #
end     position    from    to    frequency
5p      1           C       T     0.168674
5p      2           C       T     0.056428
3p      1           G       A     0.019386
...
```

- **end**: `5p` (5' end) or `3p` (3' end)
- **position**: Distance from read end (1-based)
- **from/to**: Base substitution
- **frequency**: Damage probability (0.0-1.0)

## Converting mapDamage Output

Use the provided script to convert mapDamage `misincorporation.txt` files:

```bash
python3 convert_mapdamage.py path/to/misincorporation.txt --output custom_damage.txt
```

## Expected Mutation Patterns

With these profiles, you should observe:

1. **C→T enrichment at 5' end** (first ~10-15 bases)
2. **G→A enrichment at 3' end** (last ~10-15 bases)
3. **Rapid decay** of damage signal away from ends
4. **Background mutations** distributed uniformly (if --background-rate specified)

### Example Results (10k reads, 150bp, double-stranded profile, 0.0001 background):
- C→T: ~1013 (damage + background)
- G→A: ~999 (damage + background)
- Other substitutions: ~7-19 each (background only)

## References

- **LaBrana**: Olalde et al. (2014). Derived immune and ancestral pigmentation alleles in a 7,000-year-old Mesolithic European. Nature.
- **Ust_Ishim**: Fu et al. (2014). Genome sequence of a 45,000-year-old modern human from western Siberia. Nature.
- **mapDamage**: Jónsson et al. (2013). mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. Bioinformatics.
