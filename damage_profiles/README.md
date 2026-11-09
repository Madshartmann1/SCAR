# Ancient DNA Damage Profiles

This directory contains damage profiles converted from mapDamage output for simulating ancient DNA damage patterns.

## Available Profiles

### 1. Double-Stranded DNA (Camelid - Alpaca)
**File**: `double_stranded_damage.txt`  
**Source**: ds_Camelid14-El-Olivar (Ancient alpaca from El-Olivar archaeological site)  
**Characteristics**:
- Moderate C→T deamination at 5' end (~3.5% at position 1)
- Moderate G→A deamination at 3' end (~3.6% at position 1)
- Falls to 0 by position ~25 from either end
- Symmetric damage pattern typical of double-stranded ancient DNA
- Real data from extinct/endangered mammal studies

**Usage**:
```bash
./bin/scar --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --threads 8 --seed 27
```

### 2. Single-Stranded DNA (Ccsp015)
**File**: `single_stranded_damage.txt`  
**Source**: ss_Ccsp015 ancient sample  
**Characteristics**:
- Strong C→T deamination at 5' end (~15.3% at position 1)
- Strong C→T deamination at 3' end (~14.6% at position 1, ~9.4% at position 2)
- Low C→T level (~0.3%) through the rest of the read
- C→T at both ends, typical of single-stranded ancient DNA
- Higher degradation suitable for ancient samples with poor preservation

**Usage**:
```bash
./bin/scar --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/single_stranded_damage.txt \
    --threads 8 --seed 27
```

## Optional Additional Background Mutation Rate

You can add a background mutation rate to simulate evolutionary divergence alongside damage:

```bash
./bin/scar --input reads.fastq.gz --output damaged \
    --ancient-damage damage_profiles/double_stranded_damage.txt \
    --background-rate 0.0001 \
    --threads 8 --seed 27
```

The background rate applies uniformly across all positions, while damage is position-dependent (enriched at read ends).

## Profile Format

The damage profile format is simple tab-delimited:

```
# Comment lines start with #
end     position    from    to    frequency
5p      1           C       T     0.168674
5p      2           C       T     0.056278
3p      1           G       A     0.019386
...
```

- **end**: `5p` (5' end) or `3p` (3' end)
- **position**: Distance from read end (1-based)
- **from/to**: Base substitution
- **frequency**: Damage probability (0.0-1.0)

***Note on profile length vs. read length***: If an input read is longer than the maximum position defined in your damage profile, the tool will automatically carry forward the damage frequency of the *last defined position* and apply it to all remaining positions.

## Bespoke Damage Profiles - Converting mapDamage Output

SCAR can use mapDamage output directly. Pass either the mapDamage output
directory containing `misincorporation.txt` or the `misincorporation.txt` file
itself to `--ancient-damage`.

```bash
./bin/scar --input reads.fastq.gz --output damaged --ancient-damage path/to/mapdamage_results/
```

mapDamage reports raw C→T and G→A rates relative to the reference. These rates
also include reference divergence, the sample's own variants and sequencing
error, which the input reads already contain. SCAR keeps only the damage:

- Counts are summed over every reference sequence (`Chr`) block and both strands.
- At each end and position, the opposite substitution is the undamaged control:
  C→T damage = max(0, C→T − G→A), and G→A damage = max(0, G→A − C→T). This gives
  5′ C→T and 3′ G→A for double-stranded libraries, and C→T at both ends for
  single-stranded libraries.
- Damage cannot increase away from an end, so each position is capped by the
  positions before it (using a 3-position running median). This removes the
  rise at deep positions, where mapDamage counts the other end of short molecules.

The standalone converter applies the same conversion, for inspection/export workflows.
The bundled profiles in this directory were made with it:

```bash
damage_profiles/convert_mapdamage.py path/to/misincorporation.txt --output custom_damage.txt
```

A profile file written by hand is used as given, with no subtraction.

## Expected Mutation Patterns

With these profiles, you should observe:

1. **C→T enrichment at 5' end** (first ~10-15 bases)
2. **G→A enrichment at 3' end** (last ~10-15 bases) for the double-stranded profile,
   or **C→T enrichment at 3' end** for the single-stranded profile
3. **Rapid decay** of damage signal away from ends
4. **Background mutations** distributed uniformly (if --background-rate specified)

### Example Results (10k random reads, 150bp, double-stranded profile, 0.0001 background, seed 27):
- C→T: 138 (damage + background)
- G→A: 168 (damage + background)
- Other substitutions: 6-17 each (background only)

## References

- **Double-stranded example source**: López Mendoza P, González P, Westbury MV, Saghessi D, González Venanzi L, González BA, et al. (2025). Multi-proxy analysis of El Olivar camelids (1,090-1,440 cal AD): Evaluating the presence of llamas (Lama glama, Linnaeus 1758) in the Semiarid North of Chile before the arrival of the Inca. *PLoS ONE* 20(5): e0323497. https://doi.org/10.1371/journal.pone.0323497
- **Single-stranded example source**: Westbury MV, Hartmann S, Barlow A, Preick M, Hofreiter M. (2020). Hyena paleogenomes reveal a complex evolutionary history of cross-continental gene flow between spotted and cave hyena. *Science Advances* 6(11): eaay0456. https://doi.org/10.1126/sciadv.aay0456
- **mapDamage method**: Jónsson H, Ginolhac A, Schubert M, Johnson PLF, Orlando L. (2013). mapDamage2.0: fast approximate Bayesian estimates of ancient DNA damage parameters. *Bioinformatics* 29(13): 1682-1684. https://doi.org/10.1093/bioinformatics/btt193
