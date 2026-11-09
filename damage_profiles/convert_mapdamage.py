#!/usr/bin/env python3
"""
Convert mapDamage misincorporation.txt format to simple damage profile format.

MapDamage format has columns:
Chr End Std Pos A C G T Total G>A C>T A>G T>C [other substitutions]

This is the same conversion SCAR applies when --ancient-damage is given a
mapDamage directory or misincorporation.txt directly:
- Counts are summed over every Chr block and both strands (mapDamage reports
  both strands in read orientation).
- Raw C->T and G->A also contain reference divergence, sample variants and
  sequencing error. At each end the opposite substitution is the undamaged
  control, so only the excess is kept, clamped at 0:
      C->T damage = max(0, C->T - G->A)
      G->A damage = max(0, G->A - C->T)
  This gives 5' C->T and 3' G->A for double-stranded libraries, and C->T at
  both ends for single-stranded libraries.
- Damage cannot increase away from an end, so each bin is capped by the bins
  before it. This also drops deep bins that count the far end of short molecules.
  The cap uses a 3-bin running median, so a single outlier bin cannot zero the
  rest of the profile.

Output format (simplified):
end  position  from_base  to_base  frequency
5p   1         C          T        0.1534
5p   1         G          A        0.0000
...
"""

import sys
import argparse

def parse_mapdamage(filename):
    """Sum C, G, C>T and G>A counts per end and position over all Chr blocks and strands."""
    data = {
        '5p': {},  # 5' end data by position
        '3p': {}   # 3' end data by position
    }
    
    with open(filename) as f:
        for line in f:
            # Skip comments and header
            if line.startswith('#') or line.startswith('Chr'):
                continue
            
            fields = line.split()
            if len(fields) < 13:
                continue
            
            # Parse fields
            end = fields[1]  # '5p' or '3p'
            pos = int(fields[3])
            total = float(fields[8])
            
            if end not in data:
                raise ValueError(f"Invalid mapDamage end '{end}' in {filename}")
            if total <= 0:
                continue
            
            pos_data = data[end].setdefault(pos, {'C': 0.0, 'G': 0.0, 'C>T': 0.0, 'G>A': 0.0})
            pos_data['C'] += float(fields[5])
            pos_data['G'] += float(fields[6])
            pos_data['G>A'] += float(fields[9])
            pos_data['C>T'] += float(fields[10])
    
    return data

def cap_by_earlier_bins(values):
    """Cap each bin by the smallest earlier bin, using a 3-bin running median for the cap."""
    capped = []
    cap = 1.0
    for i, value in enumerate(values):
        capped.append(min(value, cap))
        smoothed = value
        if 0 < i < len(values) - 1:
            smoothed = sorted(values[i - 1:i + 2])[1]
        cap = min(cap, smoothed)
    return capped

def calculate_frequencies(data, max_positions=None):
    """Calculate excess damage frequencies from summed mapDamage counts."""
    results = []
    
    for end in ['5p', '3p']:
        positions = []
        c_to_t_excess = []
        g_to_a_excess = []
        for pos in sorted(data[end]):
            pos_data = data[end][pos]
            if pos_data['C'] <= 0 and pos_data['G'] <= 0:
                continue
            
            c_to_t = pos_data['C>T'] / pos_data['C'] if pos_data['C'] > 0 else 0.0
            g_to_a = pos_data['G>A'] / pos_data['G'] if pos_data['G'] > 0 else 0.0
            positions.append(pos)
            c_to_t_excess.append(max(0.0, c_to_t - g_to_a))
            g_to_a_excess.append(max(0.0, g_to_a - c_to_t))
        
        c_to_t_damage = cap_by_earlier_bins(c_to_t_excess)
        g_to_a_damage = cap_by_earlier_bins(g_to_a_excess)
        for pos, c_to_t, g_to_a in zip(positions, c_to_t_damage, g_to_a_damage):
            if max_positions and pos > max_positions:
                break
            results.append({'end': end, 'position': pos, 'from': 'C', 'to': 'T',
                            'frequency': c_to_t})
            results.append({'end': end, 'position': pos, 'from': 'G', 'to': 'A',
                            'frequency': g_to_a})
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Convert mapDamage output to simplified damage profile')
    parser.add_argument('input', help='mapDamage misincorporation.txt file')
    parser.add_argument('--output', '-o', help='Output file (default: stdout)')
    parser.add_argument('--max-positions', type=int, default=None,
                       help='Optional maximum positions from read end to cap at (default: all available)')
    
    args = parser.parse_args()
    
    # Parse mapDamage file
    data = parse_mapdamage(args.input)
    
    # Calculate frequencies
    freqs = calculate_frequencies(data, args.max_positions)
    
    # Output
    if args.output:
        outf = open(args.output, 'w')
    else:
        outf = sys.stdout
    
    # Header
    outf.write("# Damage profile converted from mapDamage output\n")
    outf.write(f"# Source: {args.input}\n")
    outf.write("# Excess over the opposite substitution, summed over strands and Chr blocks\n")
    outf.write("# Format: end position from_base to_base frequency\n")
    outf.write("end\tposition\tfrom\tto\tfrequency\n")
    
    # Data
    for entry in freqs:
        outf.write(f"{entry['end']}\t{entry['position']}\t"
                  f"{entry['from']}\t{entry['to']}\t"
                  f"{entry['frequency']:.6f}\n")
    
    if args.output:
        outf.close()
        print(f"✓ Wrote damage profile to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
