#!/usr/bin/env python3
"""
Convert mapDamage misincorporation.txt format to simple damage profile format.

MapDamage format has columns:
Chr End Std Pos A C G T Total G>A C>T A>G T>C [other substitutions]

We extract position-dependent frequencies for key damage patterns:
- 5' end: C->T (deamination at 5' end)
- 3' end: G->A (deamination at 3' end, reverse complement of C->T)

Output format (simplified):
end  position  from_base  to_base  frequency
5p   1         C          T        0.1687
5p   2         C          T        0.0564
...
"""

import sys
import argparse

def parse_mapdamage(filename):
    """Parse mapDamage misincorporation.txt file."""
    data = {
        '5p': {},  # 5' end data by position
        '3p': {}   # 3' end data by position
    }
    
    with open(filename) as f:
        for line in f:
            # Skip comments and header
            if line.startswith('#') or line.startswith('Chr'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 20:
                continue
            
            # Parse fields
            chr_id = fields[0]
            end = fields[1]  # '5p' or '3p'
            strand = fields[2]  # '+' or '-'
            pos = int(fields[3])
            
            # Base counts
            a_count = int(fields[4])
            c_count = int(fields[5])
            g_count = int(fields[6])
            t_count = int(fields[7])
            total = int(fields[8])
            
            # Substitution counts
            g_to_a = int(fields[9])
            c_to_t = int(fields[10])
            a_to_g = int(fields[11])
            t_to_c = int(fields[12])
            
            # Store only positive strand data (mapDamage reports both strands)
            if strand == '+' and total > 0:
                if end not in data:
                    data[end] = {}
                if pos not in data[end]:
                    data[end][pos] = {
                        'C': c_count,
                        'G': g_count,
                        'C>T': c_to_t,
                        'G>A': g_to_a,
                        'total': total
                    }
    
    return data

def calculate_frequencies(data, max_positions=25):
    """Calculate damage frequencies from mapDamage data."""
    results = []
    
    for end in ['5p', '3p']:
        if end not in data:
            continue
        
        for pos in range(1, max_positions + 1):
            if pos not in data[end]:
                continue
            
            pos_data = data[end][pos]
            
            # C->T frequency (key aDNA damage signal)
            if pos_data['C'] > 0:
                c_to_t_freq = pos_data['C>T'] / pos_data['C']
                results.append({
                    'end': end,
                    'position': pos,
                    'from': 'C',
                    'to': 'T',
                    'frequency': c_to_t_freq
                })
            
            # G->A frequency (reverse complement)
            if pos_data['G'] > 0:
                g_to_a_freq = pos_data['G>A'] / pos_data['G']
                results.append({
                    'end': end,
                    'position': pos,
                    'from': 'G',
                    'to': 'A',
                    'frequency': g_to_a_freq
                })
    
    return results

def main():
    parser = argparse.ArgumentParser(description='Convert mapDamage output to simplified damage profile')
    parser.add_argument('input', help='mapDamage misincorporation.txt file')
    parser.add_argument('--output', '-o', help='Output file (default: stdout)')
    parser.add_argument('--max-positions', type=int, default=25, 
                       help='Maximum positions from read end (default: 25)')
    
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
