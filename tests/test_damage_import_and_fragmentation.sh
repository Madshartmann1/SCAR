#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scar_directory="$(cd "${script_directory}/.." && pwd)"
scar_executable="${SCAR_EXECUTABLE:-${scar_directory}/bin/scar}"
converter="${scar_directory}/damage_profiles/convert_mapdamage.py"
temporary_directory="$(mktemp -d "${TMPDIR:-/tmp}/scar_import_fragment_tests.XXXXXXXX")"
trap 'rm -rf "${temporary_directory}"' EXIT

assert_equals() {
    local expected="$1"
    local observed="$2"
    local description="$3"
    if [[ "${expected}" != "${observed}" ]]; then
        printf '[FAIL] %s: expected "%s", observed "%s"\n' \
            "${description}" "${expected}" "${observed}" >&2
        exit 1
    fi
    printf '[PASS] %s\n' "${description}"
}

assert_integer_between() {
    local observed="$1"
    local minimum="$2"
    local maximum="$3"
    local description="$4"
    if (( observed < minimum || observed > maximum )); then
        printf '[FAIL] %s: expected %d-%d, observed %d\n' \
            "${description}" "${minimum}" "${maximum}" "${observed}" >&2
        exit 1
    fi
    printf '[PASS] %s\n' "${description}"
}

assert_fails_with() {
    local expected_message="$1"
    local description="$2"
    shift 2
    local output
    if output="$("${scar_executable}" "$@" 2>&1)"; then
        printf '[FAIL] %s: command succeeded\n' "${description}" >&2
        exit 1
    fi
    if [[ "${output}" != *"${expected_message}"* ]]; then
        printf '[FAIL] %s: missing "%s" in "%s"\n' \
            "${description}" "${expected_message}" "${output}" >&2
        exit 1
    fi
    printf '[PASS] %s\n' "${description}"
}

run_scar() {
    local output_prefix="$1"
    shift
    "${scar_executable}" "$@" --output "${output_prefix}" \
        > "${output_prefix}.stdout" 2> "${output_prefix}.stderr"
}

write_reads() {
    local output_file="$1"
    local count="$2"
    local sequence="$3"
    local quality
    quality="$(printf 'I%.0s' $(seq 1 "${#sequence}"))"
    for read_number in $(seq 1 "${count}"); do
        printf '@read_%d\n%s\n+\n%s\n' "${read_number}" "${sequence}" "${quality}"
    done > "${output_file}"
}

# mapDamage row: chr end strand position C C>T G G>A
mapdamage_row() {
    printf '%s\t%s\t%s\t%s\t100\t%s\t%s\t100\t%s\t%s\t%s\t0\t0\n' \
        "$1" "$2" "$3" "$4" "$5" "$7" "$(( 200 + $5 + $7 ))" "$8" "$6"
}

mapdamage_header() {
    printf 'Chr\tEnd\tStd\tPos\tA\tC\tG\tT\tTotal\tG>A\tC>T\tA>G\tT>C\n'
}

# Raw C>T matched by G>A at the same end is background, not damage
background_directory="${temporary_directory}/background"
mkdir -p "${background_directory}"
{
    mapdamage_header
    mapdamage_row '*' 5p + 1 100 100 100 100
} > "${background_directory}/misincorporation.txt"
write_reads "${temporary_directory}/c_read.fastq" 1 "CAAAA"
run_scar "${temporary_directory}/background_out" \
    --input "${temporary_directory}/c_read.fastq" \
    --ancient-damage "${background_directory}" --threads 1 --seed 27
assert_equals "CAAAA" "$(sed -n '2p' "${temporary_directory}/background_out.fastq")" \
    "matching G>A background cancels mapDamage C>T"

# Counts are summed over Chr blocks and strands, not taken from the last block
blocks_directory="${temporary_directory}/blocks"
mkdir -p "${blocks_directory}"
{
    mapdamage_header
    mapdamage_row chr1 5p + 1 100 100 100 0
    mapdamage_row chr1 5p - 1 100 100 100 0
    mapdamage_row chrM 5p + 1 100 0 100 0
    mapdamage_row chrM 5p - 1 100 0 100 0
} > "${blocks_directory}/misincorporation.txt"
write_reads "${temporary_directory}/blocks_reads.fastq" 2000 "CAAAA"
run_scar "${temporary_directory}/blocks_out" \
    --input "${temporary_directory}/blocks_reads.fastq" \
    --ancient-damage "${blocks_directory}" --threads 1 --seed 27
assert_integer_between "$(wc -l < "${temporary_directory}/blocks_out.snp" | tr -d ' ')" \
    850 1150 "mapDamage Chr blocks and strands are summed"

# Damage cannot rise away from an end, but one outlier bin does not zero the rest
shape_directory="${temporary_directory}/shape"
mkdir -p "${shape_directory}"
{
    mapdamage_header
    mapdamage_row '*' 5p + 1 100 100 100 0
    mapdamage_row '*' 5p + 2 100 100 100 0
    mapdamage_row '*' 5p + 3 100 0 100 0
    mapdamage_row '*' 5p + 4 100 100 100 0
    mapdamage_row '*' 5p + 5 100 100 100 0
} > "${shape_directory}/misincorporation.txt"
write_reads "${temporary_directory}/shape_read.fastq" 1 "CCCCCCCCCCCC"
run_scar "${temporary_directory}/shape_out" \
    --input "${temporary_directory}/shape_read.fastq" \
    --ancient-damage "${shape_directory}" --threads 1 --seed 27
assert_equals "TTCTTTCCCCCC" "$(sed -n '2p' "${temporary_directory}/shape_out.fastq")" \
    "single outlier bin keeps the rest of the profile"

rise_directory="${temporary_directory}/rise"
mkdir -p "${rise_directory}"
{
    mapdamage_header
    mapdamage_row '*' 5p + 1 100 100 100 0
    mapdamage_row '*' 5p + 2 100 0 100 0
    mapdamage_row '*' 5p + 3 100 0 100 0
    mapdamage_row '*' 5p + 4 100 100 100 0
    mapdamage_row '*' 5p + 5 100 100 100 0
    mapdamage_row '*' 5p + 6 100 100 100 0
} > "${rise_directory}/misincorporation.txt"
write_reads "${temporary_directory}/rise_read.fastq" 1 "CCCCCCCCCCCCCC"
run_scar "${temporary_directory}/rise_out" \
    --input "${temporary_directory}/rise_read.fastq" \
    --ancient-damage "${rise_directory}" --threads 1 --seed 27
assert_equals "TCCCCCCCCCCCCC" "$(sed -n '2p' "${temporary_directory}/rise_out.fastq")" \
    "deep bins cannot rise above earlier bins"

# The converter writes the same profile SCAR builds from mapDamage in memory
for library in double single; do
    python3 "${converter}" \
        "${scar_directory}/damage_profiles/${library}_stranded/misincorporation.txt" \
        --output "${temporary_directory}/${library}.profile.tsv" 2> /dev/null
    write_reads "${temporary_directory}/converter_reads.fastq" 2000 \
        "$(printf 'ACGT%.0s' $(seq 1 25))"
    run_scar "${temporary_directory}/${library}_native" \
        --input "${temporary_directory}/converter_reads.fastq" \
        --ancient-damage "${temporary_directory}/${library}.profile.tsv" --threads 1 --seed 27
    run_scar "${temporary_directory}/${library}_mapdamage" \
        --input "${temporary_directory}/converter_reads.fastq" \
        --ancient-damage "${scar_directory}/damage_profiles/${library}_stranded" \
        --threads 1 --seed 27
    assert_equals "$(cat "${temporary_directory}/${library}_native.snp")" \
        "$(cat "${temporary_directory}/${library}_mapdamage.snp")" \
        "converter matches in-memory mapDamage import (${library}-stranded)"
    assert_equals "$(grep -v '^#' "${temporary_directory}/${library}.profile.tsv")" \
        "$(grep -v '^#' "${scar_directory}/damage_profiles/${library}_stranded_damage.txt")" \
        "bundled ${library}-stranded profile is up to date"
done

# Reads shorter than a sampled length are fragmented from the truncated distribution
long_reads="${temporary_directory}/long_reads.fastq"
write_reads "${long_reads}" 4000 "$(printf 'ACGTACGTAC%.0s' $(seq 1 15))"
run_scar "${temporary_directory}/truncated" --input "${long_reads}" \
    --mutation-rate 0.001 --fragment-distribution lognormal \
    --mean-length 160 --sd-length 60 --threads 1 --seed 27
assert_equals "0" \
    "$(awk '/Reads discarded/ {print $NF}' "${temporary_directory}/truncated.stdout")" \
    "no read is discarded when a shorter fragment length is possible"
first_fragment_mean="$(awk 'NR % 4 == 1 && /_frag1_/ {getline; total += length($0); count++}
    END {printf "%d", total / count}' "${temporary_directory}/truncated.fastq")"
assert_integer_between "${first_fragment_mean}" 112 118 \
    "first fragments follow the distribution truncated at the read length"

empirical_file="${temporary_directory}/empirical.txt"
printf '50\n60\n' > "${empirical_file}"
write_reads "${temporary_directory}/short_reads.fastq" 10 "$(printf 'A%.0s' $(seq 1 40))"
run_scar "${temporary_directory}/too_short" \
    --input "${temporary_directory}/short_reads.fastq" --mutation-rate 0.001 \
    --fragment-distribution empirical --fragment-distribution-file "${empirical_file}" \
    --threads 1 --seed 27
assert_equals "10" \
    "$(awk '/Reads discarded/ {print $NF}' "${temporary_directory}/too_short.stdout")" \
    "reads shorter than every allowed length are discarded"

# --max-fragment-length caps every sampled length
run_scar "${temporary_directory}/capped" --input "${long_reads}" \
    --mutation-rate 0.001 --fragment-distribution lognormal \
    --mean-length 100 --sd-length 50 --max-fragment-length 60 --threads 1 --seed 27
assert_integer_between "$(awk 'NR % 4 == 2 && length($0) > max {max = length($0)}
    END {print max}' "${temporary_directory}/capped.fastq")" 1 60 \
    "--max-fragment-length caps fragment lengths"

# --min-fragment-length is a true minimum for sampled fragments and remainders
run_scar "${temporary_directory}/floored" --input "${long_reads}" \
    --mutation-rate 0.001 --fragment-distribution exponential \
    --mean-length 40 --min-fragment-length 30 --threads 1 --seed 27
assert_integer_between "$(awk 'NR % 4 == 2 && (min == "" || length($0) < min) {min = length($0)}
    END {print min}' "${temporary_directory}/floored.fastq")" 30 150 \
    "--min-fragment-length is the shortest output fragment"

# --drop-remainders keeps only fragments with a sampled length
run_scar "${temporary_directory}/with_remainders" --input "${long_reads}" \
    --mutation-rate 0.001 --fragment-length 40 --threads 1 --seed 27
run_scar "${temporary_directory}/without_remainders" --input "${long_reads}" \
    --mutation-rate 0.001 --fragment-length 40 --drop-remainders --threads 1 --seed 27
assert_equals "30 40" \
    "$(awk 'NR % 4 == 2 {print length($0)}' "${temporary_directory}/with_remainders.fastq" |
        sort -un | tr '\n' ' ' | sed 's/ $//')" \
    "remainders are kept by default"
assert_equals "40" \
    "$(awk 'NR % 4 == 2 {print length($0)}' "${temporary_directory}/without_remainders.fastq" |
        sort -un | tr '\n' ' ' | sed 's/ $//')" \
    "--drop-remainders discards the leftover end of each read"

# Argument checks
assert_fails_with "Unknown option or unexpected argument: --max-fragment-lenght" \
    "misspelled options are rejected" \
    --input "${long_reads}" --output "${temporary_directory}/typo" -r 0.001 -f 40 \
    --max-fragment-lenght 60
assert_fails_with "--seed requires a value" "options missing a value are rejected" \
    --input "${long_reads}" --output "${temporary_directory}/missing" -r 0.001 --seed
assert_fails_with "--max-fragment-length requires fragmentation" \
    "--max-fragment-length without fragmentation is rejected" \
    --input "${long_reads}" --output "${temporary_directory}/no_fragmentation" -r 0.001 \
    --max-fragment-length 60
assert_fails_with "No fragment length between 20 and 40 bp is possible" \
    "--max-fragment-length below the static length is rejected" \
    --input "${long_reads}" --output "${temporary_directory}/static_cap" -r 0.001 \
    -f 50 --max-fragment-length 40
assert_fails_with "No fragment length of at least 20 bp is possible" \
    "static length below --min-fragment-length is rejected" \
    --input "${long_reads}" --output "${temporary_directory}/static_floor" -r 0.001 \
    -f 15

printf '[PASS] all damage import and fragmentation tests completed\n'
