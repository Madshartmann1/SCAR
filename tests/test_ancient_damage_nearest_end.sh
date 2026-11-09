#!/usr/bin/env bash
set -euo pipefail

script_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
scar_directory="$(cd "${script_directory}/.." && pwd)"
scar_executable="${scar_directory}/bin/scar"
temporary_directory="$(mktemp -d "${TMPDIR:-/tmp}/scar_nearest_end_tests.XXXXXXXX")"
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

run_scar() {
    local output_prefix="$1"
    shift
    "${scar_executable}" "$@" --output "${output_prefix}" \
        > "${output_prefix}.stdout" 2> "${output_prefix}.stderr"
}

five_prime_input="${temporary_directory}/five_prime.fastq"
five_prime_profile="${temporary_directory}/five_prime.profile.tsv"
five_prime_output="${temporary_directory}/five_prime"
{
    printf '@five_prime\nG'
    printf 'A%.0s' {1..149}
    printf '\n+\n'
    printf 'I%.0s' {1..150}
    printf '\n'
} > "${five_prime_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t1\tG\tA\t0\n'
    printf '3p\t70\tG\tA\t1\n'
} > "${five_prime_profile}"
run_scar "${five_prime_output}" --input "${five_prime_input}" \
    --ancient-damage "${five_prime_profile}" --threads 2 --seed 27
assert_equals "G" "$(sed -n '2p' "${five_prime_output}.fastq" | cut -c1)" \
    "5-prime base does not receive a distant 3-prime trial"

three_prime_input="${temporary_directory}/three_prime.fastq"
three_prime_profile="${temporary_directory}/three_prime.profile.tsv"
three_prime_output="${temporary_directory}/three_prime"
{
    printf '@three_prime\n'
    printf 'A%.0s' {1..149}
    printf 'C\n+\n'
    printf 'I%.0s' {1..150}
    printf '\n'
} > "${three_prime_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t70\tC\tT\t1\n'
    printf '3p\t1\tC\tT\t0\n'
} > "${three_prime_profile}"
run_scar "${three_prime_output}" --input "${three_prime_input}" \
    --ancient-damage "${three_prime_profile}" --threads 2 --seed 27
assert_equals "C" "$(sed -n '2p' "${three_prime_output}.fastq" | cut -c150)" \
    "3-prime base does not receive a distant 5-prime trial"

midpoint_input="${temporary_directory}/midpoint.fastq"
midpoint_profile="${temporary_directory}/midpoint.profile.tsv"
midpoint_output="${temporary_directory}/midpoint"
for read_number in $(seq 1 2000); do
    printf '@midpoint_%d\nAACAA\n+\nIIIII\n' "${read_number}"
done > "${midpoint_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t3\tC\tT\t0\n'
    printf '3p\t3\tC\tT\t1\n'
} > "${midpoint_profile}"
run_scar "${midpoint_output}" --input "${midpoint_input}" \
    --ancient-damage "${midpoint_profile}" --threads 2 --seed 27
midpoint_count="$(wc -l < "${midpoint_output}.snp" | tr -d ' ')"
assert_integer_between "${midpoint_count}" 850 1150 \
    "midpoint averages end rates and performs one trial"

carry_input="${temporary_directory}/carry.fastq"
carry_profile="${temporary_directory}/carry.profile.tsv"
carry_output="${temporary_directory}/carry"
printf '@carry\nAAGAAAA\n+\nIIIIIII\n' > "${carry_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t2\tG\tA\t1\n'
    printf '3p\t2\tG\tA\t0\n'
} > "${carry_profile}"
run_scar "${carry_output}" --input "${carry_input}" \
    --ancient-damage "${carry_profile}" --threads 2 --seed 27
assert_equals "A" "$(sed -n '2p' "${carry_output}.fastq" | cut -c3)" \
    "final profile bin carries forwards from the nearest end"

native_profile="${temporary_directory}/native.profile.tsv"
mapdamage_profile="${temporary_directory}/misincorporation.txt"
profile_input="${temporary_directory}/profile_input.fastq"
native_output="${temporary_directory}/native"
mapdamage_output="${temporary_directory}/mapdamage"
printf '@profile_input\nGAAAA\n+\nIIIII\n' > "${profile_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t1\tG\tA\t1\n'
} > "${native_profile}"
{
    printf 'Chr\tEnd\tStd\tPos\tA\tC\tG\tT\tTotal\tG>A\tC>T\tA>G\tT>C\n'
    printf '*\t5p\t+\t1\t0\t0\t1\t0\t1\t1\t0\t0\t0\n'
} > "${mapdamage_profile}"
run_scar "${native_output}" --input "${profile_input}" \
    --ancient-damage "${native_profile}" --threads 2 --seed 27
run_scar "${mapdamage_output}" --input "${profile_input}" \
    --ancient-damage "${mapdamage_profile}" --threads 2 --seed 27
assert_equals "$(sed -n '2p' "${native_output}.fastq")" \
    "$(sed -n '2p' "${mapdamage_output}.fastq")" \
    "native and mapDamage profiles have equivalent nearest-end behaviour"

thread_input="${temporary_directory}/thread_input.fastq"
thread_one_output="${temporary_directory}/thread_one"
thread_two_output="${temporary_directory}/thread_two"
for read_number in $(seq 1 20); do
    printf '@thread_%d\nGAAAA\n+\nIIIII\n' "${read_number}"
done > "${thread_input}"
run_scar "${thread_one_output}" --input "${thread_input}" \
    --ancient-damage "${native_profile}" --threads 1 --seed 27
run_scar "${thread_two_output}" --input "${thread_input}" \
    --ancient-damage "${native_profile}" --threads 2 --seed 27
assert_equals "$(awk 'NR % 4 == 2' "${thread_one_output}.fastq")" \
    "$(awk 'NR % 4 == 2' "${thread_two_output}.fastq")" \
    "single-threaded and streaming paths use nearest-end damage"

exact_input="${temporary_directory}/exact.fastq"
exact_profile="${temporary_directory}/exact.profile.tsv"
exact_output="${temporary_directory}/exact"
printf '@exact\nG\n+\nI\n' > "${exact_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t1\tA\tC\t1\n'
    printf '5p\t1\tC\tT\t1\n'
    printf '5p\t1\tT\tA\t1\n'
    printf '3p\t1\tA\tC\t1\n'
    printf '3p\t1\tC\tT\t1\n'
    printf '3p\t1\tT\tA\t1\n'
} > "${exact_profile}"
run_scar "${exact_output}" --input "${exact_input}" --num-mutations 1 \
    --ancient-damage "${exact_profile}" --threads 1 --seed 27
assert_equals "2" "$(wc -l < "${exact_output}.snp" | tr -d ' ')" \
    "exact-count mode applies one ordinary and one ancient-damage mutation"

rate_input="${temporary_directory}/rate.fastq"
rate_profile="${temporary_directory}/rate.profile.tsv"
rate_output="${temporary_directory}/rate"
for read_number in $(seq 1 5000); do
    printf '@rate_%d\nCAAAA\n+\nIIIII\n' "${read_number}"
done > "${rate_input}"
{
    printf 'end\tposition\tfrom\tto\tfrequency\n'
    printf '5p\t1\tC\tT\t0.2\n'
    printf '3p\t1\tC\tT\t0\n'
} > "${rate_profile}"
run_scar "${rate_output}" --input "${rate_input}" \
    --ancient-damage "${rate_profile}" --threads 2 --seed 27
rate_count="$(wc -l < "${rate_output}.snp" | tr -d ' ')"
assert_integer_between "${rate_count}" 850 1150 \
    "large fixed-seed simulation recovers the requested 20 percent rate"

printf '[PASS] all ancient-damage nearest-end tests completed\n'
