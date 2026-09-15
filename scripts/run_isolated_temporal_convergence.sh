#!/usr/bin/env bash
set -euo pipefail

root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$root/build/isolated_temporal_convergence"
executable="$build_dir/isolated_temporal"
results="$build_dir/results.csv"
compiler="${CC:-cc}"
grid="${GRID:-48}"
total_time="${FINAL_TIME:-0.02}"

timesteps=(0.005 0.0025 0.00125 0.000625)
core_sources=()
for source in "$root"/src/*.c; do
    [[ "$(basename -- "$source")" == "main.c" ]] || core_sources+=("$source")
done

mkdir -p "$build_dir"
trap 'rm -f "$executable"' EXIT
printf '%s\n' 'N,steps,dt,L2_ux,L2_uy,L2_uz,L2_p' > "$results"

for dt in "${timesteps[@]}"; do
    steps="$(awk -v t="$total_time" -v d="$dt" 'BEGIN {printf "%.0f",t/d}')"
    "$compiler" -std=c11 -O3 -Wall -Wextra -D_DEFAULT_SOURCE \
        -I"$root/include" -DWIDTH="$grid" -DHEIGHT="$grid" -DDEPTH="$grid" \
        -DT="$total_time" -DSTEPS="$steps" \
        "$root/test/linear_convective_man.c" "${core_sources[@]}" -lm \
        -o "$executable"
    output="$($executable)"
    printf '%s,%s,%s,%s,%s,%s,%s\n' "$grid" "$steps" "$dt" \
        "$(awk '/^  L2 error u_x:/ {print $NF}' <<< "$output")" \
        "$(awk '/^  L2 error u_y:/ {print $NF}' <<< "$output")" \
        "$(awk '/^  L2 error u_z:/ {print $NF}' <<< "$output")" \
        "$(awk '/^  L2 error p:/ {print $NF}' <<< "$output")" >> "$results"
done

awk -F, 'NR == 1 {next} NR > 2 {
  printf "dt=%-9s rates: ux=%.3f uy=%.3f uz=%.3f\n", $3,
    log(px/$4)/log(2), log(py/$5)/log(2), log(pz/$6)/log(2)
} {px=$4; py=$5; pz=$6}' "$results"
printf 'Results written to %s\n' "$results"
