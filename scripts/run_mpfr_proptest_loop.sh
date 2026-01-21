#!/usr/bin/env bash
set -euo pipefail

runs="${1:-10}"
log_dir="${LOG_DIR:-proptest-runs}"
parallel="${PARALLEL:-8}"

mkdir -p "$log_dir"

cmd=(cargo test --features mpfr)
fail_file="$(mktemp)"

cleanup() {
    local code=$?
    jobs -pr | xargs -r kill || true
    rm -f "$fail_file"
    exit "$code"
}

trap cleanup INT TERM EXIT

run_one() {
    local i="$1"
    local log="$log_dir/run_${i}.log"
    echo "=== Run $i/$runs ==="
    if ! (PROPTEST_CASES=100000 CARGO_TARGET_DIR="$log_dir/target_${i}" "${cmd[@]}" 2>&1 | tee "$log"); then
        echo "$i" >> "$fail_file"
    fi
}

running=0
for i in $(seq 1 "$runs"); do
    run_one "$i" &
    running=$((running + 1))
    if [[ "$running" -ge "$parallel" ]]; then
        wait -n
        running=$((running - 1))
    fi
done

wait

if [[ -s "$fail_file" ]]; then
    while IFS= read -r i; do
        log="$log_dir/run_${i}.log"
        echo "FAILED on run $i"
        echo "Saved full log to $log"
        echo "---- failure excerpt ----"
        if command -v rg >/dev/null 2>&1; then
            rg -n -C 2 "minimal failing input|assertion failed|expected|left:|right:|proptest|panic" "$log" || true
        else
            grep -n -E -C 2 "minimal failing input|assertion failed|expected|left:|right:|proptest|panic" "$log" || true
        fi
        echo "---- last 200 lines ----"
        tail -n 200 "$log"
    done < "$fail_file"
    rm -f "$fail_file"
    exit 1
fi

rm -f "$fail_file"
echo "All $runs runs passed."
