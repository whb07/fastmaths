#!/usr/bin/env bash
set -euo pipefail

runs="${1:-10}"
log_dir="${LOG_DIR:-proptest-runs}"
parallel="${PARALLEL:-8}"

mkdir -p "$log_dir"

cmd=(cargo test --release --features mpfr)
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
    local target_dir=""
    local -a env_vars=("PROPTEST_CASES=100000")
    if [[ "$runs" -gt 1 ]]; then
        target_dir="$log_dir/target_${i}"
        env_vars+=("CARGO_TARGET_DIR=$target_dir")
    fi
    echo "=== Run $i/$runs ==="
    if ! (env "${env_vars[@]}" "${cmd[@]}" 2>&1 | tee "$log"); then
        echo "$i" >> "$fail_file"
        # Stop any in-flight runs as soon as one fails.
        jobs -pr | xargs -r kill || true
    else
        if [[ -n "$target_dir" ]]; then
            rm -rf "$target_dir"
        fi
    fi
}

running=0
for i in $(seq 1 "$runs"); do
    if [[ -s "$fail_file" ]]; then
        break
    fi
    run_one "$i" &
    running=$((running + 1))
    if [[ "$running" -ge "$parallel" ]]; then
        wait -n
        running=$((running - 1))
        if [[ -s "$fail_file" ]]; then
            break
        fi
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
