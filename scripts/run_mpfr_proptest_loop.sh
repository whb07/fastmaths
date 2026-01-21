#!/usr/bin/env bash
set -euo pipefail

runs="${1:-10}"
log_dir="${LOG_DIR:-proptest-runs}"

mkdir -p "$log_dir"

cmd=(cargo test --features mpfr)

for i in $(seq 1 "$runs"); do
    echo "=== Run $i/$runs ==="
    log="$log_dir/run_${i}.log"

    if ! (PROPTEST_CASES=100000 "${cmd[@]}" 2>&1 | tee "$log"); then
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
        exit 1
    fi
done

echo "All $runs runs passed."
