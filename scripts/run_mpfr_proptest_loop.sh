#!/usr/bin/env bash
set -euo pipefail

runs="${1:-10}"
log_dir="${LOG_DIR:-proptest-runs}"
parallel="${PARALLEL:-8}"
proptest_cases="${PROPTEST_CASES:-100000}"

mkdir -p "$log_dir"

mode="${MODE:-release}"
build_cmd=(cargo test --features mpfr --no-run --message-format=json)
if [[ "$mode" == "release" ]]; then
    build_cmd=(cargo test --release --features mpfr --no-run --message-format=json)
fi
fail_file="$(mktemp)"
build_log="$log_dir/build.json"

cleanup() {
    local code=$?
    jobs -pr | xargs -r kill || true
    rm -f "$fail_file"
    exit "$code"
}

trap cleanup INT TERM EXIT

echo "Building test binaries (mode=$mode)..."
"${build_cmd[@]}" | tee "$build_log"

mapfile -t test_bins < <(
    python3 - "$build_log" <<'PY'
import json
import sys

path = sys.argv[1]
bins = []
with open(path, "r", encoding="utf-8") as f:
    for line in f:
        try:
            obj = json.loads(line)
        except json.JSONDecodeError:
            continue
        if obj.get("reason") != "compiler-artifact":
            continue
        prof = obj.get("profile") or {}
        if not prof.get("test"):
            continue
        exe = obj.get("executable")
        if exe:
            bins.append(exe)
for exe in bins:
    print(exe)
PY
)

if [[ "${#test_bins[@]}" -eq 0 ]]; then
    echo "No test binaries found in $build_log" >&2
    exit 1
fi

run_one() {
    local i="$1"
    local log="$log_dir/run_${i}.log"
    local -a env_vars=("PROPTEST_CASES=$proptest_cases")
    echo "=== Run $i/$runs ==="
    if ! (
        set +e
        status=0
        for exe in "${test_bins[@]}"; do
            echo ">>> $exe"
            env "${env_vars[@]}" "$exe"
            rc=$?
            if [[ "$rc" -ne 0 ]]; then
                echo ">>> $exe exited with $rc"
                status="$rc"
                break
            fi
        done
        exit "$status"
    ) 2>&1 | tee "$log"; then
        echo "$i" >> "$fail_file"
        # Stop any in-flight runs as soon as one fails.
        jobs -pr | xargs -r kill || true
    elif rg -n -q "test result: FAILED|failures:|panicked at|error: test failed" "$log" 2>/dev/null; then
        echo "$i" >> "$fail_file"
        jobs -pr | xargs -r kill || true
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
