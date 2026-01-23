# Repository Guidelines

## Project Goals (glibc Parity + Speed)

This crate aims to implement **fast, correct** replacements for glibc `libm` routines like `exp`, `log/ln`, `sin`, and `cos`. Changes should target:

- **Behavior parity with glibc** across edge cases (NaNs, `±0.0`, `±∞`, subnormals, huge arguments, and rounding where applicable).
- **Accuracy Target:** Math functions should aim for an error of < 1.0 ULP relative to the true value (verified using MPFR as the reference).
- **Performance parity or better** than the system glibc implementation on supported targets (ideally faster).
- **No System Libm:** Deferring to or use of the system's `libm` (e.g., calling `f64::sin` or similar within the implementation) is not allowed. The crate must provide its own self-contained implementations.

## Project Structure & Module Organization

- `Cargo.toml` / `Cargo.lock`: Rust crate metadata (edition **2024**) and dependency lockfile.
- `src/lib.rs`: Crate entry point and basic unit tests.
- `src/maths.rs`: Low-level `f64` math routines (table-driven/bit-level code intended to be usable in `no_std` contexts).
- `benches/`: Criterion benchmarks (e.g., `benches/math_bench.rs`).
- `MATHS_GLIBC*.md`: Reference notes comparing algorithms/assembly with glibc’s `libm`.
- `target/`: Build artifacts (ignored by Git).

## Build, Test, and Development Commands

- `cargo build` / `cargo build --release`: Compile the crate (debug/release).
- `cargo test`: Run unit tests (currently located in `src/lib.rs`).
- `cargo bench`: Run Criterion benchmarks (Cargo auto-discovers `./benches/`).
- `cargo fmt`: Auto-format Rust code with `rustfmt`.
- `cargo clippy --all-targets --all-features`: Lint for correctness/performance issues.
- `cargo doc --no-deps`: Generate local API docs.

## Benchmarking & glibc Inspection Tools

- Use **Criterion** to validate wins and catch regressions: compare `fastmaths::*` against `f64::{exp,ln,sin,cos}`.
- For building a reference glibc with AVX/FMA optimizations, use:
    - `./build_libm.sh` (builds libm from source with `-march=native`)
- When debugging parity/perf vs glibc, inspect the system `libm.so.6` or the newly built one with:
    - `nm -D /path/to/libm.so.6 | rg ' (exp|log|sin|cos)$'` (symbol discovery)
  - `objdump -d -C /path/to/libm.so.6` (assembly inspection)
- Keep investigation notes in `MATHS_GLIBC*.md`.

## Coding Style & Naming Conventions

- Formatting: run `cargo fmt` and avoid manual reflow of large constant tables.
- Indentation: 4 spaces; follow standard Rust conventions (`snake_case` fns/modules, `CamelCase` types).
- Low-level math code may use targeted `#[allow(...)]` for readability/perf; keep allowances minimal and scoped.

## Testing Guidelines

- Prefer unit tests close to the code under test (inline `#[cfg(test)]` modules).
- Name tests descriptively (e.g., `exp_matches_std_for_small_inputs`).
- When adding numeric tests, include edge cases: `±0.0`, subnormals, `NaN`, `±∞`, and large-magnitude arguments.
- **MPFR Validation:** Accuracy tests using the `mpfr` feature must be run in a loop of at least 50 iterations to ensure thorough validation across the input distribution.

## MPFR Proptest Loop Verification

- Use the loop runner to validate repeated MPFR proptest runs and fail fast on regressions:
    - `scripts/run_mpfr_proptest_loop.sh` (parallel by default; override with `PARALLEL=...`, default 8)
    - Optional: pass a run count (default 10), e.g. `scripts/run_mpfr_proptest_loop.sh 50`
    - Logs are saved under `proptest-runs/` (override with `LOG_DIR=/path`)

## Commit & Pull Request Guidelines

- This repository has no commits yet; use a consistent convention such as Conventional Commits (`feat: …`, `fix: …`, `docs: …`).
- PRs should include: a short problem statement, the approach, and how it was validated (commands + key outputs).
- For behavior/perf changes in math code, include before/after accuracy notes or benchmark results (if available).
