# fastmaths

A high-performance, `no_std` Rust math library aiming for **glibc parity**, **speed**, and **strict accuracy** (≤ 1.0 ULP).

## Features

- **Glibc Parity:** Implements core `libm` functions with behavior matching the GNU C Library (NaNs, infinity, subnormals, and edge cases).
- **High Accuracy:** Functions are verified against **MPFR** (Multiple Precision Floating-Point Reliably) to ensure an error of ≤ 1.0 ULP (Units in the Last Place).
- **Optimized for Speed:** Uses table-driven algorithms and hardware intrinsics (e.g., FMA, SSE2) where available to meet or exceed the performance of the system `libm`.
- **`no_std` Support:** Designed for embedded, kernel, or other environments where the standard library is not available.
- **Zero External Dependencies:** Self-contained implementation (does not link to the system `libm`).

## Implemented Functions

### Elementary & Power

- **Exponential:** `exp`, `exp2`, `expm1`, `exp10`
- **Logarithmic:** `ln`/`log`, `log2`, `log10`, `log1p`
- **Trigonometric:** `sin`, `cos`, `tan`, `atan`, `atan2`, `sincos`
- **Power/Root:** `pow`, `sqrt`, `cbrt`, `hypot`

### Hyperbolic & Inverse Hyperbolic

- **Hyperbolic:** `sinh`, `cosh`, `tanh`
- **Inverse Hyperbolic:** `asinh`, `acosh`, `atanh`

### Special Functions

- **Gamma family:** `lgamma`, `tgamma`
- **Error functions:** `erf`, `erfc`

### IEEE-754 Helpers & Bit-Level Utilities

- **Classification:** `fpclassify`, `isfinite`, `isinf`, `isnan`, `signbit`
- **Rounding:** `rint`, `nearbyint`, `round`, `trunc`, `floor`, `ceil`, `lrint`, `llrint`, `lround`, `llround`
- **Scaling:** `frexp`, `ldexp`, `scalbn`, `scalbln`
- **Min/Max/Delta:** `fmin`, `fmax`, `fdim`
- **Remainders:** `fmod`, `remainder`, `remquo`
- **Adjacency:** `nextafter`
- **Exponent access:** `logb`, `ilogb`
- **FMA:** `fma`
- **Decomposition:** `modf`

## Accuracy Standards

Accuracy is the primary goal of this project. Every function is tested using:
1. **Deterministic Tests:** Checking specific edge cases and known difficult values.
2. **Property-based Testing:** Using `proptest` to verify thousands of random inputs across the entire floating-point range.
3. **MPFR Verification:** Results are compared against the "ground truth" provided by the MPFR library to guarantee ≤ 1.0 ULP accuracy.

Recent precision improvements:

- **`acosh` near 1.0:** uses a log1p-based path with compensated summation for `x < 1.1192` to keep errors within ≤ 1.0 ULP.

## Performance

`fastmaths` is designed to be faster than the system `glibc` implementation. Benchmark results are machine- and toolchain-dependent, so run them locally for current numbers.

To run the full suite of benchmarks:
```bash
cargo bench
```

*Note: For maximum performance, compile with `RUSTFLAGS="-C target-cpu=native"` to enable hardware-specific optimizations like FMA and SSE2 instructions. On x86/x86_64, fastmaths **assumes FMA at compile-time** unless the `soft-fma` feature is enabled (use it for deterministic builds or legacy CPUs). For glibc comparisons, point the benches at a locally-built `libm` via `FASTMATHS_GLIBC_LIBM` or use `./build_libm.sh`.*

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
fastmaths = "0.1.0"
```

Example usage:

```rust
use fastmaths::exp;

let x = 2.0_f64;
let result = exp(x);
println!("e^{} = {}", x, result);
```

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or [MIT license](LICENSE-MIT) at your option.
