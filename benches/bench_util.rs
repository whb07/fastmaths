#![allow(dead_code)]

use criterion::{BenchmarkGroup, Criterion, black_box};
use std::sync::OnceLock;
use std::time::Duration;

const RNG_A: u64 = 6364136223846793005;
const RNG_C: u64 = 1442695040888963407;
const RNG_DENOM: f64 = (1u64 << 53) as f64;

pub fn lcg_next(state: &mut u64) -> u64 {
    *state = state.wrapping_mul(RNG_A).wrapping_add(RNG_C);
    *state
}

pub fn uniform_f64(state: &mut u64) -> f64 {
    let bits = lcg_next(state) >> 11;
    (bits as f64) / RNG_DENOM
}

pub fn gen_range(count: usize, min: f64, max: f64, seed: u64) -> Vec<f64> {
    let mut state = seed;
    let span = max - min;
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(min + uniform_f64(&mut state) * span);
    }
    values
}

pub fn gen_pairs(count: usize, min: f64, max: f64, seed: u64) -> Vec<(f64, f64)> {
    let mut state = seed;
    let span = max - min;
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        let x = min + uniform_f64(&mut state) * span;
        let y = min + uniform_f64(&mut state) * span;
        values.push((x, y));
    }
    values
}

pub fn bench_inputs<F, G>(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[f64],
    fast: F,
    glibc: G,
) where
    F: Fn(f64) -> f64 + Copy,
    G: Fn(f64) -> f64 + Copy,
{
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in inputs {
                acc += fast(black_box(x));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in inputs {
                acc += glibc(black_box(x));
            }
            black_box(acc)
        })
    });
}

pub fn bench_inputs2<F, G>(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[(f64, f64)],
    fast: F,
    glibc: G,
) where
    F: Fn(f64, f64) -> f64 + Copy,
    G: Fn(f64, f64) -> f64 + Copy,
{
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, y) in inputs {
                acc += fast(black_box(x), black_box(y));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, y) in inputs {
                acc += glibc(black_box(x), black_box(y));
            }
            black_box(acc)
        })
    });
}

pub fn configure_criterion() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(Duration::from_secs(10))
        .warm_up_time(Duration::from_secs(5))
}

struct LibmFns {
    exp: unsafe extern "C" fn(f64) -> f64,
    exp2: unsafe extern "C" fn(f64) -> f64,
    expm1: unsafe extern "C" fn(f64) -> f64,
    log: unsafe extern "C" fn(f64) -> f64,
    log2: unsafe extern "C" fn(f64) -> f64,
    log10: unsafe extern "C" fn(f64) -> f64,
    log1p: unsafe extern "C" fn(f64) -> f64,
    sin: unsafe extern "C" fn(f64) -> f64,
    cos: unsafe extern "C" fn(f64) -> f64,
    tan: unsafe extern "C" fn(f64) -> f64,
    asin: unsafe extern "C" fn(f64) -> f64,
    acos: unsafe extern "C" fn(f64) -> f64,
    atan: unsafe extern "C" fn(f64) -> f64,
    atan2: unsafe extern "C" fn(f64, f64) -> f64,
    sinh: unsafe extern "C" fn(f64) -> f64,
    cosh: unsafe extern "C" fn(f64) -> f64,
    tanh: unsafe extern "C" fn(f64) -> f64,
    hypot: unsafe extern "C" fn(f64, f64) -> f64,
    fmod: unsafe extern "C" fn(f64, f64) -> f64,
    remainder: unsafe extern "C" fn(f64, f64) -> f64,
    pow: unsafe extern "C" fn(f64, f64) -> f64,
    sqrt: unsafe extern "C" fn(f64) -> f64,
    cbrt: unsafe extern "C" fn(f64) -> f64,
}

static LIBM_FNS: OnceLock<LibmFns> = OnceLock::new();

fn libm_path() -> String {
    if let Ok(value) = std::env::var("FASTLIBM_GLIBC_LIBM") {
        let value = value.trim().to_string();
        if !value.is_empty() {
            return value;
        }
    }
    let default = "/tmp/maths/glibc-build/math/libm.so";
    if std::path::Path::new(default).exists() {
        return default.to_string();
    }
    panic!("glibc libm not found; set FASTLIBM_GLIBC_LIBM");
}

fn load_libm() -> LibmFns {
    let path = libm_path();
    let lib = unsafe { libloading::Library::new(&path).expect("load glibc libm") };
    let lib = Box::leak(Box::new(lib));
    unsafe {
        let exp: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"exp").expect("load exp");
        let exp2: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"exp2").expect("load exp2");
        let expm1: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"expm1").expect("load expm1");
        let log: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"log").expect("load log");
        let log2: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"log2").expect("load log2");
        let log10: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"log10").expect("load log10");
        let log1p: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"log1p").expect("load log1p");
        let sin: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"sin").expect("load sin");
        let cos: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"cos").expect("load cos");
        let tan: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"tan").expect("load tan");
        let asin: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"asin").expect("load asin");
        let acos: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"acos").expect("load acos");
        let atan: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"atan").expect("load atan");
        let atan2: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"atan2").expect("load atan2");
        let sinh: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"sinh").expect("load sinh");
        let cosh: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"cosh").expect("load cosh");
        let tanh: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"tanh").expect("load tanh");
        let hypot: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"hypot").expect("load hypot");
        let fmod: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"fmod").expect("load fmod");
        let remainder: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"remainder").expect("load remainder");
        let pow: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"pow").expect("load pow");
        let sqrt: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"sqrt").expect("load sqrt");
        let cbrt: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"cbrt").expect("load cbrt");
        eprintln!("Using libm from {path}");
        LibmFns {
            exp: *exp,
            exp2: *exp2,
            expm1: *expm1,
            log: *log,
            log2: *log2,
            log10: *log10,
            log1p: *log1p,
            sin: *sin,
            cos: *cos,
            tan: *tan,
            asin: *asin,
            acos: *acos,
            atan: *atan,
            atan2: *atan2,
            sinh: *sinh,
            cosh: *cosh,
            tanh: *tanh,
            hypot: *hypot,
            fmod: *fmod,
            remainder: *remainder,
            pow: *pow,
            sqrt: *sqrt,
            cbrt: *cbrt,
        }
    }
}

fn libm() -> &'static LibmFns {
    LIBM_FNS.get_or_init(load_libm)
}

#[inline(never)]
pub fn glibc_exp(x: f64) -> f64 {
    unsafe { (libm().exp)(x) }
}

#[inline(never)]
pub fn glibc_exp2(x: f64) -> f64 {
    unsafe { (libm().exp2)(x) }
}

#[inline(never)]
pub fn glibc_expm1(x: f64) -> f64 {
    unsafe { (libm().expm1)(x) }
}

#[inline(never)]
pub fn glibc_log(x: f64) -> f64 {
    unsafe { (libm().log)(x) }
}

#[inline(never)]
pub fn glibc_log2(x: f64) -> f64 {
    unsafe { (libm().log2)(x) }
}

#[inline(never)]
pub fn glibc_log10(x: f64) -> f64 {
    unsafe { (libm().log10)(x) }
}

#[inline(never)]
pub fn glibc_log1p(x: f64) -> f64 {
    unsafe { (libm().log1p)(x) }
}

#[inline(never)]
pub fn glibc_sin(x: f64) -> f64 {
    unsafe { (libm().sin)(x) }
}

#[inline(never)]
pub fn glibc_cos(x: f64) -> f64 {
    unsafe { (libm().cos)(x) }
}

#[inline(never)]
pub fn glibc_tan(x: f64) -> f64 {
    unsafe { (libm().tan)(x) }
}

#[inline(never)]
pub fn glibc_asin(x: f64) -> f64 {
    unsafe { (libm().asin)(x) }
}

#[inline(never)]
pub fn glibc_acos(x: f64) -> f64 {
    unsafe { (libm().acos)(x) }
}

#[inline(never)]
pub fn glibc_atan(x: f64) -> f64 {
    unsafe { (libm().atan)(x) }
}

#[inline(never)]
pub fn glibc_atan2(y: f64, x: f64) -> f64 {
    unsafe { (libm().atan2)(y, x) }
}

#[inline(never)]
pub fn glibc_sinh(x: f64) -> f64 {
    unsafe { (libm().sinh)(x) }
}

#[inline(never)]
pub fn glibc_cosh(x: f64) -> f64 {
    unsafe { (libm().cosh)(x) }
}

#[inline(never)]
pub fn glibc_tanh(x: f64) -> f64 {
    unsafe { (libm().tanh)(x) }
}

#[inline(never)]
pub fn glibc_hypot(x: f64, y: f64) -> f64 {
    unsafe { (libm().hypot)(x, y) }
}

#[inline(never)]
pub fn glibc_fmod(x: f64, y: f64) -> f64 {
    unsafe { (libm().fmod)(x, y) }
}

#[inline(never)]
pub fn glibc_remainder(x: f64, y: f64) -> f64 {
    unsafe { (libm().remainder)(x, y) }
}

#[inline(never)]
pub fn glibc_pow(x: f64, y: f64) -> f64 {
    unsafe { (libm().pow)(x, y) }
}

#[inline(never)]
pub fn glibc_sqrt(x: f64) -> f64 {
    unsafe { (libm().sqrt)(x) }
}

#[inline(never)]
pub fn glibc_cbrt(x: f64) -> f64 {
    unsafe { (libm().cbrt)(x) }
}
