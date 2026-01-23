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

pub fn gen_triples(count: usize, min: f64, max: f64, seed: u64) -> Vec<(f64, f64, f64)> {
    let mut state = seed;
    let span = max - min;
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        let x = min + uniform_f64(&mut state) * span;
        let y = min + uniform_f64(&mut state) * span;
        let z = min + uniform_f64(&mut state) * span;
        values.push((x, y, z));
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
    group.bench_function("fastmaths", |b| {
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

pub fn bench_inputs_i64<F, G>(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[f64],
    fast: F,
    glibc: G,
) where
    F: Fn(f64) -> i64 + Copy,
    G: Fn(f64) -> i64 + Copy,
{
    group.bench_function("fastmaths", |b| {
        b.iter(|| {
            let mut acc: i64 = 0;
            for &x in inputs {
                acc = acc.wrapping_add(fast(black_box(x)));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc: i64 = 0;
            for &x in inputs {
                acc = acc.wrapping_add(glibc(black_box(x)));
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
    group.bench_function("fastmaths", |b| {
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

pub fn bench_inputs3<F, G>(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[(f64, f64, f64)],
    fast: F,
    glibc: G,
) where
    F: Fn(f64, f64, f64) -> f64 + Copy,
    G: Fn(f64, f64, f64) -> f64 + Copy,
{
    group.bench_function("fastmaths", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, y, z) in inputs {
                acc += fast(black_box(x), black_box(y), black_box(z));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, y, z) in inputs {
                acc += glibc(black_box(x), black_box(y), black_box(z));
            }
            black_box(acc)
        })
    });
}

pub fn bench_inputs_i32_arg<F, G>(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[(f64, i32)],
    fast: F,
    glibc: G,
) where
    F: Fn(f64, i32) -> f64 + Copy,
    G: Fn(f64, i32) -> f64 + Copy,
{
    group.bench_function("fastmaths", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, n) in inputs {
                acc += fast(black_box(x), black_box(n));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, n) in inputs {
                acc += glibc(black_box(x), black_box(n));
            }
            black_box(acc)
        })
    });
}

pub fn bench_inputs_i64_arg<F, G>(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[(f64, i64)],
    fast: F,
    glibc: G,
) where
    F: Fn(f64, i64) -> f64 + Copy,
    G: Fn(f64, i64) -> f64 + Copy,
{
    group.bench_function("fastmaths", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, n) in inputs {
                acc += fast(black_box(x), black_box(n));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &(x, n) in inputs {
                acc += glibc(black_box(x), black_box(n));
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
    asinh: unsafe extern "C" fn(f64) -> f64,
    acosh: unsafe extern "C" fn(f64) -> f64,
    atanh: unsafe extern "C" fn(f64) -> f64,
    erf: unsafe extern "C" fn(f64) -> f64,
    erfc: unsafe extern "C" fn(f64) -> f64,
    exp10: unsafe extern "C" fn(f64) -> f64,
    atan: unsafe extern "C" fn(f64) -> f64,
    atan2: unsafe extern "C" fn(f64, f64) -> f64,
    sinh: unsafe extern "C" fn(f64) -> f64,
    cosh: unsafe extern "C" fn(f64) -> f64,
    tanh: unsafe extern "C" fn(f64) -> f64,
    hypot: unsafe extern "C" fn(f64, f64) -> f64,
    fmod: unsafe extern "C" fn(f64, f64) -> f64,
    fdim: unsafe extern "C" fn(f64, f64) -> f64,
    fmax: unsafe extern "C" fn(f64, f64) -> f64,
    fmin: unsafe extern "C" fn(f64, f64) -> f64,
    modf: unsafe extern "C" fn(f64, *mut f64) -> f64,
    logb: unsafe extern "C" fn(f64) -> f64,
    ilogb: unsafe extern "C" fn(f64) -> i32,
    nextafter: unsafe extern "C" fn(f64, f64) -> f64,
    remainder: unsafe extern "C" fn(f64, f64) -> f64,
    pow: unsafe extern "C" fn(f64, f64) -> f64,
    sqrt: unsafe extern "C" fn(f64) -> f64,
    cbrt: unsafe extern "C" fn(f64) -> f64,
    floor: unsafe extern "C" fn(f64) -> f64,
    ceil: unsafe extern "C" fn(f64) -> f64,
    trunc: unsafe extern "C" fn(f64) -> f64,
    round: unsafe extern "C" fn(f64) -> f64,
    rint: unsafe extern "C" fn(f64) -> f64,
    nearbyint: unsafe extern "C" fn(f64) -> f64,
    lrint: unsafe extern "C" fn(f64) -> i64,
    llrint: unsafe extern "C" fn(f64) -> i64,
    lround: unsafe extern "C" fn(f64) -> i64,
    llround: unsafe extern "C" fn(f64) -> i64,
    copysign: unsafe extern "C" fn(f64, f64) -> f64,
    fabs: unsafe extern "C" fn(f64) -> f64,
    frexp: unsafe extern "C" fn(f64, *mut i32) -> f64,
    ldexp: unsafe extern "C" fn(f64, i32) -> f64,
    scalbn: unsafe extern "C" fn(f64, i32) -> f64,
    scalbln: unsafe extern "C" fn(f64, i64) -> f64,
    fma: unsafe extern "C" fn(f64, f64, f64) -> f64,
    remquo: unsafe extern "C" fn(f64, f64, *mut i32) -> f64,
    lgamma: unsafe extern "C" fn(f64) -> f64,
    tgamma: unsafe extern "C" fn(f64) -> f64,
}

static LIBM_FNS: OnceLock<LibmFns> = OnceLock::new();

fn libm_path() -> String {
    if let Ok(value) = std::env::var("FASTMATHS_GLIBC_LIBM") {
        let value = value.trim().to_string();
        if !value.is_empty() {
            return value;
        }
    }
    let default = "/tmp/maths/glibc-build/math/libm.so";
    if std::path::Path::new(default).exists() {
        return default.to_string();
    }
    panic!("glibc libm not found; set FASTMATHS_GLIBC_LIBM");
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
        let asinh: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"asinh").expect("load asinh");
        let acosh: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"acosh").expect("load acosh");
        let atanh: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"atanh").expect("load atanh");
        let erf: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"erf").expect("load erf");
        let erfc: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"erfc").expect("load erfc");
        let exp10: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"exp10").expect("load exp10");
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
        let fdim: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"fdim").expect("load fdim");
        let fmax: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"fmax").expect("load fmax");
        let fmin: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"fmin").expect("load fmin");
        let modf: libloading::Symbol<unsafe extern "C" fn(f64, *mut f64) -> f64> =
            lib.get(b"modf").expect("load modf");
        let logb: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"logb").expect("load logb");
        let ilogb: libloading::Symbol<unsafe extern "C" fn(f64) -> i32> =
            lib.get(b"ilogb").expect("load ilogb");
        let nextafter: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"nextafter").expect("load nextafter");
        let remainder: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"remainder").expect("load remainder");
        let pow: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"pow").expect("load pow");
        let sqrt: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"sqrt").expect("load sqrt");
        let cbrt: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"cbrt").expect("load cbrt");
        let floor: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"floor").expect("load floor");
        let ceil: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"ceil").expect("load ceil");
        let trunc: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"trunc").expect("load trunc");
        let round: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"round").expect("load round");
        let rint: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"rint").expect("load rint");
        let nearbyint: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"nearbyint").expect("load nearbyint");
        let lrint: libloading::Symbol<unsafe extern "C" fn(f64) -> i64> =
            lib.get(b"lrint").expect("load lrint");
        let llrint: libloading::Symbol<unsafe extern "C" fn(f64) -> i64> =
            lib.get(b"llrint").expect("load llrint");
        let lround: libloading::Symbol<unsafe extern "C" fn(f64) -> i64> =
            lib.get(b"lround").expect("load lround");
        let llround: libloading::Symbol<unsafe extern "C" fn(f64) -> i64> =
            lib.get(b"llround").expect("load llround");
        let copysign: libloading::Symbol<unsafe extern "C" fn(f64, f64) -> f64> =
            lib.get(b"copysign").expect("load copysign");
        let fabs: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"fabs").expect("load fabs");
        let frexp: libloading::Symbol<unsafe extern "C" fn(f64, *mut i32) -> f64> =
            lib.get(b"frexp").expect("load frexp");
        let ldexp: libloading::Symbol<unsafe extern "C" fn(f64, i32) -> f64> =
            lib.get(b"ldexp").expect("load ldexp");
        let scalbn: libloading::Symbol<unsafe extern "C" fn(f64, i32) -> f64> =
            lib.get(b"scalbn").expect("load scalbn");
        let scalbln: libloading::Symbol<unsafe extern "C" fn(f64, i64) -> f64> =
            lib.get(b"scalbln").expect("load scalbln");
        let fma: libloading::Symbol<unsafe extern "C" fn(f64, f64, f64) -> f64> =
            lib.get(b"fma").expect("load fma");
        let remquo: libloading::Symbol<unsafe extern "C" fn(f64, f64, *mut i32) -> f64> =
            lib.get(b"remquo").expect("load remquo");
        let lgamma: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"lgamma").expect("load lgamma");
        let tgamma: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"tgamma").expect("load tgamma");
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
            asinh: *asinh,
            acosh: *acosh,
            atanh: *atanh,
            erf: *erf,
            erfc: *erfc,
            exp10: *exp10,
            atan: *atan,
            atan2: *atan2,
            sinh: *sinh,
            cosh: *cosh,
            tanh: *tanh,
            hypot: *hypot,
            fmod: *fmod,
            fdim: *fdim,
            fmax: *fmax,
            fmin: *fmin,
            modf: *modf,
            logb: *logb,
            ilogb: *ilogb,
            nextafter: *nextafter,
            remainder: *remainder,
            pow: *pow,
            sqrt: *sqrt,
            cbrt: *cbrt,
            floor: *floor,
            ceil: *ceil,
            trunc: *trunc,
            round: *round,
            rint: *rint,
            nearbyint: *nearbyint,
            lrint: *lrint,
            llrint: *llrint,
            lround: *lround,
            llround: *llround,
            copysign: *copysign,
            fabs: *fabs,
            frexp: *frexp,
            ldexp: *ldexp,
            scalbn: *scalbn,
            scalbln: *scalbln,
            fma: *fma,
            remquo: *remquo,
            lgamma: *lgamma,
            tgamma: *tgamma,
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
pub fn glibc_asinh(x: f64) -> f64 {
    unsafe { (libm().asinh)(x) }
}

#[inline(never)]
pub fn glibc_acosh(x: f64) -> f64 {
    unsafe { (libm().acosh)(x) }
}

#[inline(never)]
pub fn glibc_atanh(x: f64) -> f64 {
    unsafe { (libm().atanh)(x) }
}

#[inline(never)]
pub fn glibc_erf(x: f64) -> f64 {
    unsafe { (libm().erf)(x) }
}

#[inline(never)]
pub fn glibc_erfc(x: f64) -> f64 {
    unsafe { (libm().erfc)(x) }
}

#[inline(never)]
pub fn glibc_exp10(x: f64) -> f64 {
    unsafe { (libm().exp10)(x) }
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
pub fn glibc_fdim(x: f64, y: f64) -> f64 {
    unsafe { (libm().fdim)(x, y) }
}

#[inline(never)]
pub fn glibc_fmax(x: f64, y: f64) -> f64 {
    unsafe { (libm().fmax)(x, y) }
}

#[inline(never)]
pub fn glibc_fmin(x: f64, y: f64) -> f64 {
    unsafe { (libm().fmin)(x, y) }
}

#[inline(never)]
pub fn glibc_modf(x: f64) -> (f64, f64) {
    let mut ip = 0.0;
    let frac = unsafe { (libm().modf)(x, &mut ip as *mut f64) };
    (frac, ip)
}

#[inline(never)]
pub fn glibc_logb(x: f64) -> f64 {
    unsafe { (libm().logb)(x) }
}

#[inline(never)]
pub fn glibc_ilogb(x: f64) -> i32 {
    unsafe { (libm().ilogb)(x) }
}

#[inline(never)]
pub fn glibc_nextafter(x: f64, y: f64) -> f64 {
    unsafe { (libm().nextafter)(x, y) }
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

#[inline(never)]
pub fn glibc_floor(x: f64) -> f64 {
    unsafe { (libm().floor)(x) }
}

#[inline(never)]
pub fn glibc_ceil(x: f64) -> f64 {
    unsafe { (libm().ceil)(x) }
}

#[inline(never)]
pub fn glibc_trunc(x: f64) -> f64 {
    unsafe { (libm().trunc)(x) }
}

#[inline(never)]
pub fn glibc_round(x: f64) -> f64 {
    unsafe { (libm().round)(x) }
}

#[inline(never)]
pub fn glibc_rint(x: f64) -> f64 {
    unsafe { (libm().rint)(x) }
}

#[inline(never)]
pub fn glibc_nearbyint(x: f64) -> f64 {
    unsafe { (libm().nearbyint)(x) }
}

#[inline(never)]
pub fn glibc_lrint(x: f64) -> i64 {
    unsafe { (libm().lrint)(x) }
}

#[inline(never)]
pub fn glibc_llrint(x: f64) -> i64 {
    unsafe { (libm().llrint)(x) }
}

#[inline(never)]
pub fn glibc_lround(x: f64) -> i64 {
    unsafe { (libm().lround)(x) }
}

#[inline(never)]
pub fn glibc_llround(x: f64) -> i64 {
    unsafe { (libm().llround)(x) }
}

#[inline(never)]
pub fn glibc_copysign(x: f64, y: f64) -> f64 {
    unsafe { (libm().copysign)(x, y) }
}

#[inline(never)]
pub fn glibc_fabs(x: f64) -> f64 {
    unsafe { (libm().fabs)(x) }
}

#[inline(never)]
pub fn glibc_frexp(x: f64) -> (f64, i32) {
    let mut exp = 0;
    let mant = unsafe { (libm().frexp)(x, &mut exp as *mut i32) };
    (mant, exp)
}

#[inline(never)]
pub fn glibc_ldexp(x: f64, n: i32) -> f64 {
    unsafe { (libm().ldexp)(x, n) }
}

#[inline(never)]
pub fn glibc_scalbn(x: f64, n: i32) -> f64 {
    unsafe { (libm().scalbn)(x, n) }
}

#[inline(never)]
pub fn glibc_scalbln(x: f64, n: i64) -> f64 {
    unsafe { (libm().scalbln)(x, n) }
}

#[inline(never)]
pub fn glibc_fma(x: f64, y: f64, z: f64) -> f64 {
    unsafe { (libm().fma)(x, y, z) }
}

#[inline(never)]
pub fn glibc_remquo(x: f64, y: f64) -> (f64, i32) {
    let mut q = 0;
    let r = unsafe { (libm().remquo)(x, y, &mut q as *mut i32) };
    (r, q)
}

#[inline(never)]
pub fn glibc_lgamma(x: f64) -> f64 {
    unsafe { (libm().lgamma)(x) }
}

#[inline(never)]
pub fn glibc_tgamma(x: f64) -> f64 {
    unsafe { (libm().tgamma)(x) }
}
