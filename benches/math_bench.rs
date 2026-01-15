use criterion::{Criterion, black_box};
use fastmaths::fastlibm;
use std::time::Duration;

use std::sync::OnceLock;

const RNG_A: u64 = 6364136223846793005;
const RNG_C: u64 = 1442695040888963407;
const RNG_DENOM: f64 = (1u64 << 53) as f64;

fn lcg_next(state: &mut u64) -> u64 {
    *state = state.wrapping_mul(RNG_A).wrapping_add(RNG_C);
    *state
}

fn uniform_f64(state: &mut u64) -> f64 {
    let bits = lcg_next(state) >> 11;
    (bits as f64) / RNG_DENOM
}

fn gen_range(count: usize, min: f64, max: f64, seed: u64) -> Vec<f64> {
    let mut state = seed;
    let span = max - min;
    let mut values = Vec::with_capacity(count);
    for _ in 0..count {
        values.push(min + uniform_f64(&mut state) * span);
    }
    values
}

fn bench_inputs<F, G>(c: &mut Criterion, name: &str, inputs: &[f64], fast: F, glibc: G)
where
    F: Fn(f64) -> f64 + Copy,
    G: Fn(f64) -> f64 + Copy,
{
    let mut group = c.benchmark_group(name);
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
    group.finish();
}

struct LibmFns {
    exp: unsafe extern "C" fn(f64) -> f64,
    log: unsafe extern "C" fn(f64) -> f64,
    sin: unsafe extern "C" fn(f64) -> f64,
    cos: unsafe extern "C" fn(f64) -> f64,
}

static LIBM_FNS: OnceLock<LibmFns> = OnceLock::new();

fn libm_path() -> Option<String> {
    if let Ok(value) = std::env::var("FASTLIBM_GLIBC_LIBM") {
        let value = value.trim().to_string();
        if !value.is_empty() {
            return Some(value);
        }
    }
    let default = "/tmp/maths/glibc-build/math/libm.so";
    if std::path::Path::new(default).exists() {
        return Some(default.to_string());
    }
    None
}

fn load_libm() -> LibmFns {
    let path = libm_path().expect("glibc libm not found; set FASTLIBM_GLIBC_LIBM");
    let lib = unsafe { libloading::Library::new(&path).expect("load glibc libm") };
    let lib = Box::leak(Box::new(lib));
    unsafe {
        let exp: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"exp").expect("load exp");
        let log: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"log").expect("load log");
        let sin: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"sin").expect("load sin");
        let cos: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> =
            lib.get(b"cos").expect("load cos");
        eprintln!("Using libm from {path}");
        LibmFns {
            exp: *exp,
            log: *log,
            sin: *sin,
            cos: *cos,
        }
    }
}

fn libm() -> &'static LibmFns {
    LIBM_FNS.get_or_init(load_libm)
}

#[inline(never)]
fn glibc_exp(x: f64) -> f64 {
    let fns = libm();
    unsafe { (fns.exp)(x) }
}

#[inline(never)]
fn glibc_log(x: f64) -> f64 {
    let fns = libm();
    unsafe { (fns.log)(x) }
}

#[inline(never)]
fn glibc_sin(x: f64) -> f64 {
    let fns = libm();
    unsafe { (fns.sin)(x) }
}

#[inline(never)]
fn glibc_cos(x: f64) -> f64 {
    let fns = libm();
    unsafe { (fns.cos)(x) }
}

fn bench_only() -> Option<String> {
    if let Ok(value) = std::env::var("FASTLIBM_BENCH_ONLY") {
        let value = value.trim().to_lowercase();
        if !value.is_empty() {
            return Some(value);
        }
    }

    let mut iter = std::env::args().skip(1).peekable();
    while let Some(arg) = iter.next() {
        if arg == "--only" {
            if let Some(value) = iter.next() {
                let value = value.trim().to_lowercase();
                if !value.is_empty() {
                    return Some(value);
                }
            }
            continue;
        }

        // Ignore Criterion CLI options that take a value so their value doesn't get
        // misinterpreted as a positional filter (e.g. `--sample-size 50`).
        match arg.as_str() {
            "--sample-size"
            | "--warm-up-time"
            | "--measurement-time"
            | "--nresamples"
            | "--noise-threshold"
            | "--significance-level"
            | "--confidence-level"
            | "--plotting-backend"
            | "--output-directory"
            | "--baseline"
            | "--save-baseline"
            | "--load-baseline"
            | "--filter" => {
                let _ = iter.next();
                continue;
            }
            _ => {}
        }

        if !arg.starts_with('-') {
            return Some(arg.trim().to_lowercase());
        }
    }

    None
}

fn bench_enabled(name: &str) -> bool {
    match bench_only() {
        Some(filter) => name.contains(&filter),
        None => true,
    }
}

fn bench_exp(c: &mut Criterion) {
    if !bench_enabled("exp") {
        return;
    }

    let inputs = [
        -745.133_219_101_941_1,
        -100.0,
        -20.0,
        -1.0,
        -1e-6,
        0.0,
        1e-6,
        0.5,
        1.0,
        2.0,
        10.0,
        100.0,
        700.0,
    ];
    bench_inputs(c, "exp/smoke", &inputs, fastlibm::exp, glibc_exp);
}

fn bench_ln(c: &mut Criterion) {
    if !bench_enabled("ln") && !bench_enabled("log") {
        return;
    }

    let inputs = [
        f64::MIN_POSITIVE,
        1e-300,
        1e-20,
        1e-6,
        0.1,
        0.5,
        0.9,
        1.0,
        1.000_000_000_001,
        2.0,
        10.0,
        1e5,
        1e100,
    ];
    bench_inputs(c, "ln/smoke", &inputs, fastlibm::ln, glibc_log);
}

fn bench_sin(c: &mut Criterion) {
    if !bench_enabled("sin") {
        return;
    }

    let inputs = [
        0.0,
        1e-6,
        -1e-6,
        0.5,
        1.0,
        -1.0,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI,
        10.0,
        -10.0,
        1e6,
        -1e6,
    ];
    let common = gen_range(
        1024,
        -2.0 * std::f64::consts::PI,
        2.0 * std::f64::consts::PI,
        0x1357,
    );
    let medium = gen_range(1024, -1e6, 1e6, 0x2468);
    let huge = gen_range(1024, -1e20, 1e20, 0x9abc);

    bench_inputs(c, "sin/smoke", &inputs, fastlibm::sin, glibc_sin);
    bench_inputs(c, "sin/common", &common, fastlibm::sin, glibc_sin);
    bench_inputs(c, "sin/medium", &medium, fastlibm::sin, glibc_sin);
    bench_inputs(c, "sin/huge", &huge, fastlibm::sin, glibc_sin);
}

fn bench_cos(c: &mut Criterion) {
    if !bench_enabled("cos") {
        return;
    }

    let inputs = [
        0.0,
        1e-6,
        -1e-6,
        0.5,
        1.0,
        -1.0,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI,
        10.0,
        -10.0,
        1e6,
        -1e6,
    ];
    let common = gen_range(
        1024,
        -2.0 * std::f64::consts::PI,
        2.0 * std::f64::consts::PI,
        0x1357,
    );
    let medium = gen_range(1024, -1e6, 1e6, 0x2468);
    let huge = gen_range(1024, -1e20, 1e20, 0x9abc);

    bench_inputs(c, "cos/smoke", &inputs, fastlibm::cos, glibc_cos);
    bench_inputs(c, "cos/common", &common, fastlibm::cos, glibc_cos);
    bench_inputs(c, "cos/medium", &medium, fastlibm::cos, glibc_cos);
    bench_inputs(c, "cos/huge", &huge, fastlibm::cos, glibc_cos);
}

fn configure_criterion() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(Duration::from_secs(10))
        .warm_up_time(Duration::from_secs(5))
}

fn main() {
    let mut c = configure_criterion();
    bench_exp(&mut c);
    bench_ln(&mut c);
    bench_sin(&mut c);
    bench_cos(&mut c);
    c.final_summary();
}
