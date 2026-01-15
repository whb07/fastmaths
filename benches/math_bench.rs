use criterion::{Criterion, black_box};
use fastmaths::fastlibm;
use std::time::Duration;

use std::sync::OnceLock;

struct LibmFns {
    exp: unsafe extern "C" fn(f64) -> f64,
    log: unsafe extern "C" fn(f64) -> f64,
    sin: unsafe extern "C" fn(f64) -> f64,
    cos: unsafe extern "C" fn(f64) -> f64,
}

static LIBM_FNS: OnceLock<Option<LibmFns>> = OnceLock::new();

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

fn load_libm() -> Option<LibmFns> {
    let path = libm_path()?;
    let lib = unsafe { libloading::Library::new(&path).ok()? };
    let lib = Box::leak(Box::new(lib));
    unsafe {
        let exp: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> = lib.get(b"exp").ok()?;
        let log: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> = lib.get(b"log").ok()?;
        let sin: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> = lib.get(b"sin").ok()?;
        let cos: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> = lib.get(b"cos").ok()?;
        eprintln!("Using libm from {path}");
        Some(LibmFns {
            exp: *exp,
            log: *log,
            sin: *sin,
            cos: *cos,
        })
    }
}

fn libm() -> Option<&'static LibmFns> {
    LIBM_FNS.get_or_init(load_libm).as_ref()
}

fn glibc_exp(x: f64) -> f64 {
    if let Some(fns) = libm() {
        unsafe { (fns.exp)(x) }
    } else {
        x.exp()
    }
}

fn glibc_log(x: f64) -> f64 {
    if let Some(fns) = libm() {
        unsafe { (fns.log)(x) }
    } else {
        x.ln()
    }
}

fn glibc_sin(x: f64) -> f64 {
    if let Some(fns) = libm() {
        unsafe { (fns.sin)(x) }
    } else {
        x.sin()
    }
}

fn glibc_cos(x: f64) -> f64 {
    if let Some(fns) = libm() {
        unsafe { (fns.cos)(x) }
    } else {
        x.cos()
    }
}

fn bench_only() -> Option<String> {
    if let Ok(value) = std::env::var("FASTLIBM_BENCH_ONLY") {
        let value = value.trim().to_lowercase();
        if !value.is_empty() {
            return Some(value);
        }
    }

    let mut iter = std::env::args().skip(1);
    while let Some(arg) = iter.next() {
        if arg == "--only" {
            if let Some(value) = iter.next() {
                let value = value.trim().to_lowercase();
                if !value.is_empty() {
                    return Some(value);
                }
            }
        } else if !arg.starts_with('-') {
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

    let mut group = c.benchmark_group("exp");
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

    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += fastlibm::exp(black_box(x));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += glibc_exp(black_box(x));
            }
            black_box(acc)
        })
    });
    group.finish();
}

fn bench_ln(c: &mut Criterion) {
    if !bench_enabled("ln") && !bench_enabled("log") {
        return;
    }

    let mut group = c.benchmark_group("ln");
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

    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += fastlibm::ln(black_box(x));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += glibc_log(black_box(x));
            }
            black_box(acc)
        })
    });
    group.finish();
}

fn bench_sin(c: &mut Criterion) {
    if !bench_enabled("sin") {
        return;
    }

    let mut group = c.benchmark_group("sin");
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

    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += fastlibm::sin(black_box(x));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += glibc_sin(black_box(x));
            }
            black_box(acc)
        })
    });
    group.finish();
}

fn bench_cos(c: &mut Criterion) {
    if !bench_enabled("cos") {
        return;
    }

    let mut group = c.benchmark_group("cos");
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

    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += fastlibm::cos(black_box(x));
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in &inputs {
                acc += glibc_cos(black_box(x));
            }
            black_box(acc)
        })
    });
    group.finish();
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
