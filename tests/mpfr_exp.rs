#![cfg(feature = "mpfr")]

use fastmaths::fastlibm;
use rug::Float;
use std::env;

const MPFR_PREC: u32 = 256;

fn mpfr_exp_f64(x: f64) -> f64 {
    let mut v = Float::with_val(MPFR_PREC, x);
    v.exp_mut();
    v.to_f64()
}

fn ulp_size(x: f64) -> f64 {
    if x == 0.0 {
        return f64::from_bits(1);
    }
    if x.is_nan() || x.is_infinite() {
        return f64::NAN;
    }
    let next = if x.is_sign_negative() {
        x.next_down()
    } else {
        x.next_up()
    };
    (next - x).abs()
}

fn ulp_error(actual: f64, expected: f64) -> f64 {
    let diff = (actual - expected).abs();
    if diff == 0.0 {
        return 0.0;
    }
    let ulp = ulp_size(expected);
    if !ulp.is_finite() || ulp == 0.0 {
        return f64::INFINITY;
    }
    diff / ulp
}

struct LibmFns {
    exp: unsafe extern "C" fn(f64) -> f64,
}

fn glibc_exp_opt() -> Option<LibmFns> {
    let path = env::var("FASTLIBM_GLIBC_LIBM")
        .ok()
        .filter(|v| !v.trim().is_empty())
        .or_else(|| {
            let default = "/tmp/maths/glibc-build/math/libm.so";
            if std::path::Path::new(default).exists() {
                Some(default.to_string())
            } else {
                None
            }
        })?;

    let lib = unsafe { libloading::Library::new(&path).ok()? };
    let lib = Box::leak(Box::new(lib));
    unsafe {
        let exp: libloading::Symbol<unsafe extern "C" fn(f64) -> f64> = lib.get(b"exp").ok()?;
        Some(LibmFns { exp: *exp })
    }
}

fn sweep_offsets(radius: i64, stride: i64) -> Vec<i64> {
    let mut offsets = Vec::new();
    let mut off = -radius;
    while off <= radius {
        offsets.push(off);
        off = off.saturating_add(stride);
        if off == i64::MAX {
            break;
        }
    }
    offsets
}

struct ReportRow {
    label: &'static str,
    x: f64,
    mpfr: f64,
    fast: f64,
    fast_ulps: f64,
    glibc: Option<(f64, f64)>,
}

fn push_report(rows: &mut Vec<ReportRow>, label: &'static str, x: f64, glibc: Option<&LibmFns>) {
    let mpfr = mpfr_exp_f64(x);
    let fast = fastlibm::exp(x);
    let fast_ulps = ulp_error(fast, mpfr);
    let glibc_row = glibc
        .map(|g| unsafe { (g.exp)(x) })
        .map(|v| (v, ulp_error(v, mpfr)));
    rows.push(ReportRow {
        label,
        x,
        mpfr,
        fast,
        fast_ulps,
        glibc: glibc_row,
    });
}

fn print_report(rows: &[ReportRow]) {
    println!("| Case | x | mpfr bits | fast bits | fast ulp | glibc bits | glibc ulp |");
    println!("| :--- | ---: | :--- | :--- | ---: | :--- | ---: |");
    for row in rows {
        let mpfr_bits = format!("{:016x}", row.mpfr.to_bits());
        let fast_bits = format!("{:016x}", row.fast.to_bits());
        let (glibc_bits, glibc_ulps) = match row.glibc {
            Some((v, ulps)) => (format!("{:016x}", v.to_bits()), format!("{ulps:.3}")),
            None => ("n/a".to_string(), "n/a".to_string()),
        };
        println!(
            "| {} | {:.17e} | {} | {} | {:.3} | {} | {} |",
            row.label, row.x, mpfr_bits, fast_bits, row.fast_ulps, glibc_bits, glibc_ulps
        );
    }
}

#[test]
fn mpfr_exp_sweep() {
    let x0 = match env::var("FASTLIBM_MPFR_X") {
        Ok(v) => v.parse::<f64>().expect("FASTLIBM_MPFR_X must be f64"),
        Err(_) => return,
    };
    let radius = env::var("FASTLIBM_MPFR_RADIUS")
        .ok()
        .and_then(|v| v.parse::<i64>().ok())
        .unwrap_or(10_000);
    let stride = env::var("FASTLIBM_MPFR_STRIDE")
        .ok()
        .and_then(|v| v.parse::<i64>().ok())
        .unwrap_or(1);

    let glibc = glibc_exp_opt();
    let base_bits = x0.to_bits();
    let mut max_ulps = 0.0f64;
    let mut max_x = x0;
    let mut first_mismatch: Option<(f64, f64, f64)> = None;
    let mut max_glibc_ulps = 0.0f64;
    let mut max_glibc_x = x0;
    let mut report = Vec::new();

    push_report(&mut report, "x0", x0, glibc.as_ref());

    for offset in sweep_offsets(radius, stride.max(1)) {
        let bits = if offset < 0 {
            base_bits.wrapping_sub((-offset) as u64)
        } else {
            base_bits.wrapping_add(offset as u64)
        };
        let x = f64::from_bits(bits);
        let expected = mpfr_exp_f64(x);
        let actual = fastlibm::exp(x);
        let ulps = ulp_error(actual, expected);
        if ulps > max_ulps {
            max_ulps = ulps;
            max_x = x;
        }
        if first_mismatch.is_none() && ulps != 0.0 {
            first_mismatch = Some((x, actual, expected));
        }

        if let Some(ref glibc) = glibc {
            let g = unsafe { (glibc.exp)(x) };
            let gulps = ulp_error(g, expected);
            if gulps > max_glibc_ulps {
                max_glibc_ulps = gulps;
                max_glibc_x = x;
            }
        }
    }

    println!("MPFR sweep around x0={x0} (radius={radius} stride={stride})");
    println!("fastlibm max ulp error vs MPFR: ulps={max_ulps} at x={max_x}");
    if let Some((x, actual, expected)) = first_mismatch {
        println!(
            "first fastlibm mismatch: x={x} actual={actual:.17e} expected={expected:.17e} ulps={}",
            ulp_error(actual, expected)
        );
    } else {
        println!("no mismatches against MPFR in sweep range");
    }

    if glibc.is_some() {
        println!("glibc max ulp error vs MPFR: ulps={max_glibc_ulps} at x={max_glibc_x}");
    }

    push_report(&mut report, "fastlibm_max", max_x, glibc.as_ref());
    if glibc.is_some() {
        push_report(&mut report, "glibc_max", max_glibc_x, glibc.as_ref());
    }
    if let Some((x, _, _)) = first_mismatch {
        push_report(&mut report, "fastlibm_first", x, glibc.as_ref());
    }

    let report_enabled = env::var("FASTLIBM_MPFR_REPORT")
        .ok()
        .map(|v| v != "0")
        .unwrap_or(true);
    if report_enabled {
        print_report(&report);
    }
}
