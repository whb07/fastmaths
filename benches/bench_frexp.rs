use criterion::{Criterion, black_box};
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{configure_criterion, gen_range, glibc_frexp};

fn bench_frexp(c: &mut Criterion) {
    let inputs = [
        0.0, -0.0, 1.0, -1.0, 0.5, 2.0, 1e-300, -1e-300, 1e300, -1e300,
    ];
    let common = gen_range(2048, -1e300, 1e300, 0x2a2b);

    let mut group = c.benchmark_group("frexp/smoke");
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_e = 0i32;
            for &x in &inputs {
                let (m, e) = fastlibm::frexp(black_box(x));
                acc += m;
                acc_e ^= e;
            }
            black_box((acc, acc_e))
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_e = 0i32;
            for &x in &inputs {
                let (m, e) = glibc_frexp(black_box(x));
                acc += m;
                acc_e ^= e;
            }
            black_box((acc, acc_e))
        })
    });
    group.finish();

    let mut group = c.benchmark_group("frexp/common");
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_e = 0i32;
            for &x in &common {
                let (m, e) = fastlibm::frexp(black_box(x));
                acc += m;
                acc_e ^= e;
            }
            black_box((acc, acc_e))
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_e = 0i32;
            for &x in &common {
                let (m, e) = glibc_frexp(black_box(x));
                acc += m;
                acc_e ^= e;
            }
            black_box((acc, acc_e))
        })
    });
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_frexp(&mut c);
    c.final_summary();
}
