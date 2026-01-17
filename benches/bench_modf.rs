use criterion::{BenchmarkGroup, Criterion, black_box};
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{configure_criterion, gen_range, glibc_modf};

fn bench_modf_group(
    group: &mut BenchmarkGroup<'_, criterion::measurement::WallTime>,
    inputs: &[f64],
) {
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in inputs {
                let (frac, int) = fastlibm::modf(black_box(x));
                acc += frac + int;
            }
            black_box(acc)
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            for &x in inputs {
                let (frac, int) = glibc_modf(black_box(x));
                acc += frac + int;
            }
            black_box(acc)
        })
    });
}

fn bench_modf(c: &mut Criterion) {
    let inputs = [-10.5, -1.5, -0.0, 0.0, 0.5, 1.5, 10.25, 100.75];
    let common = gen_range(1024, -100.0, 100.0, 0x1b01);

    let mut group = c.benchmark_group("modf/smoke");
    bench_modf_group(&mut group, &inputs);
    group.finish();

    let mut group = c.benchmark_group("modf/common");
    bench_modf_group(&mut group, &common);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_modf(&mut c);
    c.final_summary();
}
