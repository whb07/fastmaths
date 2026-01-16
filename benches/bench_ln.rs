use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, glibc_log};

fn bench_ln(c: &mut Criterion) {
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
    let mut group = c.benchmark_group("ln/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::ln, glibc_log);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_ln(&mut c);
    c.final_summary();
}
