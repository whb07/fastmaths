use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, glibc_exp};

fn bench_exp(c: &mut Criterion) {
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
    let mut group = c.benchmark_group("exp/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::exp, glibc_exp);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_exp(&mut c);
    c.final_summary();
}
