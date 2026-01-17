use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs3, configure_criterion, gen_triples, glibc_fma};

fn bench_fma(c: &mut Criterion) {
    let inputs = [
        (0.0, 0.0, 0.0),
        (1.0, 1.0, 1.0),
        (-1.0, 2.0, -3.0),
        (1e-6, 1e-6, 1e-6),
        (1e6, -1e6, 1.0),
        (1.5, 2.5, -3.5),
    ];
    let common = gen_triples(2048, -100.0, 100.0, 0x2829);

    let mut group = c.benchmark_group("fma/smoke");
    bench_inputs3(&mut group, &inputs, fastlibm::fma, glibc_fma);
    group.finish();

    let mut group = c.benchmark_group("fma/common");
    bench_inputs3(&mut group, &common, fastlibm::fma, glibc_fma);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_fma(&mut c);
    c.final_summary();
}
