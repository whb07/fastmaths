use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_fdim};

fn bench_fdim(c: &mut Criterion) {
    let inputs = [(0.0, 0.0), (1.0, 2.0), (2.0, 1.0), (-1.0, -2.0), (5.3, 2.1)];
    let pairs = gen_pairs(1024, -100.0, 100.0, 0x1801);

    let mut group = c.benchmark_group("fdim/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::fdim, glibc_fdim);
    group.finish();

    let mut group = c.benchmark_group("fdim/common");
    bench_inputs2(&mut group, &pairs, fastlibm::fdim, glibc_fdim);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_fdim(&mut c);
    c.final_summary();
}
