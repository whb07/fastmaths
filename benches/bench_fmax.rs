use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_fmax};

fn bench_fmax(c: &mut Criterion) {
    let inputs = [(0.0, -0.0), (1.0, 2.0), (-1.0, -2.0), (5.3, 2.1)];
    let pairs = gen_pairs(1024, -100.0, 100.0, 0x1901);

    let mut group = c.benchmark_group("fmax/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::fmax, glibc_fmax);
    group.finish();

    let mut group = c.benchmark_group("fmax/common");
    bench_inputs2(&mut group, &pairs, fastlibm::fmax, glibc_fmax);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_fmax(&mut c);
    c.final_summary();
}
