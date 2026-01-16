use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_remainder};

fn bench_remainder(c: &mut Criterion) {
    let inputs = [(5.3, 2.1), (-5.3, 2.1), (1.0, 0.5), (1.5, 0.5), (1e20, 3.0)];
    let mut common = gen_pairs(1024, -100.0, 100.0, 0x401b);
    for (_, y) in &mut common {
        if y.abs() < 1e-6 {
            *y = if *y >= 0.0 { 1e-6 } else { -1e-6 };
        }
    }

    let mut group = c.benchmark_group("remainder/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::remainder, glibc_remainder);
    group.finish();

    let mut group = c.benchmark_group("remainder/common");
    bench_inputs2(&mut group, &common, fastlibm::remainder, glibc_remainder);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_remainder(&mut c);
    c.final_summary();
}
