use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_lgamma};

fn bench_lgamma(c: &mut Criterion) {
    let inputs = [
        -10.5, -3.5, -1.5, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 10.0, 50.0, 100.0,
    ];
    let mut common = gen_range(2048, -20.0, 20.0, 0x3435);
    for x in &mut common {
        if *x <= 0.0 && *x == x.trunc() {
            *x += 0.5;
        }
    }

    let mut group = c.benchmark_group("lgamma/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::lgamma, glibc_lgamma);
    group.finish();

    let mut group = c.benchmark_group("lgamma/common");
    bench_inputs(&mut group, &common, fastlibm::lgamma, glibc_lgamma);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_lgamma(&mut c);
    c.final_summary();
}
