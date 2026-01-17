use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_tgamma};

fn bench_tgamma(c: &mut Criterion) {
    let inputs = [-10.5, -3.5, -1.5, -0.5, 0.5, 1.0, 1.5, 2.0, 3.0, 10.0];
    let mut common = gen_range(2048, -10.0, 10.0, 0x3637);
    for x in &mut common {
        if *x <= 0.0 && *x == x.trunc() {
            *x += 0.5;
        }
    }

    let mut group = c.benchmark_group("tgamma/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::tgamma, glibc_tgamma);
    group.finish();

    let mut group = c.benchmark_group("tgamma/common");
    bench_inputs(&mut group, &common, fastlibm::tgamma, glibc_tgamma);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_tgamma(&mut c);
    c.final_summary();
}
