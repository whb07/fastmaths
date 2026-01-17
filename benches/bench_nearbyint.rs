use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_nearbyint};

fn bench_nearbyint(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.5, -1.5, -1.1, -1.0, -0.5, -0.0, 0.0, 0.5, 1.0, 1.1, 1.5, 10.5, 100.5,
        1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1a1b);

    let mut group = c.benchmark_group("nearbyint/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::nearbyint, glibc_nearbyint);
    group.finish();

    let mut group = c.benchmark_group("nearbyint/common");
    bench_inputs(&mut group, &common, fastlibm::nearbyint, glibc_nearbyint);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_nearbyint(&mut c);
    c.final_summary();
}
