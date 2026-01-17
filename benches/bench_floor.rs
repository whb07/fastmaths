use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_floor};

fn bench_floor(c: &mut Criterion) {
    let inputs = [
        -1e6, -100.5, -10.1, -1.0, -0.9, -0.5, -0.0, 0.0, 0.5, 0.9, 1.0, 1.1, 10.9, 100.5, 1e6,
    ];
    let common = gen_range(2048, -1e6, 1e6, 0x1011);

    let mut group = c.benchmark_group("floor/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::floor, glibc_floor);
    group.finish();

    let mut group = c.benchmark_group("floor/common");
    bench_inputs(&mut group, &common, fastlibm::floor, glibc_floor);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_floor(&mut c);
    c.final_summary();
}
