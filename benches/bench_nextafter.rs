use criterion::Criterion;
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_nextafter};

fn bench_nextafter(c: &mut Criterion) {
    let inputs = [
        (0.0, 1.0),
        (0.0, -1.0),
        (1.0, 2.0),
        (1.0, 0.0),
        (-1.0, -2.0),
    ];
    let pairs = gen_pairs(1024, -100.0, 100.0, 0x1701);

    let mut group = c.benchmark_group("nextafter/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::nextafter, glibc_nextafter);
    group.finish();

    let mut group = c.benchmark_group("nextafter/common");
    bench_inputs2(&mut group, &pairs, fastlibm::nextafter, glibc_nextafter);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_nextafter(&mut c);
    c.final_summary();
}
