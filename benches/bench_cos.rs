use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_cos};

fn bench_cos(c: &mut Criterion) {
    let inputs = [
        0.0,
        1e-6,
        -1e-6,
        0.5,
        1.0,
        -1.0,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI,
        10.0,
        -10.0,
        1e6,
        -1e6,
    ];
    let common = gen_range(
        1024,
        -2.0 * std::f64::consts::PI,
        2.0 * std::f64::consts::PI,
        0x1357,
    );
    let medium = gen_range(1024, -1e6, 1e6, 0x2468);
    let huge = gen_range(1024, -1e20, 1e20, 0x9abc);

    let mut group = c.benchmark_group("cos/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::cos, glibc_cos);
    group.finish();

    let mut group = c.benchmark_group("cos/common");
    bench_inputs(&mut group, &common, fastlibm::cos, glibc_cos);
    group.finish();

    let mut group = c.benchmark_group("cos/medium");
    bench_inputs(&mut group, &medium, fastlibm::cos, glibc_cos);
    group.finish();

    let mut group = c.benchmark_group("cos/huge");
    bench_inputs(&mut group, &huge, fastlibm::cos, glibc_cos);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_cos(&mut c);
    c.final_summary();
}
