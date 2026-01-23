use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_tan};

fn bench_tan(c: &mut Criterion) {
    let inputs = [
        0.0,
        1e-6,
        -1e-6,
        0.5,
        1.0,
        -1.0,
        std::f64::consts::FRAC_PI_4,
        std::f64::consts::FRAC_PI_2 - 1e-12,
        10.0,
        -10.0,
    ];
    let common = gen_range(
        1024,
        -2.0 * std::f64::consts::PI,
        2.0 * std::f64::consts::PI,
        0x1357,
    );
    let medium = gen_range(1024, -1e6, 1e6, 0x2468);
    let huge = gen_range(1024, -1e20, 1e20, 0x9abc);

    let mut group = c.benchmark_group("tan/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::tan, glibc_tan);
    group.finish();

    let mut group = c.benchmark_group("tan/common");
    bench_inputs(&mut group, &common, fastmaths::tan, glibc_tan);
    group.finish();

    let mut group = c.benchmark_group("tan/medium");
    bench_inputs(&mut group, &medium, fastmaths::tan, glibc_tan);
    group.finish();

    let mut group = c.benchmark_group("tan/huge");
    bench_inputs(&mut group, &huge, fastmaths::tan, glibc_tan);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_tan(&mut c);
    c.final_summary();
}
