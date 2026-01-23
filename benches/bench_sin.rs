use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_sin};

fn bench_sin(c: &mut Criterion) {
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

    let mut group = c.benchmark_group("sin/smoke");
    bench_inputs(&mut group, &inputs, fastmaths::sin, glibc_sin);
    group.finish();

    let mut group = c.benchmark_group("sin/common");
    bench_inputs(&mut group, &common, fastmaths::sin, glibc_sin);
    group.finish();

    let mut group = c.benchmark_group("sin/medium");
    bench_inputs(&mut group, &medium, fastmaths::sin, glibc_sin);
    group.finish();

    let mut group = c.benchmark_group("sin/huge");
    bench_inputs(&mut group, &huge, fastmaths::sin, glibc_sin);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_sin(&mut c);
    c.final_summary();
}
