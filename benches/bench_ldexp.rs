use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs_i32_arg, configure_criterion, gen_range, glibc_ldexp};

fn bench_ldexp(c: &mut Criterion) {
    let inputs = [
        (0.0, 0),
        (1.0, 1),
        (1.0, -1),
        (-1.0, 5),
        (1e-300, 10),
        (-1e-300, -10),
        (1e300, -10),
        (-1e300, 10),
    ];
    let xs = gen_range(2048, -1e300, 1e300, 0x2c2d);
    let n_values = [-1024, -100, -10, -1, 0, 1, 10, 100, 1024];
    let mut common = Vec::with_capacity(xs.len());
    for (i, x) in xs.iter().enumerate() {
        let n = n_values[i % n_values.len()];
        common.push((*x, n));
    }

    let mut group = c.benchmark_group("ldexp/smoke");
    bench_inputs_i32_arg(&mut group, &inputs, fastlibm::ldexp, glibc_ldexp);
    group.finish();

    let mut group = c.benchmark_group("ldexp/common");
    bench_inputs_i32_arg(&mut group, &common, fastlibm::ldexp, glibc_ldexp);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_ldexp(&mut c);
    c.final_summary();
}
