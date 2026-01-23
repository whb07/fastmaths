use criterion::Criterion;

mod bench_util;
use bench_util::{bench_inputs_i64_arg, configure_criterion, gen_range, glibc_scalbln};

fn bench_scalbln(c: &mut Criterion) {
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
    let xs = gen_range(2048, -1e300, 1e300, 0x3031);
    let n_values = [-1024_i64, -100, -10, -1, 0, 1, 10, 100, 1024];
    let mut common = Vec::with_capacity(xs.len());
    for (i, x) in xs.iter().enumerate() {
        let n = n_values[i % n_values.len()];
        common.push((*x, n));
    }

    let mut group = c.benchmark_group("scalbln/smoke");
    bench_inputs_i64_arg(&mut group, &inputs, fastmaths::scalbln, glibc_scalbln);
    group.finish();

    let mut group = c.benchmark_group("scalbln/common");
    bench_inputs_i64_arg(&mut group, &common, fastmaths::scalbln, glibc_scalbln);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_scalbln(&mut c);
    c.final_summary();
}
