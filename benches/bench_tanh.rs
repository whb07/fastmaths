use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs, configure_criterion, gen_range, glibc_tanh};

fn bench_tanh(c: &mut Criterion) {
    let inputs = [-20.0, -10.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 10.0, 20.0];
    let common = gen_range(1024, -5.0, 5.0, 0x358a);
    let wide = gen_range(1024, -20.0, 20.0, 0x358b);

    let mut group = c.benchmark_group("tanh/smoke");
    bench_inputs(&mut group, &inputs, fastlibm::tanh, glibc_tanh);
    group.finish();

    let mut group = c.benchmark_group("tanh/common");
    bench_inputs(&mut group, &common, fastlibm::tanh, glibc_tanh);
    group.finish();

    let mut group = c.benchmark_group("tanh/wide");
    bench_inputs(&mut group, &wide, fastlibm::tanh, glibc_tanh);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_tanh(&mut c);
    c.final_summary();
}
