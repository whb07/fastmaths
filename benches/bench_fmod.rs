use criterion::Criterion;
use fastmaths as fastlibm;

mod bench_util;
use bench_util::{bench_inputs2, configure_criterion, gen_pairs, glibc_fmod};

fn bench_fmod(c: &mut Criterion) {
    let inputs = [(5.3, 2.1), (-5.3, 2.1), (1.0, 0.5), (1.5, 0.5), (1e20, 3.0)];
    let mut common = gen_pairs(1024, -100.0, 100.0, 0x401a);
    for (_, y) in &mut common {
        if y.abs() < 1e-6 {
            *y = if *y >= 0.0 { 1e-6 } else { -1e-6 };
        }
    }

    let mut group = c.benchmark_group("fmod/smoke");
    bench_inputs2(&mut group, &inputs, fastlibm::fmod, glibc_fmod);
    group.finish();

    let mut group = c.benchmark_group("fmod/common");
    bench_inputs2(&mut group, &common, fastlibm::fmod, glibc_fmod);
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_fmod(&mut c);
    c.final_summary();
}
