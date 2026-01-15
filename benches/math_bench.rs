use criterion::{Criterion, black_box, criterion_group, criterion_main};
use fastlibm_bench::fastlibm;
use std::time::Duration;

fn bench_exp(c: &mut Criterion) {
    let mut group = c.benchmark_group("exp");
    let input = 1.23456f64;

    group.bench_function("fastlibm", |b| b.iter(|| fastlibm::exp(black_box(input))));
    group.bench_function("glibc", |b| b.iter(|| black_box(input).exp()));
    group.finish();
}

fn bench_ln(c: &mut Criterion) {
    let mut group = c.benchmark_group("ln");
    let input = 1.23456f64;

    group.bench_function("fastlibm", |b| b.iter(|| fastlibm::ln(black_box(input))));
    group.bench_function("glibc", |b| b.iter(|| black_box(input).ln()));
    group.finish();
}

fn bench_sin(c: &mut Criterion) {
    let mut group = c.benchmark_group("sin");
    let input = 1.23456f64;

    group.bench_function("fastlibm", |b| b.iter(|| fastlibm::sin(black_box(input))));
    group.bench_function("glibc", |b| b.iter(|| black_box(input).sin()));
    group.finish();
}

fn bench_cos(c: &mut Criterion) {
    let mut group = c.benchmark_group("cos");
    let input = 1.23456f64;

    group.bench_function("fastlibm", |b| b.iter(|| fastlibm::cos(black_box(input))));
    group.bench_function("glibc", |b| b.iter(|| black_box(input).cos()));
    group.finish();
}

fn configure_criterion() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(Duration::from_secs(10))
        .warm_up_time(Duration::from_secs(5))
}

criterion_group! {
    name = benches;
    config = configure_criterion();
    targets = bench_exp, bench_ln, bench_sin, bench_cos
}
criterion_main!(benches);
