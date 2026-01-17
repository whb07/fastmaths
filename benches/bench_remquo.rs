use criterion::{Criterion, black_box};
use fastmaths::fastlibm;

mod bench_util;
use bench_util::{configure_criterion, gen_pairs, glibc_remquo};

fn bench_remquo(c: &mut Criterion) {
    let inputs = vec![
        (5.3, 2.1),
        (-5.3, 2.1),
        (1.0, 0.5),
        (1.0, -0.5),
        (1e6, 3.0),
        (-1e6, 3.0),
    ];
    let mut common = gen_pairs(2048, -1e6, 1e6, 0x3233);
    for (_x, y) in &mut common {
        if *y == 0.0 {
            *y = 1.0;
        }
    }

    let mut group = c.benchmark_group("remquo/smoke");
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_q = 0i32;
            for &(x, y) in &inputs {
                let (r, q) = fastlibm::remquo(black_box(x), black_box(y));
                acc += r;
                acc_q ^= q;
            }
            black_box((acc, acc_q))
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_q = 0i32;
            for &(x, y) in &inputs {
                let (r, q) = glibc_remquo(black_box(x), black_box(y));
                acc += r;
                acc_q ^= q;
            }
            black_box((acc, acc_q))
        })
    });
    group.finish();

    let mut group = c.benchmark_group("remquo/common");
    group.bench_function("fastlibm", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_q = 0i32;
            for &(x, y) in &common {
                let (r, q) = fastlibm::remquo(black_box(x), black_box(y));
                acc += r;
                acc_q ^= q;
            }
            black_box((acc, acc_q))
        })
    });
    group.bench_function("glibc", |b| {
        b.iter(|| {
            let mut acc = 0.0;
            let mut acc_q = 0i32;
            for &(x, y) in &common {
                let (r, q) = glibc_remquo(black_box(x), black_box(y));
                acc += r;
                acc_q ^= q;
            }
            black_box((acc, acc_q))
        })
    });
    group.finish();
}

fn main() {
    let mut c = configure_criterion();
    bench_remquo(&mut c);
    c.final_summary();
}
