//! lgamma/tgamma implementation (core-math port).
//!
//! Implements lgamma via table-driven log/sinpi reductions, reflection formula,
//! and asymptotic expansions. Uses double-double arithmetic (ddcoremath) to keep
//! error below 1 ULP in difficult regions. Constants and tables are sourced from
//! glibc/core-math (see glibc/ directory).

use super::{asdouble, exp, expm1, fasttwosum, fma_internal, roundeven_finite};
use super::{floor, rint};

// === ddcoremath helpers (ported) ===

#[inline(always)]
fn fasttwosub(x: f64, y: f64, e: &mut f64) -> f64 {
    let s = x - y;
    let z = x - s;
    *e = z - y;
    s
}

#[inline(always)]
fn twosum(x: f64, y: f64, e: &mut f64) -> f64 {
    if x.abs() > y.abs() {
        fasttwosum(x, y, e)
    } else {
        fasttwosum(y, x, e)
    }
}

#[inline(always)]
fn fastsum(xh: f64, xl: f64, yh: f64, yl: f64, e: &mut f64) -> f64 {
    let mut sl = 0.0;
    let sh = fasttwosum(xh, yh, &mut sl);
    *e = (xl + yl) + sl;
    sh
}

#[inline(always)]
fn sumdd(xh: f64, xl: f64, yh: f64, yl: f64, e: &mut f64) -> f64 {
    let mut sl = 0.0;
    let sh = if xh.abs() > yh.abs() {
        fasttwosum(xh, yh, &mut sl)
    } else {
        fasttwosum(yh, xh, &mut sl)
    };
    *e = (xl + yl) + sl;
    sh
}

#[inline(always)]
fn adddd(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let s = xh + ch;
    let d = s - xh;
    *l = ((ch - d) + (xh + (d - s))) + (xl + cl);
    s
}

#[inline(always)]
fn muldd_acc(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahlh = ch * xl;
    let alhh = cl * xh;
    let ahhh = ch * xh;
    let mut ahhl = fma_internal(ch, xh, -ahhh);
    ahhl += alhh + ahlh;
    let chh = ahhh + ahhl;
    *l = (ahhh - chh) + ahhl;
    chh
}

#[inline(always)]
fn muldd_acc2(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahlh = ch * xl;
    let alhh = cl * xh;
    let ahhh = ch * xh;
    let mut ahhl = fma_internal(ch, xh, -ahhh);
    ahhl += alhh + ahlh;
    fasttwosum(ahhh, ahhl, l)
}

#[inline(always)]
fn muldd2(xh: f64, xl: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahhh = ch * xh;
    *l = (ch * xl + cl * xh) + fma_internal(ch, xh, -ahhh);
    ahhh
}

#[inline(always)]
fn muldd3(xh: f64, xl: f64, yh: f64, yl: f64, l: &mut f64) -> f64 {
    let ch = xh * yh;
    let cl1 = fma_internal(xh, yh, -ch);
    let tl0 = xl * yl;
    let tl1 = tl0 + xh * yl;
    let cl2 = tl1 + xl * yh;
    let cl3 = cl1 + cl2;
    fasttwosum(ch, cl3, l)
}

#[inline(always)]
fn mulddd(xh: f64, xl: f64, ch: f64, l: &mut f64) -> f64 {
    let ahlh = ch * xl;
    let ahhh = ch * xh;
    let mut ahhl = fma_internal(ch, xh, -ahhh);
    ahhl += ahlh;
    let chh = ahhh + ahhl;
    *l = (ahhh - chh) + ahhl;
    chh
}

#[inline(always)]
fn mulddd2(x: f64, ch: f64, cl: f64, l: &mut f64) -> f64 {
    let ahhh = ch * x;
    *l = cl * x + fma_internal(ch, x, -ahhh);
    ahhh
}

#[inline(always)]
fn mulddd3(xh: f64, xl: f64, ch: f64, l: &mut f64) -> f64 {
    let hh = xh * ch;
    *l = fma_internal(ch, xh, -hh) + xl * ch;
    hh
}

#[inline(always)]
fn polydd(xh: f64, xl: f64, n: usize, c: &[[f64; 2]], l: &mut f64) -> f64 {
    let mut i = n - 1;
    let mut cl = 0.0;
    let mut ch = fasttwosum(c[i][0], *l, &mut cl);
    cl += c[i][1];
    while i > 0 {
        i -= 1;
        ch = muldd_acc(xh, xl, ch, cl, &mut cl);
        let th = ch + c[i][0];
        let tl = (c[i][0] - th) + ch;
        ch = th;
        cl += tl + c[i][1];
    }
    *l = cl;
    ch
}

#[inline(always)]
fn polydd2(xh: f64, xl: f64, n: usize, c: &[[f64; 2]], l: &mut f64) -> f64 {
    let mut i = n - 1;
    let mut cl = 0.0;
    let mut ch = fasttwosum(c[i][0], *l, &mut cl);
    cl += c[i][1];
    while i > 0 {
        i -= 1;
        ch = muldd2(xh, xl, ch, cl, &mut cl);
        ch = fastsum(c[i][0], c[i][1], ch, cl, &mut cl);
    }
    *l = cl;
    ch
}

#[inline(always)]
fn polydd3(xh: f64, xl: f64, n: usize, c: &[[f64; 2]], l: &mut f64) -> f64 {
    let mut i = n - 1;
    let mut cl = 0.0;
    let mut ch = fasttwosum(c[i][0], *l, &mut cl);
    cl += c[i][1];
    while i > 0 {
        i -= 1;
        ch = muldd_acc2(xh, xl, ch, cl, &mut cl);
        let mut tl = 0.0;
        let th = fasttwosum(c[i][0], ch, &mut tl);
        ch = th;
        cl += tl + c[i][1];
    }
    *l = cl;
    ch
}

#[inline(always)]
fn polyddd(x: f64, n: usize, c: &[[f64; 2]], l: &mut f64) -> f64 {
    let mut i = n - 1;
    let mut cl = 0.0;
    let mut ch = fasttwosum(c[i][0], *l, &mut cl);
    cl += c[i][1];
    while i > 0 {
        i -= 1;
        ch = mulddd2(x, ch, cl, &mut cl);
        ch = sumdd(c[i][0], c[i][1], ch, cl, &mut cl);
    }
    *l = cl;
    ch
}

#[inline(always)]
fn polydddfst(x: f64, n: usize, c: &[[f64; 2]], l: &mut f64) -> f64 {
    let mut i = n - 1;
    let mut cl = 0.0;
    let mut ch = fasttwosum(c[i][0], *l, &mut cl);
    cl += c[i][1];
    while i > 0 {
        i -= 1;
        ch = mulddd2(x, ch, cl, &mut cl);
        ch = fastsum(c[i][0], c[i][1], ch, cl, &mut cl);
    }
    *l = cl;
    ch
}

#[inline(always)]
fn polyd(x: f64, n: usize, c: &[[f64; 2]]) -> f64 {
    let mut i = n - 1;
    let mut ch = c[i][0];
    while i > 0 {
        i -= 1;
        ch = c[i][0] + x * ch;
    }
    ch
}

const UBRD: [u32; 20] = [
    0x1ff0000, 0x1ff146c, 0x1ff2b7b, 0x1ff4532, 0x1ff614c, 0x1ff8310, 0x1ff93f7, 0x1ffa880,
    0x1ffc05e, 0x1ffdb73, 0x1fff8a5, 0x2001147, 0x2002703, 0x20041ac, 0x200622a, 0x20084d9,
    0x2009ce7, 0x200ba2c, 0x200ddd7, 0x20104ba,
];

#[derive(Clone, Copy)]
struct LogB {
    c0: u16,
    c1: i16,
}

const LOG_B: [LogB; 32] = [
    LogB { c0: 301, c1: 27565 },
    LogB {
        c0: 7189,
        c1: 24786,
    },
    LogB {
        c0: 13383,
        c1: 22167,
    },
    LogB {
        c0: 18923,
        c1: 19696,
    },
    LogB {
        c0: 23845,
        c1: 17361,
    },
    LogB {
        c0: 28184,
        c1: 15150,
    },
    LogB {
        c0: 31969,
        c1: 13054,
    },
    LogB {
        c0: 35231,
        c1: 11064,
    },
    LogB {
        c0: 37996,
        c1: 9173,
    },
    LogB {
        c0: 40288,
        c1: 7372,
    },
    LogB {
        c0: 42129,
        c1: 5657,
    },
    LogB {
        c0: 43542,
        c1: 4020,
    },
    LogB {
        c0: 44546,
        c1: 2457,
    },
    LogB { c0: 45160, c1: 962 },
    LogB {
        c0: 45399,
        c1: -468,
    },
    LogB {
        c0: 45281,
        c1: -1838,
    },
    LogB {
        c0: 44821,
        c1: -3151,
    },
    LogB {
        c0: 44032,
        c1: -4412,
    },
    LogB {
        c0: 42929,
        c1: -5622,
    },
    LogB {
        c0: 41522,
        c1: -6786,
    },
    LogB {
        c0: 39825,
        c1: -7905,
    },
    LogB {
        c0: 37848,
        c1: -8982,
    },
    LogB {
        c0: 35602,
        c1: -10020,
    },
    LogB {
        c0: 33097,
        c1: -11020,
    },
    LogB {
        c0: 30341,
        c1: -11985,
    },
    LogB {
        c0: 27345,
        c1: -12916,
    },
    LogB {
        c0: 24115,
        c1: -13816,
    },
    LogB {
        c0: 20661,
        c1: -14685,
    },
    LogB {
        c0: 16989,
        c1: -15526,
    },
    LogB {
        c0: 13107,
        c1: -16339,
    },
    LogB {
        c0: 9022,
        c1: -17126,
    },
    LogB {
        c0: 4740,
        c1: -17889,
    },
];

const LGAMMA_DB: [[f64; 3]; 19] = [
    [
        f64::from_bits(0xc021b649eb4316fb),
        f64::from_bits(0xc0250332a035af1f),
        f64::from_bits(0xbcc0000000000000),
    ],
    [
        f64::from_bits(0xc01808d3e2f56b4f),
        f64::from_bits(0xbffd779a9ab6cbff),
        f64::from_bits(0x3c90000000000000),
    ],
    [
        f64::from_bits(0xc00d02b2008d5bf5),
        f64::from_bits(0xbff67cf93db4863e),
        f64::from_bits(0xbc90000000000000),
    ],
    [
        f64::from_bits(0xc00a123d403647d7),
        f64::from_bits(0xbfe530824ff0740a),
        f64::from_bits(0x3c80000000000000),
    ],
    [
        f64::from_bits(0xc0090404f46978b6),
        f64::from_bits(0x3fc189e211cebce7),
        f64::from_bits(0xbc60000000000000),
    ],
    [
        f64::from_bits(0xc0075692939b3def),
        f64::from_bits(0x3fea11a8b4dd58bb),
        f64::from_bits(0x3c80000000000000),
    ],
    [
        f64::from_bits(0xc0064cc652d46a2d),
        f64::from_bits(0x3fb7c2fd07e26d78),
        f64::from_bits(0xbc50000000000000),
    ],
    [
        f64::from_bits(0xc00593ae533139c1),
        f64::from_bits(0xbfb316f7a9bb3e4e),
        f64::from_bits(0x3c50000000000000),
    ],
    [
        f64::from_bits(0xc00549cde4f1fd16),
        f64::from_bits(0xbfbab1280c638ada),
        f64::from_bits(0xbc50000000000000),
    ],
    [
        f64::from_bits(0xc0053712a3a51156),
        f64::from_bits(0xbfbbedf4564c976d),
        f64::from_bits(0x3c50000000000000),
    ],
    [
        f64::from_bits(0xc0051333e5a494f7),
        f64::from_bits(0xbfbd93a67f2ad756),
        f64::from_bits(0x3c50000000000000),
    ],
    [
        f64::from_bits(0xc00504f0b5b9d1f9),
        f64::from_bits(0xbfbdfa224bc3435e),
        f64::from_bits(0xbc50000000000000),
    ],
    [
        f64::from_bits(0xc004c213b7fa02dc),
        f64::from_bits(0xbfbe051ae181baab),
        f64::from_bits(0x3c50000000000000),
    ],
    [
        f64::from_bits(0xbff75824f4f0c7c8),
        f64::from_bits(0x3fecb23f1fc0296e),
        f64::from_bits(0x3c80000000000000),
    ],
    [
        f64::from_bits(0xbff667cb87bcfa0f),
        f64::from_bits(0x3fef478a478204c5),
        f64::from_bits(0x3c80000000000000),
    ],
    [
        f64::from_bits(0xbff51668b6af122d),
        f64::from_bits(0x3ff2728422577e00),
        f64::from_bits(0x3c90000000000000),
    ],
    [
        f64::from_bits(0xbff50b41410187c1),
        f64::from_bits(0x3ff290131f9d5fea),
        f64::from_bits(0xbc90000000000000),
    ],
    [
        f64::from_bits(0xbfd307eb80d114af),
        f64::from_bits(0x3ff7870d113b3feb),
        f64::from_bits(0xbc90000000000000),
    ],
    [
        f64::from_bits(0xbfd26923ac1f7c10),
        f64::from_bits(0x3ff7df2c32cf08a3),
        f64::from_bits(0x3c90000000000000),
    ],
];

const LGAMMA_C0: [[f64; 2]; 34] = [
    [
        f64::from_bits(0xbfe2788cfc6fb619),
        f64::from_bits(0x3c56cb90701fbf92),
    ],
    [
        f64::from_bits(0x3fea51a6625307d3),
        f64::from_bits(0x3c71873d891220df),
    ],
    [
        f64::from_bits(0xbfd9a4d55beab2d7),
        f64::from_bits(0x3c44c26d1b4d5994),
    ],
    [
        f64::from_bits(0x3fd151322ac7d848),
        f64::from_bits(0x3c6b5f9120fdaca6),
    ],
    [
        f64::from_bits(0xbfca8b9c17aa6149),
        f64::from_bits(0xbc52e826b9f2e3f2),
    ],
    [
        f64::from_bits(0x3fc5b40cb100c306),
        f64::from_bits(0x3c44a79a5f96e16f),
    ],
    [
        f64::from_bits(0xbfc2703a1dcea3ae),
        f64::from_bits(0xbc663066fbe425d6),
    ],
    [
        f64::from_bits(0x3fc010b36af86397),
        f64::from_bits(0xbc47438c20a423a8),
    ],
    [
        f64::from_bits(0xbfbc806706d57db4),
        f64::from_bits(0xbc55a891c586f017),
    ],
    [
        f64::from_bits(0x3fb9a01e385d5f8f),
        f64::from_bits(0x3c49eb2a3ca8e809),
    ],
    [
        f64::from_bits(0xbfb748c33114c6d5),
        f64::from_bits(0xbc5505ae0736c854),
    ],
    [
        f64::from_bits(0x3fb556ad63243bc2),
        f64::from_bits(0xbc50cfd5c4f86f62),
    ],
    [
        f64::from_bits(0xbfb3b1d971fc59e2),
        f64::from_bits(0xbc29a19a5ba46076),
    ],
    [
        f64::from_bits(0x3fb2496df8320d56),
        f64::from_bits(0x3c551588f44e1564),
    ],
    [
        f64::from_bits(0xbfb11133476e5f97),
        f64::from_bits(0xbc4fb27c7d952db7),
    ],
    [
        f64::from_bits(0x3fb00010064c948a),
        f64::from_bits(0x3c3c2a9309d327a6),
    ],
    [
        f64::from_bits(0xbfae1e2d312e9ab0),
        f64::from_bits(0x3c4beee6edc135cf),
    ],
    [
        f64::from_bits(0x3fac71ce3a414da2),
        f64::from_bits(0x3c17d6e71837bc1d),
    ],
    [
        f64::from_bits(0xbfaaf28a18673c39),
        f64::from_bits(0x3c3df2c138661f9b),
    ],
    [
        f64::from_bits(0x3fa9999b2dfca9bb),
        f64::from_bits(0xbc4afb6a810fa7c2),
    ],
    [
        f64::from_bits(0xbfa861874180c0f2),
        f64::from_bits(0xbc4c7900b8893e2b),
    ],
    [
        f64::from_bits(0x3fa745d2798dac6a),
        f64::from_bits(0x3c2e962ae5d1c788),
    ],
    [
        f64::from_bits(0xbfa642be2fd1947d),
        f64::from_bits(0x3c294490a4e3d6b4),
    ],
    [
        f64::from_bits(0x3fa55545da95a50b),
        f64::from_bits(0x3c0b56c4cc4dd96b),
    ],
    [
        f64::from_bits(0xbfa47ba80f965799),
        f64::from_bits(0x3c33bcac9fdf0d1a),
    ],
    [
        f64::from_bits(0x3fa3b24eb343ccef),
        f64::from_bits(0xbc283fb9b692b013),
    ],
    [
        f64::from_bits(0xbfa2eba18fc96750),
        f64::from_bits(0xbc33069a0cc7e12d),
    ],
    [
        f64::from_bits(0x3fa23b27295dc8b4),
        f64::from_bits(0x3c3a6d42b2188101),
    ],
    [
        f64::from_bits(0xbfa2136ba2a595c4),
        f64::from_bits(0xbc4c2acae5604f9d),
    ],
    [
        f64::from_bits(0x3fa191f0bce33ba9),
        f64::from_bits(0xbc28470cbefa841a),
    ],
    [
        f64::from_bits(0xbf9b86662fc0fb81),
        f64::from_bits(0xbc3236cd67b781c2),
    ],
    [
        f64::from_bits(0x3f99d6d718c336ab),
        f64::from_bits(0x3c3c632d9ae925d4),
    ],
    [
        f64::from_bits(0xbfa9ea8d7c263a80),
        f64::from_bits(0xbc4bdcf945212655),
    ],
    [
        f64::from_bits(0x3fa9f461fd74cc30),
        f64::from_bits(0x3c35e5f99eb89ba4),
    ],
];

const LGAMMA_B: [[f64; 2]; 30] = [
    [
        f64::from_bits(0xbfbeeb95b094c191),
        f64::from_bits(0xbc5346863f58b074),
    ],
    [
        f64::from_bits(0x3fa2aed059bd608a),
        f64::from_bits(0x3c0cd3d2ca77b58c),
    ],
    [
        f64::from_bits(0x3fdde9e64df22ef3),
        f64::from_bits(0xbc66d48ec99346d4),
    ],
    [
        f64::from_bits(0xbfc1ae55b180726c),
        f64::from_bits(0xbc4959aeebbdcd5b),
    ],
    [
        f64::from_bits(0x3fae0f840dad61da),
        f64::from_bits(0xbc4599fc3f5e3cc0),
    ],
    [
        f64::from_bits(0xbf9da59d5374a543),
        f64::from_bits(0xbc00628c37280b13),
    ],
    [
        f64::from_bits(0x3f8f9ca39daa929c),
        f64::from_bits(0xbbc69f468b73dadd),
    ],
    [
        f64::from_bits(0xbf81a8ba4f0ea597),
        f64::from_bits(0xbc27d1edde63f06c),
    ],
    [
        f64::from_bits(0x3f7456f1ad666a3b),
        f64::from_bits(0xbc131d73594fa521),
    ],
    [
        f64::from_bits(0xbf67edb812f6426e),
        f64::from_bits(0xbbf5f45d321e0ae2),
    ],
    [
        f64::from_bits(0x3f5c9735ae9db2c1),
        f64::from_bits(0xbbf827b962d41f07),
    ],
    [
        f64::from_bits(0xbf5148a319eec639),
        f64::from_bits(0x3bfa64298294f3a5),
    ],
    [
        f64::from_bits(0x3f4517c5a1579f40),
        f64::from_bits(0xbbc055841dae38a6),
    ],
    [
        f64::from_bits(0xbf39eff1d1c8be2d),
        f64::from_bits(0xbbcfba56c797e258),
    ],
    [
        f64::from_bits(0x3f300c41c13e3352),
        f64::from_bits(0xbb69c7d468e4cbd4),
    ],
    [
        f64::from_bits(0xbf23f6dff22a8ffe),
        f64::from_bits(0x3bb79e9e76cea415),
    ],
    [
        f64::from_bits(0x3f18f36195536e0c),
        f64::from_bits(0xbb99de6ad54f6997),
    ],
    [
        f64::from_bits(0xbf0f4ea079eaf734),
        f64::from_bits(0xbbafed2cdd5ce819),
    ],
    [
        f64::from_bits(0x3f03b5e73a6bf0ce),
        f64::from_bits(0xbb9b8bf8cedf2f5d),
    ],
    [
        f64::from_bits(0xbef8e583400720b6),
        f64::from_bits(0xbb958007f65918dc),
    ],
    [
        f64::from_bits(0x3eef88ed0371a528),
        f64::from_bits(0x3b4ca095e5525ba6),
    ],
    [
        f64::from_bits(0xbee40597e1beb94d),
        f64::from_bits(0x3b899bf6a1da80a1),
    ],
    [
        f64::from_bits(0x3ed97b2f3649dc6a),
        f64::from_bits(0x3b739d919c548f05),
    ],
    [
        f64::from_bits(0xbed03fa905165ffe),
        f64::from_bits(0x3b70e61588355356),
    ],
    [
        f64::from_bits(0x3ec4c8df9a3d3680),
        f64::from_bits(0xbb67a6689d0b7dd2),
    ],
    [
        f64::from_bits(0xbeba9bafae87dbd3),
        f64::from_bits(0x3b55777bb7ae251a),
    ],
    [
        f64::from_bits(0x3eb0b2f2a9e11fe8),
        f64::from_bits(0xbb4cf209a5b02c67),
    ],
    [
        f64::from_bits(0xbea567358802c384),
        f64::from_bits(0x3b4c44014666f1c9),
    ],
    [
        f64::from_bits(0x3ea1219755496b17),
        f64::from_bits(0xbb4b113ac02e4428),
    ],
    [
        f64::from_bits(0xbe963804ffdcd648),
        f64::from_bits(0x3b3ba1993b4ec9b7),
    ],
];

const LGAMMA_OFFS: [f64; 19] = [
    f64::from_bits(0x3fe146cd80000000),
    f64::from_bits(0x3fe3fe8980000000),
    f64::from_bits(0x3fe70aea80000000),
    f64::from_bits(0x3fea67fc00000000),
    f64::from_bits(0x3fee76db80000000),
    f64::from_bits(0x3ff1708380000000),
    f64::from_bits(0x3ff3c78a80000000),
    f64::from_bits(0x3ff68df200000000),
    f64::from_bits(0x3ff9bd1400000000),
    f64::from_bits(0x3ffd418680000000),
    f64::from_bits(0x4000d9a640000000),
    f64::from_bits(0x400384b800000000),
    f64::from_bits(0x40068b0600000000),
    f64::from_bits(0x400a3d6d00000000),
    f64::from_bits(0x400ebdd900000000),
    f64::from_bits(0x40121c1000000000),
    f64::from_bits(0x4015713680000000),
    f64::from_bits(0x4019803e80000000),
    f64::from_bits(0x401e74cc80000000),
];

const LGAMMA_CL: [[f64; 8]; 19] = [
    [
        f64::from_bits(0xc0118ad63ca097e9),
        f64::from_bits(0x401af8e15b715c51),
        f64::from_bits(0xc0256213b7191ba4),
        f64::from_bits(0x403151f165a9425f),
        f64::from_bits(0xc03c826426e4b7cd),
        f64::from_bits(0x4047c313095e4b75),
        f64::from_bits(0xc0544f3d7d848e78),
        f64::from_bits(0x40613384c97ea99d),
    ],
    [
        f64::from_bits(0xc000f58e76c8d235),
        f64::from_bits(0x40067c3f6b7124f6),
        f64::from_bits(0xc00ec78d7d8185a3),
        f64::from_bits(0x401588d6487de574),
        f64::from_bits(0xc01e9fbe8564220d),
        f64::from_bits(0x40260dd913b80b5e),
        f64::from_bits(0xc030465db7c895a6),
        f64::from_bits(0x4037ca34d903fc30),
    ],
    [
        f64::from_bits(0xbff0c505555a86b2),
        f64::from_bits(0x3ff33d3a22d1bb51),
        f64::from_bits(0xbff6d2a2457f05d4),
        f64::from_bits(0x3ffbb1c77fad8b03),
        f64::from_bits(0xc00115210a553746),
        f64::from_bits(0x400558e305fd6940),
        f64::from_bits(0xc00b4f9a0654679f),
        f64::from_bits(0x4011489e269cbf39),
    ],
    [
        f64::from_bits(0xbfe1170ead9585ba),
        f64::from_bits(0x3fe10c67b04495d7),
        f64::from_bits(0xbfe19de3af9dd349),
        f64::from_bits(0x3fe2a34cd6e66472),
        f64::from_bits(0xbfe40e2f93066eb8),
        f64::from_bits(0x3fe5ddc21a559735),
        f64::from_bits(0xbfe860a742d837aa),
        f64::from_bits(0x3feace7238771e7f),
    ],
    [
        f64::from_bits(0xbfcbe31df8f7d605),
        f64::from_bits(0x3fc8e33a32c94cf7),
        f64::from_bits(0xbfc6c62efd534bf3),
        f64::from_bits(0x3fc53719e404a7d6),
        f64::from_bits(0xbfc4074d1d083331),
        f64::from_bits(0x3fc31c6c5226f2b3),
        f64::from_bits(0xbfc2b062b8eedd9f),
        f64::from_bits(0x3fc219431f82fbfc),
    ],
    [
        f64::from_bits(0xbfc168e45409b785),
        f64::from_bits(0x3fba04e5759477f0),
        f64::from_bits(0xbfb43c620bb1d77f),
        f64::from_bits(0x3fb027d79414ff7d),
        f64::from_bits(0xbfaa46afc0776356),
        f64::from_bits(0x3fa5a92c3ddb75f5),
        f64::from_bits(0xbfa23d37b3b3e3b6),
        f64::from_bits(0x3f9f66a6169fe8ef),
    ],
    [
        f64::from_bits(0xbfb2db051283fb7a),
        f64::from_bits(0x3fa8afce072c9222),
        f64::from_bits(0xbfa0dcc84e49a658),
        f64::from_bits(0x3f97af09a263459b),
        f64::from_bits(0xbf90f51524188551),
        f64::from_bits(0x3f88a1f8bbd73d04),
        f64::from_bits(0xbf824e389a08ab23),
        f64::from_bits(0x3f7b5c2228d32783),
    ],
    [
        f64::from_bits(0xbfa3fbb4a9e75e6d),
        f64::from_bits(0x3f96c40332da72c0),
        f64::from_bits(0xbf8b235e2a6ed724),
        f64::from_bits(0x3f80a9023dcf81a5),
        f64::from_bits(0xbf74e128ec28e27a),
        f64::from_bits(0x3f6a90e4421d59e7),
        f64::from_bits(0xbf614e6becdb3889),
        f64::from_bits(0x3f568e6a1727f763),
    ],
    [
        f64::from_bits(0xbf95365e61675f08),
        f64::from_bits(0x3f84fd143859dc2c),
        f64::from_bits(0xbf75cb1d911bcabb),
        f64::from_bits(0x3f675a869d793508),
        f64::from_bits(0xbf59941834996ea5),
        f64::from_bits(0x3f4c782075b40e2c),
        f64::from_bits(0xbf403ab2e70c9df2),
        f64::from_bits(0x3f326c9636a06719),
    ],
    [
        f64::from_bits(0xbf87145b3bd2da75),
        f64::from_bits(0x3f73e6749a0fe630),
        f64::from_bits(0xbf620ea7ae3d0208),
        f64::from_bits(0x3f50f17a25bb48f3),
        f64::from_bits(0xbf4045c4de3c2101),
        f64::from_bits(0x3f2fcc9430345441),
        f64::from_bits(0xbf1fccc917990854),
        f64::from_bits(0x3f0f41fbb5026b50),
    ],
    [
        f64::from_bits(0x3f6253f3fc844189),
        f64::from_bits(0xbf4cadf5cc04da1b),
        f64::from_bits(0x3f373e9dbf6ed988),
        f64::from_bits(0xbf234f75abb9acfd),
        f64::from_bits(0x3f105502104fc072),
        f64::from_bits(0xbefc015daee9145b),
        f64::from_bits(0x3ee8af9ccda4578c),
        f64::from_bits(0xbed5b973bad98b6b),
    ],
    [
        f64::from_bits(0xbf67f80bfa6d705e),
        f64::from_bits(0x3f4e416c7e5d3bb3),
        f64::from_bits(0xbf34361a69711e0f),
        f64::from_bits(0x3f1c0beed7451d56),
        f64::from_bits(0xbf03fc0c220552b8),
        f64::from_bits(0x3eed09a0850c9ad6),
        f64::from_bits(0xbed5b403dca46450),
        f64::from_bits(0x3ec07c1349c4989a),
    ],
    [
        f64::from_bits(0xbf58bd8d36b6b680),
        f64::from_bits(0x3f3ab590101636b6),
        f64::from_bits(0xbf1e980083cba776),
        f64::from_bits(0x3f023c1bb53d55c4),
        f64::from_bits(0xbee65be06922ac50),
        f64::from_bits(0x3ecbfd7d06279b09),
        f64::from_bits(0xbeb2127a54c7c981),
        f64::from_bits(0x3e979ae1f9c24de8),
    ],
    [
        f64::from_bits(0xbf48db1211cc179c),
        f64::from_bits(0x3f26c24488bdc8e9),
        f64::from_bits(0xbf062882b2ca3c96),
        f64::from_bits(0x3ee67e36a9a5f89c),
        f64::from_bits(0xbec785b8294e4cf2),
        f64::from_bits(0x3ea925c40ccc4611),
        f64::from_bits(0xbe8bccb0b78110bc),
        f64::from_bits(0x3e6f05d365676624),
    ],
    [
        f64::from_bits(0xbf3879379df4c28c),
        f64::from_bits(0x3f12e193030ccfd1),
        f64::from_bits(0xbeef0900c0b3c1fc),
        f64::from_bits(0x3ecaa304a80f3ce4),
        f64::from_bits(0xbea795cdc082b2df),
        f64::from_bits(0x3e856025546a45a8),
        f64::from_bits(0xbe64137abecddaa9),
        f64::from_bits(0x3e4303d852294977),
    ],
    [
        f64::from_bits(0xbf27b7a8dbd38635),
        f64::from_bits(0x3efeab932219e072),
        f64::from_bits(0xbed528291d0efb42),
        f64::from_bits(0x3eae86163ded7066),
        f64::from_bits(0xbe86be3b6446506a),
        f64::from_bits(0x3e615d4d64fd204a),
        f64::from_bits(0xbe3b8820763832ef),
        f64::from_bits(0x3e15ff781e6e4e19),
    ],
    [
        f64::from_bits(0xbf16b1f37f261621),
        f64::from_bits(0x3ee87d88455d6443),
        f64::from_bits(0xbebc3a69901a7d7c),
        f64::from_bits(0x3e9107e1d8456499),
        f64::from_bits(0xbe653f4b8e0d8be8),
        f64::from_bits(0x3e3b3016ba9fadff),
        f64::from_bits(0xbe12175c3eb445d8),
        f64::from_bits(0x3de841e8df7f84e7),
    ],
    [
        f64::from_bits(0xbf057cce7fdc9fe7),
        f64::from_bits(0x3ed347c6b65ace16),
        f64::from_bits(0xbea27eb9df26d911),
        f64::from_bits(0x3e7296c0dcf0b476),
        f64::from_bits(0xbe4354fd38659786),
        f64::from_bits(0x3e14a2dbe1c4af19),
        f64::from_bits(0xbde6f1651636db10),
        f64::from_bits(0x3db9b1457d56445a),
    ],
    [
        f64::from_bits(0xbef423d487c99e54),
        f64::from_bits(0x3ebdf4d1f34022a3),
        f64::from_bits(0xbe87d54e9cd7e7ea),
        f64::from_bits(0x3e53e131c44c6382),
        f64::from_bits(0xbe212afb6bfa8c14),
        f64::from_bits(0x3dee7412bd9ebd87),
        f64::from_bits(0xbdbc2ab005ebc13b),
        f64::from_bits(0x3d8a3bbacb6ee6b7),
    ],
];

const LGAMMA_CH: [[[f64; 2]; 5]; 19] = [
    [
        [
            f64::from_bits(0x3fdfdbd7c56b02b5),
            f64::from_bits(0xbc79f8c66985b6f3),
        ],
        [
            f64::from_bits(0xbffc771ed8981f3e),
            f64::from_bits(0x3c98d8b72ce9b19d),
        ],
        [
            f64::from_bits(0x4001558ba7c0144d),
            f64::from_bits(0x3ca4fc1fa0f0451c),
        ],
        [
            f64::from_bits(0xc001fa938f4d4b53),
            f64::from_bits(0xbcaf29beb3ca3738),
        ],
        [
            f64::from_bits(0x4007f7469f6781ef),
            f64::from_bits(0xbcab59ce1aa03545),
        ],
    ],
    [
        [
            f64::from_bits(0x3fd71c14e711391e),
            f64::from_bits(0x3c32ad5eb4fb4f59),
        ],
        [
            f64::from_bits(0xbff740c890bd54d3),
            f64::from_bits(0xbc86978dab8a1160),
        ],
        [
            f64::from_bits(0x3ffb38de2e957c18),
            f64::from_bits(0xbc8aba2b91749902),
        ],
        [
            f64::from_bits(0xbff7ab358c51c087),
            f64::from_bits(0x3c846a8f1bc5883b),
        ],
        [
            f64::from_bits(0x3ffaf1b63b322b6d),
            f64::from_bits(0x3c82d98d261df8f3),
        ],
    ],
    [
        [
            f64::from_bits(0x3fce53b12b3407e2),
            f64::from_bits(0x3c697cb2965d31b5),
        ],
        [
            f64::from_bits(0xbff2a144e9a8b92e),
            f64::from_bits(0xbc9bbf90d2717ba5),
        ],
        [
            f64::from_bits(0x3ff5adc4ef58621e),
            f64::from_bits(0x3c9d41b3282f1d5b),
        ],
        [
            f64::from_bits(0xbfefb259e2817239),
            f64::from_bits(0x3c8a19b744867ccb),
        ],
        [
            f64::from_bits(0x3feee43256a6bfd3),
            f64::from_bits(0x3c8880c7ca4d6687),
        ],
    ],
    [
        [
            f64::from_bits(0x3fc0719312af823c),
            f64::from_bits(0xbc677ca1d8b99601),
        ],
        [
            f64::from_bits(0xbfed11f75dc5be7d),
            f64::from_bits(0x3c6997295e7f58d5),
        ],
        [
            f64::from_bits(0x3ff18a58180335dd),
            f64::from_bits(0xbc59e4f675e9e244),
        ],
        [
            f64::from_bits(0xbfe5aea0e9166a08),
            f64::from_bits(0x3c84799eb996a78b),
        ],
        [
            f64::from_bits(0x3fe22b448094c052),
            f64::from_bits(0xbc7221db12561423),
        ],
    ],
    [
        [
            f64::from_bits(0x3fe3c3b637596f8d),
            f64::from_bits(0xbc7051b18f5744ba),
        ],
        [
            f64::from_bits(0xbfeb9ccef0d71197),
            f64::from_bits(0xbc8cf98e73bfb3d7),
        ],
        [
            f64::from_bits(0x3fdc55517304ef35),
            f64::from_bits(0x3c6dfe2299217a1a),
        ],
        [
            f64::from_bits(0xbfd4230fb2a20b13),
            f64::from_bits(0xbc68eb1c5690348f),
        ],
        [
            f64::from_bits(0x3fd03aa1691c1841),
            f64::from_bits(0x3c6a0e14e4b5a96c),
        ],
    ],
    [
        [
            f64::from_bits(0xbfa752403c835a4d),
            f64::from_bits(0x3c4a3a43faf6eccc),
        ],
        [
            f64::from_bits(0xbfdc0be76051e3a5),
            f64::from_bits(0xbc6c737cd3ea73d9),
        ],
        [
            f64::from_bits(0x3fe73c36ef7bf402),
            f64::from_bits(0xbc740c4dff8e4c1e),
        ],
        [
            f64::from_bits(0xbfd458cec1d1393d),
            f64::from_bits(0xbc7c7f148cf356ef),
        ],
        [
            f64::from_bits(0x3fc8ec2d305516c4),
            f64::from_bits(0x3c69566535c9eab0),
        ],
    ],
    [
        [
            f64::from_bits(0xbfb85361b993719f),
            f64::from_bits(0xbc5dc41ac35a716f),
        ],
        [
            f64::from_bits(0xbfcf3e2bae2cdf7d),
            f64::from_bits(0x3c66d5cae27956a4),
        ],
        [
            f64::from_bits(0x3fe3745220b46975),
            f64::from_bits(0x3c356d68f9018bb8),
        ],
        [
            f64::from_bits(0xbfcd29172b1a4407),
            f64::from_bits(0xbc593fc4238117bd),
        ],
        [
            f64::from_bits(0x3fbef0f914e4a75b),
            f64::from_bits(0xbc56f0339a5cbb3a),
        ],
    ],
    [
        [
            f64::from_bits(0xbfbec2ab5aa5843a),
            f64::from_bits(0xbc1adc658df2c1c1),
        ],
        [
            f64::from_bits(0xbfaa6243a7f3534c),
            f64::from_bits(0x3c40dc0b707b85ab),
        ],
        [
            f64::from_bits(0x3fe04116f85f23a3),
            f64::from_bits(0x3c6517c0b25b9233),
        ],
        [
            f64::from_bits(0xbfc4bb33f1abe408),
            f64::from_bits(0x3c5cc0c1f637cea4),
        ],
        [
            f64::from_bits(0x3fb2ecfafae59f8f),
            f64::from_bits(0x3c53c57c7651ae8a),
        ],
    ],
    [
        [
            f64::from_bits(0xbfbc8928613eb4f5),
            f64::from_bits(0x3c155f36a43c02bc),
        ],
        [
            f64::from_bits(0x3fc1151b40dad4e9),
            f64::from_bits(0x3c667907a753aa66),
        ],
        [
            f64::from_bits(0x3fdb46b0b78660ac),
            f64::from_bits(0xbc79bcdfa3bbcd41),
        ],
        [
            f64::from_bits(0xbfbd9cd6009ac89d),
            f64::from_bits(0xbc469c4d18a5c993),
        ],
        [
            f64::from_bits(0x3fa73b079d35c37c),
            f64::from_bits(0xbc44d3891ecef09e),
        ],
    ],
    [
        [
            f64::from_bits(0xbfb00ad2093da6e4),
            f64::from_bits(0xbc5cbf7cf8850330),
        ],
        [
            f64::from_bits(0x3fd391f431d39831),
            f64::from_bits(0x3c78fb94bb0e7df5),
        ],
        [
            f64::from_bits(0x3fd71d5a6e677f1c),
            f64::from_bits(0x3c4d1dc12aaa3806),
        ],
        [
            f64::from_bits(0xbfb57f6fbf9108c1),
            f64::from_bits(0x3c24e341fb4cef78),
        ],
        [
            f64::from_bits(0x3f9d1e33efae7a1d),
            f64::from_bits(0x3c3c4938a6deffbe),
        ],
    ],
    [
        [
            f64::from_bits(0x3fdd344dabcc201e),
            f64::from_bits(0x3c7574f453e55614),
        ],
        [
            f64::from_bits(0x3fd3c3a02b015763),
            f64::from_bits(0xbc7342e3d6a27dfa),
        ],
        [
            f64::from_bits(0xbfaf5d49f62ecfd6),
            f64::from_bits(0xbc307444b43ab601),
        ],
        [
            f64::from_bits(0x3f922abe7bbdf628),
            f64::from_bits(0x3c02cb184651725a),
        ],
        [
            f64::from_bits(0xbf78b52066552f48),
            f64::from_bits(0x3c1bc2dbb1b8365d),
        ],
    ],
    [
        [
            f64::from_bits(0x3fcf22e8b160e053),
            f64::from_bits(0x3c689c03c62a66d7),
        ],
        [
            f64::from_bits(0x3fe58ae0ae321620),
            f64::from_bits(0xbc7594df075ee813),
        ],
        [
            f64::from_bits(0x3fd028e87f2859fd),
            f64::from_bits(0x3c5bf1ead4dde3d4),
        ],
        [
            f64::from_bits(0xbfa55b4949f3971a),
            f64::from_bits(0xbc42cfd594571487),
        ],
        [
            f64::from_bits(0x3f84cfe08a2baa09),
            f64::from_bits(0x3c1495ab3aeecaf0),
        ],
    ],
    [
        [
            f64::from_bits(0x3fe104861734d948),
            f64::from_bits(0x3c732e74856dbad8),
        ],
        [
            f64::from_bits(0x3feb2445e9d82006),
            f64::from_bits(0x3c86e48e474ddfbf),
        ],
        [
            f64::from_bits(0x3fcb352d20042182),
            f64::from_bits(0x3c3a8ac4f9b7c938),
        ],
        [
            f64::from_bits(0xbf9e6c5b3585790e),
            f64::from_bits(0xbc331a8ef26cbf2e),
        ],
        [
            f64::from_bits(0x3f793111b206dab4),
            f64::from_bits(0xbc0aa3ae79b17070),
        ],
    ],
    [
        [
            f64::from_bits(0x3feeed49cf014c0b),
            f64::from_bits(0x3c8bca14c01f79ae),
        ],
        [
            f64::from_bits(0x3ff0718fe597659b),
            f64::from_bits(0x3c87d14012138c17),
        ],
        [
            f64::from_bits(0x3fc6c89e19ff8e58),
            f64::from_bits(0x3c412dfe29d6e296),
        ],
        [
            f64::from_bits(0xbf956a9890298c3a),
            f64::from_bits(0xbc22181516eb15d6),
        ],
        [
            f64::from_bits(0x3f6deaa0ec93f6d9),
            f64::from_bits(0x3c0e4d7a3e816168),
        ],
    ],
    [
        [
            f64::from_bits(0x3ff990530fe5fa37),
            f64::from_bits(0xbc7cf639a3a54f76),
        ],
        [
            f64::from_bits(0x3ff35e029ece68d0),
            f64::from_bits(0xbc3e3db2cbb514eb),
        ],
        [
            f64::from_bits(0x3fc301f23426a05f),
            f64::from_bits(0xbc6b5ec346a456bc),
        ],
        [
            f64::from_bits(0xbf8de5b0dd5127b5),
            f64::from_bits(0x3c2371374acf777f),
        ],
        [
            f64::from_bits(0x3f61843ded6af0f6),
            f64::from_bits(0xbbf7779056d71400),
        ],
    ],
    [
        [
            f64::from_bits(0x4003ef64cb5ced7b),
            f64::from_bits(0x3c9c3c21b0562715),
        ],
        [
            f64::from_bits(0x3ff654a3f497c726),
            f64::from_bits(0x3c93331f28ee09bb),
        ],
        [
            f64::from_bits(0x3fbf9f5117f295a1),
            f64::from_bits(0xbc5ee3d2bb334106),
        ],
        [
            f64::from_bits(0xbf84bb07b47ebf8d),
            f64::from_bits(0xbc264c2c019b90b5),
        ],
        [
            f64::from_bits(0x3f5449d9854bac59),
            f64::from_bits(0xbbfd0a2827bf2270),
        ],
    ],
    [
        [
            f64::from_bits(0x400de185c1178ad9),
            f64::from_bits(0x3c8d477f1a273bfc),
        ],
        [
            f64::from_bits(0x3ff9539397e34b21),
            f64::from_bits(0x3c99743cc0cd10f2),
        ],
        [
            f64::from_bits(0x3fba3e2c09f7886d),
            f64::from_bits(0xbc417f6c25e05338),
        ],
        [
            f64::from_bits(0xbf7c98eb5fc97ce2),
            f64::from_bits(0x3c00a5104a9f402d),
        ],
        [
            f64::from_bits(0x3f474b50213890ab),
            f64::from_bits(0x3beff0ae56647ad0),
        ],
    ],
    [
        [
            f64::from_bits(0x4015c2be39a4c6fd),
            f64::from_bits(0x3cbff2814687494c),
        ],
        [
            f64::from_bits(0x3ffc59e5d40889c7),
            f64::from_bits(0x3c8299ee0827992a),
        ],
        [
            f64::from_bits(0x3fb5bbf97b18270e),
            f64::from_bits(0xbc3d04ddc6346897),
        ],
        [
            f64::from_bits(0xbf73a2d0322cf70e),
            f64::from_bits(0x3be53fe131154027),
        ],
        [
            f64::from_bits(0x3f3a8c6d657c0cfd),
            f64::from_bits(0xbbdb402fb82b45ef),
        ],
    ],
    [
        [
            f64::from_bits(0x401f07834a362b11),
            f64::from_bits(0xbcb738a86a953af8),
        ],
        [
            f64::from_bits(0x3fff68034cafc0d3),
            f64::from_bits(0x3c7b8d6c9e2cd7d4),
        ],
        [
            f64::from_bits(0x3fb1f68e6efd00fa),
            f64::from_bits(0xbc26083738e28e87),
        ],
        [
            f64::from_bits(0xbf6ad889b8da1552),
            f64::from_bits(0x3bf1325e8a48689d),
        ],
        [
            f64::from_bits(0x3f2e0ae44f526429),
            f64::from_bits(0xbbc997df9412e4aa),
        ],
    ],
];

const LGAMMA_SMALL_C0: [[f64; 2]; 4] = [
    [
        f64::from_bits(0xbfe2788cfc6fb619),
        f64::from_bits(0x3c56cb9a4ff7c53b),
    ],
    [
        f64::from_bits(0x3fea51a6625307d3),
        f64::from_bits(0x3c718722054895e9),
    ],
    [
        f64::from_bits(0xbfd9a4d55beab2d7),
        f64::from_bits(0xbc074ded0474fe66),
    ],
    [
        f64::from_bits(0x3fd151322ac7d848),
        f64::from_bits(0x3c7825b3df1d5722),
    ],
];

const LGAMMA_SMALL_Q: [f64; 8] = [
    f64::from_bits(0xbfca8b9c17aa5d3d),
    f64::from_bits(0x3fc5b40cb100b9bf),
    f64::from_bits(0xbfc2703a1e13bcbc),
    f64::from_bits(0x3fc010b36b6afdc1),
    f64::from_bits(0xbfbc8062dd09ec62),
    f64::from_bits(0x3fb9a018c7345316),
    f64::from_bits(0xbfb7578ea8068cc4),
    f64::from_bits(0x3fb566b51c990008),
];

const LGAMMA_ASYM_C: [[f64; 2]; 2] = [
    [
        f64::from_bits(0x3fdacfe390c97d6a),
        f64::from_bits(0xbc51d9792ced423a),
    ],
    [
        f64::from_bits(0x3fb55555555554c1),
        f64::from_bits(0xbc40143af34001bd),
    ],
];

const LGAMMA_ASYM_Q: [f64; 5] = [
    f64::from_bits(0xbf66c16c1697de08),
    f64::from_bits(0x3f4a019f47b230fd),
    f64::from_bits(0xbf4380aab821e42e),
    f64::from_bits(0x3f4b617d2c5b5b66),
    f64::from_bits(0xbf5a7fd66a05ccfc),
];

const LOG_R1: [f64; 33] = [
    f64::from_bits(0x3ff0000000000000),
    f64::from_bits(0x3fef508000000000),
    f64::from_bits(0x3feea4a000000000),
    f64::from_bits(0x3fedfca000000000),
    f64::from_bits(0x3fed582000000000),
    f64::from_bits(0x3fecb72000000000),
    f64::from_bits(0x3fec19a000000000),
    f64::from_bits(0x3feb7f8000000000),
    f64::from_bits(0x3feae8a000000000),
    f64::from_bits(0x3fea550000000000),
    f64::from_bits(0x3fe9c4a000000000),
    f64::from_bits(0x3fe9374000000000),
    f64::from_bits(0x3fe8ace000000000),
    f64::from_bits(0x3fe8258000000000),
    f64::from_bits(0x3fe7a12000000000),
    f64::from_bits(0x3fe71f8000000000),
    f64::from_bits(0x3fe6a0a000000000),
    f64::from_bits(0x3fe6248000000000),
    f64::from_bits(0x3fe5ab0000000000),
    f64::from_bits(0x3fe5342000000000),
    f64::from_bits(0x3fe4bfe000000000),
    f64::from_bits(0x3fe44e0000000000),
    f64::from_bits(0x3fe3dea000000000),
    f64::from_bits(0x3fe371a000000000),
    f64::from_bits(0x3fe3070000000000),
    f64::from_bits(0x3fe29ea000000000),
    f64::from_bits(0x3fe2388000000000),
    f64::from_bits(0x3fe1d48000000000),
    f64::from_bits(0x3fe172c000000000),
    f64::from_bits(0x3fe1130000000000),
    f64::from_bits(0x3fe0b56000000000),
    f64::from_bits(0x3fe059c000000000),
    f64::from_bits(0x3fe0000000000000),
];

const LOG_R2: [f64; 33] = [
    f64::from_bits(0x3ff0000000000000),
    f64::from_bits(0x3feffa7000000000),
    f64::from_bits(0x3feff4f000000000),
    f64::from_bits(0x3fefef6000000000),
    f64::from_bits(0x3fefe9e000000000),
    f64::from_bits(0x3fefe45000000000),
    f64::from_bits(0x3fefded000000000),
    f64::from_bits(0x3fefd94000000000),
    f64::from_bits(0x3fefd3c000000000),
    f64::from_bits(0x3fefce4000000000),
    f64::from_bits(0x3fefc8c000000000),
    f64::from_bits(0x3fefc34000000000),
    f64::from_bits(0x3fefbdc000000000),
    f64::from_bits(0x3fefb84000000000),
    f64::from_bits(0x3fefb2c000000000),
    f64::from_bits(0x3fefad4000000000),
    f64::from_bits(0x3fefa7c000000000),
    f64::from_bits(0x3fefa24000000000),
    f64::from_bits(0x3fef9cd000000000),
    f64::from_bits(0x3fef975000000000),
    f64::from_bits(0x3fef91e000000000),
    f64::from_bits(0x3fef8c6000000000),
    f64::from_bits(0x3fef86f000000000),
    f64::from_bits(0x3fef817000000000),
    f64::from_bits(0x3fef7c0000000000),
    f64::from_bits(0x3fef769000000000),
    f64::from_bits(0x3fef711000000000),
    f64::from_bits(0x3fef6ba000000000),
    f64::from_bits(0x3fef663000000000),
    f64::from_bits(0x3fef60c000000000),
    f64::from_bits(0x3fef5b5000000000),
    f64::from_bits(0x3fef55e000000000),
    f64::from_bits(0x3fef507000000000),
];

const LOG_L1: [[f64; 2]; 33] = [
    [
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
    ],
    [
        f64::from_bits(0x3da9f5e440f128db),
        f64::from_bits(0x3f962d07ab000000),
    ],
    [
        f64::from_bits(0xbda527d64b444fa3),
        f64::from_bits(0x3fa62f483d000000),
    ],
    [
        f64::from_bits(0x3d83aff57187d0cf),
        f64::from_bits(0x3fb0a26721400000),
    ],
    [
        f64::from_bits(0xbd64634c201e2b9c),
        f64::from_bits(0x3fb62e04bc000000),
    ],
    [
        f64::from_bits(0xbdbd46364a8017c7),
        f64::from_bits(0x3fbbb9db70800000),
    ],
    [
        f64::from_bits(0xbdb882b6acb3f696),
        f64::from_bits(0x3fc0a29f69c00000),
    ],
    [
        f64::from_bits(0x3da5a5833aeff542),
        f64::from_bits(0x3fc368507da00000),
    ],
    [
        f64::from_bits(0xbdb3876d32b0cbf5),
        f64::from_bits(0x3fc62e4116c00000),
    ],
    [
        f64::from_bits(0x3daf5712171380e6),
        f64::from_bits(0x3fc8f41d56800000),
    ],
    [
        f64::from_bits(0x3dbfc0b2e87a92c1),
        f64::from_bits(0x3fcbb98bc4c00000),
    ],
    [
        f64::from_bits(0x3db44c7ceb2f93f2),
        f64::from_bits(0x3fce7f71f0800000),
    ],
    [
        f64::from_bits(0x3daa147c39e44eba),
        f64::from_bits(0x3fd0a2bfe2c00000),
    ],
    [
        f64::from_bits(0x3da36d8fc46707d1),
        f64::from_bits(0x3fd205afe0300000),
    ],
    [
        f64::from_bits(0xbda0fd8155ea5850),
        f64::from_bits(0x3fd3685b58900000),
    ],
    [
        f64::from_bits(0x3da8954f1c1b010f),
        f64::from_bits(0x3fd4cb42e1900000),
    ],
    [
        f64::from_bits(0xbdb5d0bcd7fa4afa),
        f64::from_bits(0x3fd62e3e78c00000),
    ],
    [
        f64::from_bits(0xbdbb0a96458bf187),
        f64::from_bits(0x3fd7912364700000),
    ],
    [
        f64::from_bits(0x3dbc543eab5348b9),
        f64::from_bits(0x3fd8f42299600000),
    ],
    [
        f64::from_bits(0xbda15143e5c177e1),
        f64::from_bits(0x3fda5711c7e00000),
    ],
    [
        f64::from_bits(0x3d93be09bf52475c),
        f64::from_bits(0x3fdbb9c3ceb00000),
    ],
    [
        f64::from_bits(0xbd79b3b32e71e21d),
        f64::from_bits(0x3fdd1cd255b00000),
    ],
    [
        f64::from_bits(0xbd98f02175f93786),
        f64::from_bits(0x3fde7fb067100000),
    ],
    [
        f64::from_bits(0xbdbc5fb374b7ddcf),
        f64::from_bits(0x3fdfe2980ec00000),
    ],
    [
        f64::from_bits(0xbdb8e174c5571bbd),
        f64::from_bits(0x3fe0a2aef3500000),
    ],
    [
        f64::from_bits(0x3dbfa33ff819b3ec),
        f64::from_bits(0x3fe15420d4900000),
    ],
    [
        f64::from_bits(0x3d923d2634096ca6),
        f64::from_bits(0x3fe2058ca7900000),
    ],
    [
        f64::from_bits(0x3d9c8afc264146b2),
        f64::from_bits(0x3fe2b7156ff00000),
    ],
    [
        f64::from_bits(0x3dae21780abaa301),
        f64::from_bits(0x3fe3686c62000000),
    ],
    [
        f64::from_bits(0x3db3d67aee28cdc4),
        f64::from_bits(0x3fe419f01cd00000),
    ],
    [
        f64::from_bits(0x3dbccd8a77731be8),
        f64::from_bits(0x3fe4cb504d680000),
    ],
    [
        f64::from_bits(0x3da0cc7dc4dbbcfd),
        f64::from_bits(0x3fe57cb333b80000),
    ],
    [
        f64::from_bits(0x3db1cf79abc9e3b4),
        f64::from_bits(0x3fe62e42fef80000),
    ],
];

const LOG_L2: [[f64; 2]; 32] = [
    [
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
    ],
    [
        f64::from_bits(0x3db2ccace5b018a7),
        f64::from_bits(0x3f4641ef40000000),
    ],
    [
        f64::from_bits(0xbdb88a5cd275513a),
        f64::from_bits(0x3f5623d3f0000000),
    ],
    [
        f64::from_bits(0xbd90006a77b80a2d),
        f64::from_bits(0x3f60a45310000000),
    ],
    [
        f64::from_bits(0x3d881a0ebe451dd0),
        f64::from_bits(0x3f6627a998000000),
    ],
    [
        f64::from_bits(0x3da4297627f3b4ac),
        f64::from_bits(0x3f6bbc0150000000),
    ],
    [
        f64::from_bits(0x3dbafb8521676db1),
        f64::from_bits(0x3f70a0a0c8000000),
    ],
    [
        f64::from_bits(0x3db080b4c8bf43ca),
        f64::from_bits(0x3f736bc4a0000000),
    ],
    [
        f64::from_bits(0x3dabbcdc4ef244f5),
        f64::from_bits(0x3f762f5a48000000),
    ],
    [
        f64::from_bits(0xbdbeb4354215c794),
        f64::from_bits(0x3f78f36a44000000),
    ],
    [
        f64::from_bits(0xbdb020c23741371b),
        f64::from_bits(0x3f7bb7f4b8000000),
    ],
    [
        f64::from_bits(0xbdbff1289819f095),
        f64::from_bits(0x3f7e7cf9d4000000),
    ],
    [
        f64::from_bits(0x3db9b1383de1a3f4),
        f64::from_bits(0x3f80a13cde000000),
    ],
    [
        f64::from_bits(0x3d97fa1788e44213),
        f64::from_bits(0x3f82043a52000000),
    ],
    [
        f64::from_bits(0x3daf22c04fdaaa38),
        f64::from_bits(0x3f83677558000000),
    ],
    [
        f64::from_bits(0xbdbd9828c736de23),
        f64::from_bits(0x3f84caee08000000),
    ],
    [
        f64::from_bits(0xbdb41b9d79946440),
        f64::from_bits(0x3f862ea474000000),
    ],
    [
        f64::from_bits(0x3dbaced891ec8e07),
        f64::from_bits(0x3f879298b2000000),
    ],
    [
        f64::from_bits(0xbdbe4c5365e893ff),
        f64::from_bits(0x3f88f2be4e000000),
    ],
    [
        f64::from_bits(0x3dbe8d32fe7dc6f5),
        f64::from_bits(0x3f8a572dbe000000),
    ],
    [
        f64::from_bits(0xbdb3bf707f8ee0b7),
        f64::from_bits(0x3f8bb7cd50000000),
    ],
    [
        f64::from_bits(0xbdbc60b95a619b91),
        f64::from_bits(0x3f8d1cb84a000000),
    ],
    [
        f64::from_bits(0x3da64cbc2b83b45c),
        f64::from_bits(0x3f8e7dd224000000),
    ],
    [
        f64::from_bits(0xbda68543f75f32c6),
        f64::from_bits(0x3f8fe338fc000000),
    ],
    [
        f64::from_bits(0x3db421a5be17c2ec),
        f64::from_bits(0x3f90a266bb000000),
    ],
    [
        f64::from_bits(0xbdaf6b329a1da537),
        f64::from_bits(0x3f91534f84000000),
    ],
    [
        f64::from_bits(0xbda9a217f361c264),
        f64::from_bits(0x3f92065ff9000000),
    ],
    [
        f64::from_bits(0x3d73856e49eac8dd),
        f64::from_bits(0x3f92b78651000000),
    ],
    [
        f64::from_bits(0x3d8c6c3e72a945a1),
        f64::from_bits(0x3f9368cb54000000),
    ],
    [
        f64::from_bits(0xbdbdcfdc96d3a1c6),
        f64::from_bits(0x3f941a2f0d000000),
    ],
    [
        f64::from_bits(0x3d9a0da4813daf37),
        f64::from_bits(0x3f94cbb185000000),
    ],
    [
        f64::from_bits(0x3dbb4d3084ac1ad1),
        f64::from_bits(0x3f957d52c8000000),
    ],
];

const LOG_C: [f64; 4] = [
    f64::from_bits(0xbfdfffffffffffd3),
    f64::from_bits(0x3fd55555555543d5),
    f64::from_bits(0xbfd000002bb2d74e),
    f64::from_bits(0x3fc999a692c56e4e),
];

const LOG_H1: [[f64; 3]; 33] = [
    [
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
    ],
    [
        f64::from_bits(0x3a07052459b95dfc),
        f64::from_bits(0x3d4d79103c4a36a4),
        f64::from_bits(0x3f962d07ab330000),
    ],
    [
        f64::from_bits(0x39c6eaee7553baa0),
        f64::from_bits(0x3d460a6d2eec1732),
        f64::from_bits(0x3fa62f483cea8000),
    ],
    [
        f64::from_bits(0x3a0a79804e7ea455),
        f64::from_bits(0x3d4aff57187d0cec),
        f64::from_bits(0x3fb0a26721424000),
    ],
    [
        f64::from_bits(0x3a03d01ad8f4bdb1),
        f64::from_bits(0x3d3ce59eff0ea320),
        f64::from_bits(0x3fb62e04bbff4000),
    ],
    [
        f64::from_bits(0x39e28aff4e83414b),
        f64::from_bits(0x3d4ce4dabff41cb0),
        f64::from_bits(0x3fbbb9db70628000),
    ],
    [
        f64::from_bits(0x39f947079fb56ea3),
        f64::from_bits(0x3d4ea4a9a604b4f2),
        f64::from_bits(0x3fc0a29f69b3a000),
    ],
    [
        f64::from_bits(0x3a0071b8f347abdc),
        f64::from_bits(0x3d32c19d77faa108),
        f64::from_bits(0x3fc368507da56000),
    ],
    [
        f64::from_bits(0x3a060dc9eae5e8c7),
        f64::from_bits(0x3d4c4966a79a05a8),
        f64::from_bits(0x3fc62e4116b62000),
    ],
    [
        f64::from_bits(0x3a019de1930ac389),
        f64::from_bits(0x3d45c485c4e03992),
        f64::from_bits(0x3fc8f41d5687c000),
    ],
    [
        f64::from_bits(0x3a038f75cf02d618),
        f64::from_bits(0x3ce65d0f52581680),
        f64::from_bits(0x3fcbb98bc4cfe000),
    ],
    [
        f64::from_bits(0x39c66d0ece8cc055),
        f64::from_bits(0x3d28f9d65f27e3c0),
        f64::from_bits(0x3fce7f71f08a2000),
    ],
    [
        f64::from_bits(0x39e110c0a2d0d02a),
        f64::from_bits(0x3d247c39e44eb9a0),
        f64::from_bits(0x3fd0a2bfe2c34000),
    ],
    [
        f64::from_bits(0x3a0c2e1c7624937c),
        f64::from_bits(0x3d4b63f119c1f42a),
        f64::from_bits(0x3fd205afe0326000),
    ],
    [
        f64::from_bits(0x39d80a68afc2eec0),
        f64::from_bits(0x3cf3f550ad3d7d80),
        f64::from_bits(0x3fd3685b588de000),
    ],
    [
        f64::from_bits(0x3a02cb8203dc0bb6),
        f64::from_bits(0x3d254f1c1b010eb8),
        f64::from_bits(0x3fd4cb42e1931000),
    ],
    [
        f64::from_bits(0x3a0e234ee6f272d2),
        f64::from_bits(0x3d47a19402da833c),
        f64::from_bits(0x3fd62e3e78ba8000),
    ],
    [
        f64::from_bits(0x39f3d19c1b2c67cb),
        f64::from_bits(0x3d4ab4dd3a073c68),
        f64::from_bits(0x3fd7912364693000),
    ],
    [
        f64::from_bits(0x3a062e1c0308ab55),
        f64::from_bits(0x3d343eab5348b948),
        f64::from_bits(0x3fd8f42299671000),
    ],
    [
        f64::from_bits(0x3a09b5e4e30c5f0f),
        f64::from_bits(0x3d375e0d1f440f44),
        f64::from_bits(0x3fda5711c7ddd000),
    ],
    [
        f64::from_bits(0x39c32ba52f147867),
        f64::from_bits(0x3d47c137ea48eb80),
        f64::from_bits(0x3fdbb9c3ceb13000),
    ],
    [
        f64::from_bits(0x3a032f3b1a0f392d),
        f64::from_bits(0x3d4262668c70ef16),
        f64::from_bits(0x3fdd1cd255af9000),
    ],
    [
        f64::from_bits(0x39f078661f35a007),
        f64::from_bits(0x3d0fbd140d90f480),
        f64::from_bits(0x3fde7fb0670e7000),
    ],
    [
        f64::from_bits(0x39f0a4b6082eb02f),
        f64::from_bits(0x3d402645a4111858),
        f64::from_bits(0x3fdfe2980eb8e000),
    ],
    [
        f64::from_bits(0x3a08f3673c008d9d),
        f64::from_bits(0x3d3e8b3aa8e4433c),
        f64::from_bits(0x3fe0a2aef34ce000),
    ],
    [
        f64::from_bits(0x39fe4109a84a717e),
        f64::from_bits(0x3d419ffc0cd9f620),
        f64::from_bits(0x3fe15420d493f000),
    ],
    [
        f64::from_bits(0x39fb400798743957),
        f64::from_bits(0x3d2e931a04b652c0),
        f64::from_bits(0x3fe2058ca7909000),
    ],
    [
        f64::from_bits(0x39f3046edc6794a0),
        f64::from_bits(0x3d415f84c828d642),
        f64::from_bits(0x3fe2b7156ff0e000),
    ],
    [
        f64::from_bits(0x3a05cb5072b4b776),
        f64::from_bits(0x3d30bc055d518094),
        f64::from_bits(0x3fe3686c6201e000),
    ],
    [
        f64::from_bits(0x39e051804915de87),
        f64::from_bits(0x3d367aee28cdc41c),
        f64::from_bits(0x3fe419f01cd27800),
    ],
    [
        f64::from_bits(0x39d748ebb6888301),
        f64::from_bits(0x3d2b14eee637d048),
        f64::from_bits(0x3fe4cb504d6b9800),
    ],
    [
        f64::from_bits(0x39e55098cdd851ed),
        f64::from_bits(0x3d431f7136ef3f30),
        f64::from_bits(0x3fe57cb333b90800),
    ],
    [
        f64::from_bits(0x398f97b57a079a19),
        f64::from_bits(0x3d2ef35793c76730),
        f64::from_bits(0x3fe62e42fefa3800),
    ],
];

const LOG_H2: [[f64; 3]; 33] = [
    [
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
    ],
    [
        f64::from_bits(0x3a04826401258afc),
        f64::from_bits(0x3d2959cb60314e00),
        f64::from_bits(0x3f4641ef49600000),
    ],
    [
        f64::from_bits(0x3a09c71d4cdcb54b),
        f64::from_bits(0x3d4ad196c55762e2),
        f64::from_bits(0x3f5623d3e9d00000),
    ],
    [
        f64::from_bits(0x39f012c2a9a87c83),
        f64::from_bits(0x3d4ff2b108feba50),
        f64::from_bits(0x3f60a4530f780000),
    ],
    [
        f64::from_bits(0x3a0437bb5bf9fc34),
        f64::from_bits(0x3d0a0ebe451dd000),
        f64::from_bits(0x3f6627a998600000),
    ],
    [
        f64::from_bits(0x39f1c54e00c02c66),
        f64::from_bits(0x3d34bb13f9da5614),
        f64::from_bits(0x3f6bbc0151400000),
    ],
    [
        f64::from_bits(0x3a0b959e749a5ea5),
        f64::from_bits(0x3d4dc290b3b6d866),
        f64::from_bits(0x3f70a0a0c9ac0000),
    ],
    [
        f64::from_bits(0x3a09d74597e3c084),
        f64::from_bits(0x3ce69917e8794b00),
        f64::from_bits(0x3f736bc4a1080000),
    ],
    [
        f64::from_bits(0x39dee2ae1f36e35a),
        f64::from_bits(0x3d3e6e2779227a9c),
        f64::from_bits(0x3f762f5a48dc0000),
    ],
    [
        f64::from_bits(0x3a06236c8dcf0cfd),
        f64::from_bits(0x3d27957bd470d7f0),
        f64::from_bits(0x3f78f36a42140000),
    ],
    [
        f64::from_bits(0x3a0a291a656020e5),
        f64::from_bits(0x3d3f3dc8bec8e520),
        f64::from_bits(0x3f7bb7f4b6fc0000),
    ],
    [
        f64::from_bits(0x3a0a34920afb6399),
        f64::from_bits(0x3d2daecfcc1ed628),
        f64::from_bits(0x3f7e7cf9d2000000),
    ],
    [
        f64::from_bits(0x3a07b470c6ab05c6),
        f64::from_bits(0x3d489c1ef0d1fa32),
        f64::from_bits(0x3f80a13cdecc0000),
    ],
    [
        f64::from_bits(0x3a02c071f819920f),
        f64::from_bits(0x3d4f42f11c88426e),
        f64::from_bits(0x3f82043a522e0000),
    ],
    [
        f64::from_bits(0x3a03d45aed9f22a0),
        f64::from_bits(0x3d316027ed551c38),
        f64::from_bits(0x3f836775587c0000),
    ],
    [
        f64::from_bits(0x3a0c4ec608241822),
        f64::from_bits(0x3d43eb9c6490ee9e),
        f64::from_bits(0x3f84caee07120000),
    ],
    [
        f64::from_bits(0x3a05c705f1ca65ef),
        f64::from_bits(0x3d42314335cde020),
        f64::from_bits(0x3f862ea4735e0000),
    ],
    [
        f64::from_bits(0x39fc6e86fe6c904c),
        f64::from_bits(0x3d2db123d91c0e60),
        f64::from_bits(0x3f879298b2d60000),
    ],
    [
        f64::from_bits(0x3a0cf83c85a09c0d),
        f64::from_bits(0x3d49d64d0bb600a6),
        f64::from_bits(0x3f88f2be4d0c0000),
    ],
    [
        f64::from_bits(0x39f8a899399986f8),
        f64::from_bits(0x3d2a65fcfb8de998),
        f64::from_bits(0x3f8a572dbef40000),
    ],
    [
        f64::from_bits(0x39a75dd2d25f587c),
        f64::from_bits(0x3ce1f00e23e91000),
        f64::from_bits(0x3f8bb7cd4f620000),
    ],
    [
        f64::from_bits(0x3a066dc5826a7f9d),
        f64::from_bits(0x3d3f46a59e646ed8),
        f64::from_bits(0x3f8d1cb8491c0000),
    ],
    [
        f64::from_bits(0x3a00bdd376e8742a),
        f64::from_bits(0x3d432f0ae0ed1704),
        f64::from_bits(0x3f8e7dd224580000),
    ],
    [
        f64::from_bits(0x3a03bcc642260db7),
        f64::from_bits(0x3d4eaf0228334e7e),
        f64::from_bits(0x3f8fe338fba40000),
    ],
    [
        f64::from_bits(0x3a0fe76520a42aec),
        f64::from_bits(0x3d40d2df0be17626),
        f64::from_bits(0x3f90a266bb500000),
    ],
    [
        f64::from_bits(0x39e3f0eff67b7e50),
        f64::from_bits(0x3d24cd65e25ac960),
        f64::from_bits(0x3f91534f83c10000),
    ],
    [
        f64::from_bits(0x39c824751a1e6c5b),
        f64::from_bits(0x3d477a03278f66f6),
        f64::from_bits(0x3f92065ff8cc0000),
    ],
    [
        f64::from_bits(0x3a0506cfc4639acc),
        f64::from_bits(0x3d4c2b724f5646e6),
        f64::from_bits(0x3f92b78651040000),
    ],
    [
        f64::from_bits(0x39ca7ac1fdfa37e0),
        f64::from_bits(0x3d2b0f9caa516828),
        f64::from_bits(0x3f9368cb540e0000),
    ],
    [
        f64::from_bits(0x39fb61c029607481),
        f64::from_bits(0x3d4811b4962f1cf2),
        f64::from_bits(0x3f941a2f0c880000),
    ],
    [
        f64::from_bits(0x39e2879f319383ad),
        f64::from_bits(0x3d0b49027b5e6de0),
        f64::from_bits(0x3f94cbb1851a0000),
    ],
    [
        f64::from_bits(0x39ef27f4e668c317),
        f64::from_bits(0x3d2a61095835a298),
        f64::from_bits(0x3f957d52c86d0000),
    ],
    [
        f64::from_bits(0x39fc301f232c0e74),
        f64::from_bits(0x3d4b6fc11defa4a8),
        f64::from_bits(0x3f962f12e1320000),
    ],
];

const LOG_C_ACC: [[f64; 2]; 9] = [
    [
        f64::from_bits(0x3ff0000000000000),
        f64::from_bits(0x389a193d7f59d80a),
    ],
    [
        f64::from_bits(0xbfe0000000000000),
        f64::from_bits(0x39c8c7d7a8733406),
    ],
    [
        f64::from_bits(0x3fd5555555555555),
        f64::from_bits(0x3c755555554f571d),
    ],
    [
        f64::from_bits(0xbfd0000000000000),
        f64::from_bits(0xbb6e516b5d7b8c15),
    ],
    [
        f64::from_bits(0x3fc999999999999a),
        f64::from_bits(0xbc697f2898534175),
    ],
    [
        f64::from_bits(0xbfc55555555554b5),
        f64::from_bits(0xbc43834d62d64ec4),
    ],
    [
        f64::from_bits(0x3fc249249248dbdc),
        f64::from_bits(0xbc58aa032979ebed),
    ],
    [
        f64::from_bits(0xbfc000004e71581b),
        f64::from_bits(0x3c62e7f17d1c0e63),
    ],
    [
        f64::from_bits(0x3fbc71e5ec7051f6),
        f64::from_bits(0x3c5217ec3dcb2f03),
    ],
];

const STPI: [[f64; 2]; 65] = [
    [
        f64::from_bits(0x0000000000000000),
        f64::from_bits(0x0000000000000000),
    ],
    [
        f64::from_bits(0x3bfc14eff99a3ff1),
        f64::from_bits(0x3f7fff2d746c8895),
    ],
    [
        f64::from_bits(0xbc18c4d4c1bbe38b),
        f64::from_bits(0x3f8ffcb5e52d1f36),
    ],
    [
        f64::from_bits(0xbc208ef2408930eb),
        f64::from_bits(0x3f97fa7329846feb),
    ],
    [
        f64::from_bits(0xbc314daa07929354),
        f64::from_bits(0x3f9ff2d8cc5320c7),
    ],
    [
        f64::from_bits(0x3c3d845cf264d016),
        f64::from_bits(0x3fa3f3289bb44643),
    ],
    [
        f64::from_bits(0xbc343aa63f69acea),
        f64::from_bits(0x3fa7e9d144d37f33),
    ],
    [
        f64::from_bits(0xbc4bc90382ed68a4),
        f64::from_bits(0x3fabdcc9ea69fc93),
    ],
    [
        f64::from_bits(0x3c30fbc215a3c756),
        f64::from_bits(0x3fafcb76a6ecccab),
    ],
    [
        f64::from_bits(0x3c572b75e84ab5e2),
        f64::from_bits(0x3fb1da9e1f36c497),
    ],
    [
        f64::from_bits(0xbc420d100fccf991),
        f64::from_bits(0x3fb3ccc01b453709),
    ],
    [
        f64::from_bits(0xbc0f7aac846eccfd),
        f64::from_bits(0x3fb5bbd477204be0),
    ],
    [
        f64::from_bits(0xbc417799578a6651),
        f64::from_bits(0x3fb7a78edace5e27),
    ],
    [
        f64::from_bits(0x3c50c85deb5bb812),
        f64::from_bits(0x3fb98fa372a35c37),
    ],
    [
        f64::from_bits(0xbc367d2eb81bbf36),
        f64::from_bits(0x3fbb73c6faf2275c),
    ],
    [
        f64::from_bits(0xbc014b2141507a9d),
        f64::from_bits(0x3fbd53aecba7bf00),
    ],
    [
        f64::from_bits(0xbc58939cffeb036c),
        f64::from_bits(0x3fbf2f10e3ce6d42),
    ],
    [
        f64::from_bits(0xbc62f3fbb178d1c5),
        f64::from_bits(0x3fc082d1fa7b9738),
    ],
    [
        f64::from_bits(0x3c608479c62d3d77),
        f64::from_bits(0x3fc16b8fb743c879),
    ],
    [
        f64::from_bits(0xbc6894149dc3b5f7),
        f64::from_bits(0x3fc2519dc47527b3),
    ],
    [
        f64::from_bits(0xbc344dad213ab344),
        f64::from_bits(0x3fc334d8a850758d),
    ],
    [
        f64::from_bits(0xbc52d415416bae28),
        f64::from_bits(0x3fc4151d589a490f),
    ],
    [
        f64::from_bits(0xbc62e0d0b51ed237),
        f64::from_bits(0x3fc4f24940025067),
    ],
    [
        f64::from_bits(0x3c6e8045a3cf3213),
        f64::from_bits(0x3fc5cc3a43788a30),
    ],
    [
        f64::from_bits(0x3c6be4e50e1bf91f),
        f64::from_bits(0x3fc6a2cec76fa4b0),
    ],
    [
        f64::from_bits(0x3c1b1e18c1f7f635),
        f64::from_bits(0x3fc775e5b50bb365),
    ],
    [
        f64::from_bits(0xbc010946c1f6f484),
        f64::from_bits(0x3fc8455e7f3c6e5a),
    ],
    [
        f64::from_bits(0x3c4291a88889a4e6),
        f64::from_bits(0x3fc9111927c231cf),
    ],
    [
        f64::from_bits(0xbc6bedd6f9a25da4),
        f64::from_bits(0x3fc9d8f6441cf80b),
    ],
    [
        f64::from_bits(0xbc6ce108006670c7),
        f64::from_bits(0x3fca9cd702648a97),
    ],
    [
        f64::from_bits(0x3c24c65624119572),
        f64::from_bits(0x3fcb5c9d2e092baa),
    ],
    [
        f64::from_bits(0xbc6a26e2a2682111),
        f64::from_bits(0x3fcc182b347bfc21),
    ],
    [
        f64::from_bits(0x3c4fce159c2bb59b),
        f64::from_bits(0x3fcccf6429be6621),
    ],
    [
        f64::from_bits(0xbc559b1cffa69603),
        f64::from_bits(0x3fcd822bccd7d86e),
    ],
    [
        f64::from_bits(0x3c6677083288397a),
        f64::from_bits(0x3fce30668c31224e),
    ],
    [
        f64::from_bits(0x3c69a49696faa0ec),
        f64::from_bits(0x3fced9f989d4c415),
    ],
    [
        f64::from_bits(0x3c5ca323e77a3345),
        f64::from_bits(0x3fcf7eca9f938c6f),
    ],
    [
        f64::from_bits(0xbc6c702625d3863b),
        f64::from_bits(0x3fd00f6031866f76),
    ],
    [
        f64::from_bits(0x3c7180cf0e52237d),
        f64::from_bits(0x3fd05ce114cd024a),
    ],
    [
        f64::from_bits(0xbc74be56fec860b9),
        f64::from_bits(0x3fd0a7dc060df5ee),
    ],
    [
        f64::from_bits(0x3c5b5d970e5d9d07),
        f64::from_bits(0x3fd0f045755560d9),
    ],
    [
        f64::from_bits(0x3c54c32e06c67499),
        f64::from_bits(0x3fd1361238136929),
    ],
    [
        f64::from_bits(0xbc7b512d49aedaa1),
        f64::from_bits(0x3fd179378ad51274),
    ],
    [
        f64::from_bits(0xbc5161478130996d),
        f64::from_bits(0x3fd1b9ab12ed2518),
    ],
    [
        f64::from_bits(0xbc425feb091e921f),
        f64::from_bits(0x3fd1f762e00ced83),
    ],
    [
        f64::from_bits(0x3c73750bc95dae67),
        f64::from_bits(0x3fd232556dcc945f),
    ],
    [
        f64::from_bits(0xbc7257966a1044c5),
        f64::from_bits(0x3fd26a79a522d332),
    ],
    [
        f64::from_bits(0xbc6ac6af78c05e44),
        f64::from_bits(0x3fd29fc6ddcbcb72),
    ],
    [
        f64::from_bits(0x3c771dbd64ba4f95),
        f64::from_bits(0x3fd2d234df9ec8c9),
    ],
    [
        f64::from_bits(0x3c6020107d2c17b0),
        f64::from_bits(0x3fd301bbe3d2b9c7),
    ],
    [
        f64::from_bits(0xbc79d4016f0b15c4),
        f64::from_bits(0x3fd32e5496312cfc),
    ],
    [
        f64::from_bits(0x3c7f557b51b587cc),
        f64::from_bits(0x3fd357f81637a329),
    ],
    [
        f64::from_bits(0x3c6c88cee9bad9f9),
        f64::from_bits(0x3fd37e9ff82709ec),
    ],
    [
        f64::from_bits(0x3c7ebde6bb284e87),
        f64::from_bits(0x3fd3a246460134f7),
    ],
    [
        f64::from_bits(0xbc78b1d8c40ffea3),
        f64::from_bits(0x3fd3c2e580742eda),
    ],
    [
        f64::from_bits(0x3c7de48797b477f2),
        f64::from_bits(0x3fd3e0789fb33cf7),
    ],
    [
        f64::from_bits(0xbc78bc6105a80fa5),
        f64::from_bits(0x3fd3fafb143d754b),
    ],
    [
        f64::from_bits(0xbc79f4cc680744f3),
        f64::from_bits(0x3fd41268c791c743),
    ],
    [
        f64::from_bits(0x3c32ed295e9d0ef2),
        f64::from_bits(0x3fd426be1cd05c06),
    ],
    [
        f64::from_bits(0xbc34a98b72ed3789),
        f64::from_bits(0x3fd437f7f1493531),
    ],
    [
        f64::from_bits(0xbc75c080cdd72ddf),
        f64::from_bits(0x3fd446139cf7f413),
    ],
    [
        f64::from_bits(0xbc6b00c622ae015e),
        f64::from_bits(0x3fd4510ef2ecb654),
    ],
    [
        f64::from_bits(0x3c78dd5ec4960646),
        f64::from_bits(0x3fd458e841a1f7da),
    ],
    [
        f64::from_bits(0x3c7e1f89d1adcbc6),
        f64::from_bits(0x3fd45d9e533f6cac),
    ],
    [
        f64::from_bits(0xbc76b01ec5417056),
        f64::from_bits(0x3fd45f306dc9c883),
    ],
];

const SINPID_C_NEAR: [f64; 2] = [
    f64::from_bits(0xbffa51a6625307d3),
    f64::from_bits(0xbc816cc8f2044a4a),
];

const SINPID_CL_NEAR: [f64; 3] = [
    f64::from_bits(0x3fe9f9cb402bc42a),
    f64::from_bits(0xbfc86a8e46ddf78d),
    f64::from_bits(0x3f9ac644e7aa33e6),
];

const SINPID_C_COEFF: [f64; 4] = [
    f64::from_bits(0xbf33bd3cc9be45de),
    f64::from_bits(0x3e503c1f081b5ac4),
    f64::from_bits(0xbd555d3c7e3bd8bf),
    f64::from_bits(0x3c4e1f4826790653),
];

const SINPID_S_COEFF: [f64; 4] = [
    f64::from_bits(0x3f9921fb54442d18),
    f64::from_bits(0xbec4abbce625be53),
    f64::from_bits(0x3dd466bc67748efc),
    f64::from_bits(0xbcd32d26e446373a),
];

const SINPID_ACC_C: [[f64; 2]; 5] = [
    [
        f64::from_bits(0xbf33bd3cc9be45de),
        f64::from_bits(0xbbd692b71366c792),
    ],
    [
        f64::from_bits(0x3e503c1f081b5ac4),
        f64::from_bits(0xbaf32b342c5da62e),
    ],
    [
        f64::from_bits(0xbd555d3c7e3cbffa),
        f64::from_bits(0x39d17fcb68af9aae),
    ],
    [
        f64::from_bits(0x3c4e1f506890e556),
        f64::from_bits(0xb8d56c9658a7cfd5),
    ],
    [
        f64::from_bits(0xbb3a6d193ca649a0),
        f64::from_bits(0x37961669b56f0275),
    ],
];

const SINPID_ACC_S: [[f64; 2]; 6] = [
    [
        f64::from_bits(0x3f9921fb54442d18),
        f64::from_bits(0x3c31a62633145c07),
    ],
    [
        f64::from_bits(0xbec4abbce625be53),
        f64::from_bits(0x3b605511c6847960),
    ],
    [
        f64::from_bits(0x3dd466bc6775aae2),
        f64::from_bits(0xba66dc0d14c26b21),
    ],
    [
        f64::from_bits(0xbcd32d2cce62bd86),
        f64::from_bits(0x3972ba6dbc4b37c0),
    ],
    [
        f64::from_bits(0x3bc50783487e6b5c),
        f64::from_bits(0xb85b77efed0f9c1f),
    ],
    [
        f64::from_bits(0xbaae306ec8cf7c02),
        f64::from_bits(0x3745d2601f85289a),
    ],
];

const LGAMMA_ASYM_ACC_C1: [[f64; 2]; 8] = [
    [
        f64::from_bits(0x3fdacfe390c97d69),
        f64::from_bits(0x3c73494bc9001766),
    ],
    [
        f64::from_bits(0x3fb5555555555555),
        f64::from_bits(0x3c555555554133c7),
    ],
    [
        f64::from_bits(0xbf66c16c16c16c17),
        f64::from_bits(0x3bff4c03199d8517),
    ],
    [
        f64::from_bits(0x3f4a01a01a01a016),
        f64::from_bits(0x3bdfb315e77b4883),
    ],
    [
        f64::from_bits(0xbf4381381380c1b0),
        f64::from_bits(0xbbec4cc418316ed1),
    ],
    [
        f64::from_bits(0x3f4b951e22dd8dfc),
        f64::from_bits(0x3bee8da392824ecf),
    ],
    [
        f64::from_bits(0xbf5f6a875bb1ab7b),
        f64::from_bits(0x3bf27c5fcbab6b5d),
    ],
    [
        f64::from_bits(0x3f7a0a6926f49920),
        f64::from_bits(0xbc01f355cbf82229),
    ],
];

const LGAMMA_ASYM_ACC_C2: [[f64; 2]; 12] = [
    [
        f64::from_bits(0x3fdacfe390c97d69),
        f64::from_bits(0x3c73494bc9007f28),
    ],
    [
        f64::from_bits(0x3fb5555555555555),
        f64::from_bits(0x3c55555554ad7655),
    ],
    [
        f64::from_bits(0xbf66c16c16c16c17),
        f64::from_bits(0x3bff4b73ea546bd4),
    ],
    [
        f64::from_bits(0x3f4a01a01a01a01a),
        f64::from_bits(0xbbe464e31b1a609a),
    ],
    [
        f64::from_bits(0xbf43813813813692),
        f64::from_bits(0x3bd63676fd3c851e),
    ],
    [
        f64::from_bits(0x3f4b951e2b143b5e),
        f64::from_bits(0x3becfe48021143d1),
    ],
    [
        f64::from_bits(0xbf5f6ab0d459cadb),
        f64::from_bits(0x3bf770df5ee0beea),
    ],
    [
        f64::from_bits(0x3f7a41a211f4e098),
        f64::from_bits(0xbc0e00b3c1619519),
    ],
    [
        f64::from_bits(0xbf9e41f97ad4634d),
        f64::from_bits(0x3c305ffdbc72560f),
    ],
    [
        f64::from_bits(0x3fc6f15ef2b47719),
        f64::from_bits(0xbc65665cba69dbaa),
    ],
    [
        f64::from_bits(0xbff5762c9f49fe25),
        f64::from_bits(0x3c9afeff5294ad13),
    ],
    [
        f64::from_bits(0x4022ea102098f818),
        f64::from_bits(0x3cc609db97f1bc89),
    ],
];

const LGAMMA_ASYM_ACC_C3: [[f64; 2]; 28] = [
    [
        f64::from_bits(0x3fdacfe390c97d69),
        f64::from_bits(0x3c73494bce9b5c50),
    ],
    [
        f64::from_bits(0x3fb5555555555555),
        f64::from_bits(0x3c555551a0d18a1d),
    ],
    [
        f64::from_bits(0xbf66c16c16c16c17),
        f64::from_bits(0x3c007171d4a61bb9),
    ],
    [
        f64::from_bits(0x3f4a01a01a01a008),
        f64::from_bits(0xbbe49af39145d6e0),
    ],
    [
        f64::from_bits(0xbf438138138124cc),
        f64::from_bits(0xbbdb8f71068f5292),
    ],
    [
        f64::from_bits(0x3f4b951e2b09b070),
        f64::from_bits(0xbbd3c8a4c099ae72),
    ],
    [
        f64::from_bits(0xbf5f6ab0d4de4a5c),
        f64::from_bits(0xbbefbb4ed542d17e),
    ],
    [
        f64::from_bits(0x3f7a41a384fafd09),
        f64::from_bits(0x3bf1510cc5dec148),
    ],
    [
        f64::from_bits(0xbf9e4277d1a2b2e4),
        f64::from_bits(0x3c233cdd66b223e7),
    ],
    [
        f64::from_bits(0x3fc6fdf731005805),
        f64::from_bits(0xbc6a5a3066922140),
    ],
    [
        f64::from_bits(0xbff641de6a9b2f1c),
        f64::from_bits(0x3c727cfac68728a3),
    ],
    [
        f64::from_bits(0x402aa463d5553a9f),
        f64::from_bits(0x3cbaeba9a5b525f3),
    ],
    [
        f64::from_bits(0xc06312d3aa56b6a2),
        f64::from_bits(0xbcf2b2c4e643e5b9),
    ],
    [
        f64::from_bits(0x409f3e8a4cd3d268),
        f64::from_bits(0xbd2ec8b086b0215e),
    ],
    [
        f64::from_bits(0xc0dba9362f228307),
        f64::from_bits(0x3d7b8ecae7dba56f),
    ],
    [
        f64::from_bits(0x4118ed4df00421e4),
        f64::from_bits(0x3dbc620cd11f44a5),
    ],
    [
        f64::from_bits(0xc155af889993596d),
        f64::from_bits(0xbdf23a78dbbfe515),
    ],
    [
        f64::from_bits(0x4191789aff6ddc93),
        f64::from_bits(0x3e2fe7e161d8bdb1),
    ],
    [
        f64::from_bits(0xc1c94353def0e5f8),
        f64::from_bits(0xbe664a798491c211),
    ],
    [
        f64::from_bits(0x41fffb9253861b9c),
        f64::from_bits(0xbe8469d6acfea485),
    ],
    [
        f64::from_bits(0xc23159c1244c9ee4),
        f64::from_bits(0xbecf7491d4abcc7b),
    ],
    [
        f64::from_bits(0x425f991ba9776c6a),
        f64::from_bits(0x3eccf2c83cfe2e55),
    ],
    [
        f64::from_bits(0xc28795f45ae3fa72),
        f64::from_bits(0xbf2a777148d59de4),
    ],
    [
        f64::from_bits(0x42ac05845a910b48),
        f64::from_bits(0x3f4433e08525d4bd),
    ],
    [
        f64::from_bits(0xc2c96c399c145d71),
        f64::from_bits(0x3f4d495c9370dc8a),
    ],
    [
        f64::from_bits(0x42e083bb9c79529b),
        f64::from_bits(0xbf71640e2de8cf9c),
    ],
    [
        f64::from_bits(0xc2eb53232e707f4b),
        f64::from_bits(0x3f76d239acd71b5c),
    ],
    [
        f64::from_bits(0x42e59bad61bd81d5),
        f64::from_bits(0x3f80e3f2ea42a0e0),
    ],
];
const LGAMMA_ACC_X0_0: [f64; 3] = [
    f64::from_bits(0x4003a7fc9600f86c),
    f64::from_bits(0x3c855f64f98af8d0),
    f64::from_bits(0x391c4b0cd201366a),
];

const LGAMMA_ACC_C_0: [[f64; 2]; 26] = [
    [
        f64::from_bits(0x3ff83fe966af535f),
        f64::from_bits(0xbc8775909a36a68c),
    ],
    [
        f64::from_bits(0x3fe36eebb002f55d),
        f64::from_bits(0xbc88d4b2124a39f8),
    ],
    [
        f64::from_bits(0x3f9694a6058a7858),
        f64::from_bits(0xbc21d8c8b9b4d80e),
    ],
    [
        f64::from_bits(0x3f91718d7ca09e5b),
        f64::from_bits(0x3c383195b07fd250),
    ],
    [
        f64::from_bits(0x3f57339fe04b2764),
        f64::from_bits(0xbbf48648f4a4bf9e),
    ],
    [
        f64::from_bits(0x3f48d32f682aa0bd),
        f64::from_bits(0xbbd90953dadfba01),
    ],
    [
        f64::from_bits(0x3f1809f04ee6e0fa),
        f64::from_bits(0xbbb6100a3d177e7c),
    ],
    [
        f64::from_bits(0x3f048eaa81657361),
        f64::from_bits(0x3b942b86f83f623d),
    ],
    [
        f64::from_bits(0x3ed9297adb2def5a),
        f64::from_bits(0xbb60c295492288fd),
    ],
    [
        f64::from_bits(0x3ec286fb8cbaebb5),
        f64::from_bits(0xbb2ebf1f6a584b62),
    ],
    [
        f64::from_bits(0x3e9a92e0a5de4bc9),
        f64::from_bits(0xbb3688db240a9a82),
    ],
    [
        f64::from_bits(0x3e81a9d4d8c6284e),
        f64::from_bits(0x3b108624c7186639),
    ],
    [
        f64::from_bits(0x3e5c4cd2594e91c9),
        f64::from_bits(0x3aea225e59cfef5d),
    ],
    [
        f64::from_bits(0x3e418737ec8e68aa),
        f64::from_bits(0x3ad6421ab9e4c553),
    ],
    [
        f64::from_bits(0x3e1e6028795be5c1),
        f64::from_bits(0xbabb4bd2f36bb0dd),
    ],
    [
        f64::from_bits(0x3e01eacaeb800afd),
        f64::from_bits(0xba8cedc75cbd8cc4),
    ],
    [
        f64::from_bits(0x3de06bce9e1f6b90),
        f64::from_bits(0xba80fb524fda64eb),
    ],
    [
        f64::from_bits(0x3dc2bb110a516b79),
        f64::from_bits(0x3a52d7a224941f92),
    ],
    [
        f64::from_bits(0x3da1e0024589b848),
        f64::from_bits(0x3a31a84fa416703d),
    ],
    [
        f64::from_bits(0x3d83ec6ebb4ba730),
        f64::from_bits(0xba1eb9cd73c631d9),
    ],
    [
        f64::from_bits(0x3d63936af34a0150),
        f64::from_bits(0xba083b442df73674),
    ],
    [
        f64::from_bits(0x3d457d4ece262198),
        f64::from_bits(0x39e20b160178fca3),
    ],
    [
        f64::from_bits(0x3d25b530d802ffff),
        f64::from_bits(0x3978bb2ee582dbd1),
    ],
    [
        f64::from_bits(0x3d07c323d20053d0),
        f64::from_bits(0xb9832d5dd28cac13),
    ],
    [
        f64::from_bits(0x3ce6131339e2b76a),
        f64::from_bits(0x3974422b26549563),
    ],
    [
        f64::from_bits(0x3cbb4f6aaa0d8886),
        f64::from_bits(0x395796e66ee63a34),
    ],
];

const LGAMMA_ACC_X0_1: [f64; 3] = [
    f64::from_bits(0x4005fb410a1bd901),
    f64::from_bits(0xbc9a19a96d2e6f85),
    f64::from_bits(0xb93140b4ff4b7d60),
];

const LGAMMA_ACC_C_1: [[f64; 2]; 27] = [
    [
        f64::from_bits(0xbffea12da904b18c),
        f64::from_bits(0xbc9220130f99b2cb),
    ],
    [
        f64::from_bits(0x3fe3267f3c265a52),
        f64::from_bits(0xbc81c630ff19db83),
    ],
    [
        f64::from_bits(0xbfb4185ac30c8bf2),
        f64::from_bits(0x3c4f161263693e12),
    ],
    [
        f64::from_bits(0x3f8f504accc9f19b),
        f64::from_bits(0xbc1eacc0226c7208),
    ],
    [
        f64::from_bits(0xbf68588458207eac),
        f64::from_bits(0x3c04b51668f3dff4),
    ],
    [
        f64::from_bits(0x3f44373f7cc709b3),
        f64::from_bits(0xbbd2474e20e777ae),
    ],
    [
        f64::from_bits(0xbf212239bdd6c013),
        f64::from_bits(0x3ba46df69aa3032c),
    ],
    [
        f64::from_bits(0x3efdba65e27421c4),
        f64::from_bits(0x3b836bfe12004625),
    ],
    [
        f64::from_bits(0xbeda2d2504d7e987),
        f64::from_bits(0x3b4a6b6dfe9fa6f4),
    ],
    [
        f64::from_bits(0x3eb7581739ee6087),
        f64::from_bits(0xbb223100aff1ab78),
    ],
    [
        f64::from_bits(0xbe9506c65fad6188),
        f64::from_bits(0xbb217d1749eb738f),
    ],
    [
        f64::from_bits(0x3e7318ef724f7814),
        f64::from_bits(0xbad6517edf27abc2),
    ],
    [
        f64::from_bits(0xbe517767260d9825),
        f64::from_bits(0xbad03e683dfe1d87),
    ],
    [
        f64::from_bits(0x3e3011e34454ade3),
        f64::from_bits(0xbad2a844fd6e7622),
    ],
    [
        f64::from_bits(0xbe0db8b9e6b1e0e1),
        f64::from_bits(0xbaaf73094d7d4e40),
    ],
    [
        f64::from_bits(0x3deb9bab1fca3210),
        f64::from_bits(0xba789b7151db3448),
    ],
    [
        f64::from_bits(0xbdc9bed46f4d0aa2),
        f64::from_bits(0xba51af984f4b8e40),
    ],
    [
        f64::from_bits(0x3da81780d4242119),
        f64::from_bits(0xba47576147eed5dc),
    ],
    [
        f64::from_bits(0xbd869d3ce15a3f87),
        f64::from_bits(0x3a29823fc1b6be9e),
    ],
    [
        f64::from_bits(0x3d65494a34ce8470),
        f64::from_bits(0xba00ec96907a37b4),
    ],
    [
        f64::from_bits(0xbd44160eccd5982e),
        f64::from_bits(0x39d2b4a92d7b0c76),
    ],
    [
        f64::from_bits(0x3d22fe0a640fb1a0),
        f64::from_bits(0xb95b1f41cf9689ed),
    ],
    [
        f64::from_bits(0xbd01ff8659402df5),
        f64::from_bits(0x397511b4097b0fcc),
    ],
    [
        f64::from_bits(0x3ce1357a3e9e1cf1),
        f64::from_bits(0xb95bda0e37acdd97),
    ],
    [
        f64::from_bits(0xbcc0b145c5c5b9ab),
        f64::from_bits(0xb96b52b4002dbc0c),
    ],
    [
        f64::from_bits(0x3c9dc1888c7036ca),
        f64::from_bits(0x3922a69b093c9a8b),
    ],
    [
        f64::from_bits(0xbc7060b52b5d6f68),
        f64::from_bits(0x391d17de397a2b10),
    ],
];

const LGAMMA_ACC_X0_2: [f64; 3] = [
    f64::from_bits(0x4009260dbc9e59af),
    f64::from_bits(0x3caf717cd335a7b3),
    f64::from_bits(0x394d32a2a65bfd63),
];

const LGAMMA_ACC_C_2: [[f64; 2]; 20] = [
    [
        f64::from_bits(0x401f20a65f2fac55),
        f64::from_bits(0xbca1d258e4b0beb5),
    ],
    [
        f64::from_bits(0x3fd9d4d2977150ef),
        f64::from_bits(0x3c7a04089578ae88),
    ],
    [
        f64::from_bits(0x3f9c1137124d5c5b),
        f64::from_bits(0x3c2d6c922d2512d6),
    ],
    [
        f64::from_bits(0x3f6267203d776b0e),
        f64::from_bits(0xbc0aa6081d790e05),
    ],
    [
        f64::from_bits(0x3f299a6337da39dd),
        f64::from_bits(0x3bb49aeda9400147),
    ],
    [
        f64::from_bits(0x3ef293c3f78d3bdb),
        f64::from_bits(0x3b6ee59e48e8b181),
    ],
    [
        f64::from_bits(0x3ebbb97aa0b71e45),
        f64::from_bits(0xbb1592602ea73795),
    ],
    [
        f64::from_bits(0x3e851ea3345f5349),
        f64::from_bits(0x3b0ef9552f50f84e),
    ],
    [
        f64::from_bits(0x3e5057f65c64b210),
        f64::from_bits(0xbaf36ca6c4396475),
    ],
    [
        f64::from_bits(0x3e199c8650e3a4f9),
        f64::from_bits(0x3ab7d3b9d2266a10),
    ],
    [
        f64::from_bits(0x3de44520c3a4f841),
        f64::from_bits(0x3a68c284070d32df),
    ],
    [
        f64::from_bits(0x3db02d221961d8b5),
        f64::from_bits(0xba50d56a1938ace2),
    ],
    [
        f64::from_bits(0x3d79ffcd97899d22),
        f64::from_bits(0xba0d6040c8c03da8),
    ],
    [
        f64::from_bits(0x3d450494b6c20faa),
        f64::from_bits(0xb9e644c84757ab17),
    ],
    [
        f64::from_bits(0x3d1113fe94711cb9),
        f64::from_bits(0xb9bf7b0588af7045),
    ],
    [
        f64::from_bits(0x3cdbe09d37bab523),
        f64::from_bits(0x3950f8dbcb173a93),
    ],
    [
        f64::from_bits(0x3ca6d6b80656d502),
        f64::from_bits(0xb919e5c21b32eb00),
    ],
    [
        f64::from_bits(0x3c72ca70c7ef211e),
        f64::from_bits(0x3906511b75c60833),
    ],
    [
        f64::from_bits(0x3c3f9264b4df26e4),
        f64::from_bits(0xb86aefcabcc01fd2),
    ],
    [
        f64::from_bits(0x3c09262efb579250),
        f64::from_bits(0xb8a1196c645b5736),
    ],
];

const LGAMMA_ACC_X0_3: [f64; 3] = [
    f64::from_bits(0x400fa471547c2fe5),
    f64::from_bits(0x3c770d4561291237),
    f64::from_bits(0xb909e6fadbbc171a),
];

const LGAMMA_ACC_C_3: [[f64; 2]; 20] = [
    [
        f64::from_bits(0xc034b99d966c5647),
        f64::from_bits(0x3cd9cba2450b0003),
    ],
    [
        f64::from_bits(0x3fef76deae0436be),
        f64::from_bits(0xbc85af99a1af2d4b),
    ],
    [
        f64::from_bits(0xbfad25359d4b2f38),
        f64::from_bits(0x3c310c02bb27b3e8),
    ],
    [
        f64::from_bits(0x3f6e8f829f141aa5),
        f64::from_bits(0x3be4b3ff24054d7d),
    ],
    [
        f64::from_bits(0xbf3116f7806d26d3),
        f64::from_bits(0xbbba2efbebc9bc45),
    ],
    [
        f64::from_bits(0x3ef3e8f3ab9fc1f4),
        f64::from_bits(0x3b9e3a4940cdb9f6),
    ],
    [
        f64::from_bits(0xbeb7dbbe062ffd9e),
        f64::from_bits(0xbb52986222077472),
    ],
    [
        f64::from_bits(0x3e7d2f76de7bd027),
        f64::from_bits(0xbb16c0ccf4d8669c),
    ],
    [
        f64::from_bits(0xbe42225fe4f84932),
        f64::from_bits(0x3ada38a5afc8eeba),
    ],
    [
        f64::from_bits(0x3e06d12ae1936ba4),
        f64::from_bits(0xba9de21fabe64230),
    ],
    [
        f64::from_bits(0xbdccffc2a8f6492d),
        f64::from_bits(0x3a464881462031d6),
    ],
    [
        f64::from_bits(0x3d9294e1bdd838bc),
        f64::from_bits(0x3a3dcbb67f373f7d),
    ],
    [
        f64::from_bits(0xbd57fab625b2f2ff),
        f64::from_bits(0xb9be8527e1dc6098),
    ],
    [
        f64::from_bits(0x3d1f211abe57d760),
        f64::from_bits(0x399d65f8610deeee),
    ],
    [
        f64::from_bits(0xbce44f2a05f33c9d),
        f64::from_bits(0x3953fb27bfc54f94),
    ],
    [
        f64::from_bits(0x3caa9e4622460be6),
        f64::from_bits(0x3942f5912b9785c0),
    ],
    [
        f64::from_bits(0xbc7182776b50f5cf),
        f64::from_bits(0xb91d30d2e8d300a5),
    ],
    [
        f64::from_bits(0x3c3722b65d64ddc0),
        f64::from_bits(0xb8d0ee84a196d89c),
    ],
    [
        f64::from_bits(0xbbff31305b4f263e),
        f64::from_bits(0x38940fbdbf9d128a),
    ],
    [
        f64::from_bits(0x3bc3db9842b31607),
        f64::from_bits(0x385ee1ac155bab3a),
    ],
];

const LGAMMA_ACC_X0_4: [f64; 3] = [
    f64::from_bits(0x4010284e78599581),
    f64::from_bits(0xbcae78c1e9e43cfe),
    f64::from_bits(0x3932ac17bfd6be92),
];

const LGAMMA_ACC_C_4: [[f64; 2]; 24] = [
    [
        f64::from_bits(0x403aca5cf4921642),
        f64::from_bits(0x3cca46a2e0d8fdfd),
    ],
    [
        f64::from_bits(0x40044415cd813f8e),
        f64::from_bits(0x3c7afdc2672876f4),
    ],
    [
        f64::from_bits(0x3fd559b11b2a9c7c),
        f64::from_bits(0x3c617b8ada982c18),
    ],
    [
        f64::from_bits(0x3fa96d18e21aebdb),
        f64::from_bits(0xbc1c2f2da6bea629),
    ],
    [
        f64::from_bits(0x3f80261eb5732e40),
        f64::from_bits(0x3c23910ed6591782),
    ],
    [
        f64::from_bits(0x3f555e3dbf99eb3d),
        f64::from_bits(0xbbfe2d04c05513dc),
    ],
    [
        f64::from_bits(0x3f2d14fe49c4e437),
        f64::from_bits(0xbbcd3009aa4e3bb9),
    ],
    [
        f64::from_bits(0x3f0433dce282da6e),
        f64::from_bits(0xbb96a33c9bb3cb85),
    ],
    [
        f64::from_bits(0x3edc8399c7588ccf),
        f64::from_bits(0x3b66b47c3c6f3b7a),
    ],
    [
        f64::from_bits(0x3eb45fbe666d9415),
        f64::from_bits(0x3b49aaef873df6f7),
    ],
    [
        f64::from_bits(0x3e8d68d794caf0cf),
        f64::from_bits(0xbaf62f893b4b1c92),
    ],
    [
        f64::from_bits(0x3e656729dc75d1c6),
        f64::from_bits(0xbb09e28ba9f5e686),
    ],
    [
        f64::from_bits(0x3e3f5ec3352af509),
        f64::from_bits(0x3ad63d11f1c09f18),
    ],
    [
        f64::from_bits(0x3e17205756353fb8),
        f64::from_bits(0xbaae4315fc28841b),
    ],
    [
        f64::from_bits(0x3df122e7755f935b),
        f64::from_bits(0x3a5e94449e4ce532),
    ],
    [
        f64::from_bits(0x3dc982503f30184c),
        f64::from_bits(0x3a6db16eb0f7f6e4),
    ],
    [
        f64::from_bits(0x3da30f8c448c65bc),
        f64::from_bits(0xba0ba5e1bf8c5346),
    ],
    [
        f64::from_bits(0x3d7c957f8be0461c),
        f64::from_bits(0x3a078aa87d33801b),
    ],
    [
        f64::from_bits(0x3d557fe1f7e3454f),
        f64::from_bits(0xb9f61a8ccec45dca),
    ],
    [
        f64::from_bits(0x3d3035fcc6efbe53),
        f64::from_bits(0x39d454e6ac634253),
    ],
    [
        f64::from_bits(0x3d087856bbe691ed),
        f64::from_bits(0xb9692c14fb552740),
    ],
    [
        f64::from_bits(0x3ce2b877971688d2),
        f64::from_bits(0x397c72aa6bdd8e43),
    ],
    [
        f64::from_bits(0x3cbe30abef69866b),
        f64::from_bits(0x3944e95354428576),
    ],
    [
        f64::from_bits(0x3c93f5a80dfa357e),
        f64::from_bits(0x393dc7c678b87e2e),
    ],
];

const LGAMMA_ACC_X0_5: [f64; 3] = [
    f64::from_bits(0x4013f7577a6eeafd),
    f64::from_bits(0xbca5de5eab7f12cf),
    f64::from_bits(0x3914075f5e0494a2),
];

const LGAMMA_ACC_C_5: [[f64; 2]; 19] = [
    [
        f64::from_bits(0xc05d224a3ef9e41f),
        f64::from_bits(0xbcf9be272a13bcb6),
    ],
    [
        f64::from_bits(0x401b533c678a3956),
        f64::from_bits(0xbca37da6a2b62e3c),
    ],
    [
        f64::from_bits(0xbfe0d3f7fee65d34),
        f64::from_bits(0x3c8e68bf720db9be),
    ],
    [
        f64::from_bits(0x3fa752a6f5ac2726),
        f64::from_bits(0xbc116f2fc46d2a34),
    ],
    [
        f64::from_bits(0xbf713d5d163bd3f7),
        f64::from_bits(0xbc1813df737e6e21),
    ],
    [
        f64::from_bits(0x3f3a8c5c53458ca5),
        f64::from_bits(0x3bd02a507fa1273d),
    ],
    [
        f64::from_bits(0xbf05068b3ed69409),
        f64::from_bits(0xbb78fc69d833eb31),
    ],
    [
        f64::from_bits(0x3ed0ffa575ea7fe3),
        f64::from_bits(0x3b4ab6c92ac78cc4),
    ],
    [
        f64::from_bits(0xbe9bec12dd78a230),
        f64::from_bits(0x3b381392cb5ee3a0),
    ],
    [
        f64::from_bits(0x3e67382570f0b3c5),
        f64::from_bits(0x3b06c165015a71c3),
    ],
    [
        f64::from_bits(0xbe3380ebf6161cce),
        f64::from_bits(0xbacefc3933cb1d2c),
    ],
    [
        f64::from_bits(0x3e0084de43d1a947),
        f64::from_bits(0xbaa1fca46f4984d6),
    ],
    [
        f64::from_bits(0xbdcc2d90dead7a34),
        f64::from_bits(0x3a6b9c200207dd00),
    ],
    [
        f64::from_bits(0x3d982d0afe1a76a0),
        f64::from_bits(0xba2958789fbe4df9),
    ],
    [
        f64::from_bits(0xbd64d93cbb2d04c5),
        f64::from_bits(0xba0411bc9bf3a7cc),
    ],
    [
        f64::from_bits(0x3d320ea504f0c629),
        f64::from_bits(0xb9b1c1a63ab9dadd),
    ],
    [
        f64::from_bits(0xbcff6cb7ebb79bc2),
        f64::from_bits(0x395f2f3b682d06fe),
    ],
    [
        f64::from_bits(0x3ccbe565af661ea4),
        f64::from_bits(0x393b54b6ecaea8d5),
    ],
    [
        f64::from_bits(0xbc978510a36bad6c),
        f64::from_bits(0xb926935e8745afe8),
    ],
];

const LGAMMA_ACC_X0_6: [f64; 3] = [
    f64::from_bits(0x4014086a57f0b6d9),
    f64::from_bits(0x3c895262b72ca9ca),
    f64::from_bits(0x392bd98d5e0861aa),
];

const LGAMMA_ACC_C_6: [[f64; 2]; 19] = [
    [
        f64::from_bits(0x405ed72e0829ae02),
        f64::from_bits(0xbcdfdc1859ae9c53),
    ],
    [
        f64::from_bits(0x401cecc32ec22f9b),
        f64::from_bits(0x3cab6ecc779b86fc),
    ],
    [
        f64::from_bits(0x3fe253d8563f7264),
        f64::from_bits(0xbc85cd273fbd9f05),
    ],
    [
        f64::from_bits(0x3faa225df2da6e63),
        f64::from_bits(0xbc4fe9d0b1b4c660),
    ],
    [
        f64::from_bits(0x3f73e01773762671),
        f64::from_bits(0xbc1f0e1184c71c78),
    ],
    [
        f64::from_bits(0x3f3f7d8d5bdcb186),
        f64::from_bits(0xbbdd3eb6b34618fb),
    ],
    [
        f64::from_bits(0x3f09a8d00c77a92c),
        f64::from_bits(0xbb9107ec223e7938),
    ],
    [
        f64::from_bits(0x3ed557fd8c490b43),
        f64::from_bits(0x3b7d8ff4f1780ac1),
    ],
    [
        f64::from_bits(0x3ea209221a624132),
        f64::from_bits(0xbb4990b9a98394af),
    ],
    [
        f64::from_bits(0x3e6edc98d3bbeedf),
        f64::from_bits(0x3b058dea3f0bc5b9),
    ],
    [
        f64::from_bits(0x3e3aabd28e6c9a81),
        f64::from_bits(0xbadda694af0e0ff0),
    ],
    [
        f64::from_bits(0x3e073de2dd1cfd61),
        f64::from_bits(0x3aa4f223cad85a45),
    ],
    [
        f64::from_bits(0x3dd4651832f82e2e),
        f64::from_bits(0x3a398ad9bf65ef61),
    ],
    [
        f64::from_bits(0x3da200d84be7e1fc),
        f64::from_bits(0xba25a13d3eb64f54),
    ],
    [
        f64::from_bits(0x3d6ff278dbfe5630),
        f64::from_bits(0x39b7119b25262c2d),
    ],
    [
        f64::from_bits(0x3d3c77e7644c5411),
        f64::from_bits(0x39cc21cc9e02459c),
    ],
    [
        f64::from_bits(0x3d097c73d128923a),
        f64::from_bits(0xb99155bb925610d9),
    ],
    [
        f64::from_bits(0x3cd74798cd59653b),
        f64::from_bits(0xb9790790da400208),
    ],
    [
        f64::from_bits(0x3ca436d37cbd1ac2),
        f64::from_bits(0x39434f1a6c134053),
    ],
];

const LGAMMA_ACC_X0_7: [f64; 3] = [
    f64::from_bits(0x4017fe92f591f40d),
    f64::from_bits(0x3cb7dd4ed62cbd32),
    f64::from_bits(0xb932071c071a2146),
];

const LGAMMA_ACC_C_7: [[f64; 2]; 24] = [
    [
        f64::from_bits(0xc08661f6a43a5e12),
        f64::from_bits(0xbd20c437b83bc0dd),
    ],
    [
        f64::from_bits(0x404f79dcb794f26f),
        f64::from_bits(0xbcbada8018ade6ef),
    ],
    [
        f64::from_bits(0xc01d6e8088a19ffe),
        f64::from_bits(0xbca2c08708631269),
    ],
    [
        f64::from_bits(0x3feef5d308dbfc97),
        f64::from_bits(0x3c587cd5c0d349fa),
    ],
    [
        f64::from_bits(0xbfc15ea6b0ab529e),
        f64::from_bits(0x3bd4ac7004538519),
    ],
    [
        f64::from_bits(0x3f944d54e9fe2397),
        f64::from_bits(0x3c1f1ca84e1a96c0),
    ],
    [
        f64::from_bits(0xbf68684e40cebb3d),
        f64::from_bits(0xbbf19803482caec7),
    ],
    [
        f64::from_bits(0x3f3df44c1d81c723),
        f64::from_bits(0xbbc0ed2e1ad825e8),
    ],
    [
        f64::from_bits(0xbf12ac3053f4ee18),
        f64::from_bits(0xbbbfe50dd710dd6b),
    ],
    [
        f64::from_bits(0x3ee79226ae04a847),
        f64::from_bits(0xbb4379839e0000d6),
    ],
    [
        f64::from_bits(0xbebe0dffb5f77ccb),
        f64::from_bits(0x3b5a4015df3482b4),
    ],
    [
        f64::from_bits(0x3e935217890708e7),
        f64::from_bits(0x3b12875bd051d55e),
    ],
    [
        f64::from_bits(0xbe6903aa9af5cf03),
        f64::from_bits(0xbaca619f5e70dc34),
    ],
    [
        f64::from_bits(0x3e404a103323cd6a),
        f64::from_bits(0x3aca260aa452be23),
    ],
    [
        f64::from_bits(0xbe1552efc1ecd295),
        f64::from_bits(0xbab88a83f178c250),
    ],
    [
        f64::from_bits(0x3dec0a12b600b141),
        f64::from_bits(0xba6c3732154b52ac),
    ],
    [
        f64::from_bits(0xbdc281ce75ae916f),
        f64::from_bits(0xba6a96d8e0798081),
    ],
    [
        f64::from_bits(0x3d988411a65d3af4),
        f64::from_bits(0xba30e8c83eb281a9),
    ],
    [
        f64::from_bits(0xbd7049e79a4b5baa),
        f64::from_bits(0xba164caa9cfd44a2),
    ],
    [
        f64::from_bits(0x3d45b0f16f419f7b),
        f64::from_bits(0x39974e160750d876),
    ],
    [
        f64::from_bits(0xbd1ce70ccad93bc1),
        f64::from_bits(0x39b1b340ff57e754),
    ],
    [
        f64::from_bits(0x3cf3a7c53f75547f),
        f64::from_bits(0xb9988a8d5ad25f3f),
    ],
    [
        f64::from_bits(0xbccc41c8d3e4c2c9),
        f64::from_bits(0x396e4dd9e959c88b),
    ],
    [
        f64::from_bits(0x3c9f82da393ec32a),
        f64::from_bits(0xb8c65c64b7f7e325),
    ],
];

const LGAMMA_ACC_X0_8: [f64; 3] = [
    f64::from_bits(0x4018016b25897c8d),
    f64::from_bits(0xbc927e0f49a4ba72),
    f64::from_bits(0x39172e1ab15a4d03),
];

const LGAMMA_ACC_C_8: [[f64; 2]; 24] = [
    [
        f64::from_bits(0x40869de49e3af2aa),
        f64::from_bits(0x3d0954b690943afc),
    ],
    [
        f64::from_bits(0x404fce23484cfd10),
        f64::from_bits(0x3ce8266e7580f08c),
    ],
    [
        f64::from_bits(0x401de503a3c37c40),
        f64::from_bits(0x3ca9fa7459d62427),
    ],
    [
        f64::from_bits(0x3fef9c7b52558abb),
        f64::from_bits(0x3c8b6896749474c1),
    ],
    [
        f64::from_bits(0x3fc1d3d50714416a),
        f64::from_bits(0x3c556020c7693337),
    ],
    [
        f64::from_bits(0x3f94f21e2fb9e060),
        f64::from_bits(0x3c39e1cfc6a7c130),
    ],
    [
        f64::from_bits(0x3f69500994cd8a9e),
        f64::from_bits(0xbc0ce1777b097bd6),
    ],
    [
        f64::from_bits(0x3f3f3a2c23c19d79),
        f64::from_bits(0xbb9972494ba4cd4d),
    ],
    [
        f64::from_bits(0x3f139152652eb3aa),
        f64::from_bits(0x3ba8f85c4250ca80),
    ],
    [
        f64::from_bits(0x3ee8d45f8be891de),
        f64::from_bits(0x3b841e12761b265a),
    ],
    [
        f64::from_bits(0x3ebfd3214a702b29),
        f64::from_bits(0xbb5140076f8887f9),
    ],
    [
        f64::from_bits(0x3e9490b476815d46),
        f64::from_bits(0xbae3c92eb9376d39),
    ],
    [
        f64::from_bits(0x3e6ac3b965280088),
        f64::from_bits(0xbad468539ffe88b5),
    ],
    [
        f64::from_bits(0x3e41851c43c40d43),
        f64::from_bits(0xbad6918e4fb64554),
    ],
    [
        f64::from_bits(0x3e170dfb529f2c5a),
        f64::from_bits(0xbabb5dfdc2a8890d),
    ],
    [
        f64::from_bits(0x3dee791ef4ef8602),
        f64::from_bits(0x3a5b20188eabd9e6),
    ],
    [
        f64::from_bits(0x3dc437e60f850abe),
        f64::from_bits(0x3a6b5786abafa12e),
    ],
    [
        f64::from_bits(0x3d9aec29bcf6c8c1),
        f64::from_bits(0x3a3774aed905f5cf),
    ],
    [
        f64::from_bits(0x3d71fb22834000f1),
        f64::from_bits(0x3a0705349812bcba),
    ],
    [
        f64::from_bits(0x3d4811ce4ec7b869),
        f64::from_bits(0x39c8be5e35bce704),
    ],
    [
        f64::from_bits(0x3d201e72b00fd763),
        f64::from_bits(0x39cc8c867f1090e3),
    ],
    [
        f64::from_bits(0x3cf60a0a9eee14d8),
        f64::from_bits(0x395ed307bf965e72),
    ],
    [
        f64::from_bits(0x3ccfdce60428ba59),
        f64::from_bits(0x395268eaf9e15016),
    ],
    [
        f64::from_bits(0x3ca1dcf4e3ac0033),
        f64::from_bits(0x393bdbdd62565840),
    ],
];

const LGAMMA_ACC_X0_9: [f64; 3] = [
    f64::from_bits(0x401bffcbf76b86f0),
    f64::from_bits(0xbc6853b29347b806),
    f64::from_bits(0x3900fa018051dd41),
];

const LGAMMA_ACC_C_9: [[f64; 2]; 25] = [
    [
        f64::from_bits(0xc0b3abf7a5cea91b),
        f64::from_bits(0xbd58257b8abd02ce),
    ],
    [
        f64::from_bits(0x4088349a2550422d),
        f64::from_bits(0xbd2c6f2ef4105b99),
    ],
    [
        f64::from_bits(0xc063d91dadc98428),
        f64::from_bits(0x3cf46605ff1fc41c),
    ],
    [
        f64::from_bits(0x40424f3d636f3339),
        f64::from_bits(0x3ce5966a8447ea56),
    ],
    [
        f64::from_bits(0xc0220427df1b3492),
        f64::from_bits(0xbcbe9801a53bb06e),
    ],
    [
        f64::from_bits(0x4002775e857fb69c),
        f64::from_bits(0x3ca88bd494f9ea51),
    ],
    [
        f64::from_bits(0xbfe377e70b463c13),
        f64::from_bits(0xbc82524c9cdd93db),
    ],
    [
        f64::from_bits(0x3fc4f3d28edba5cd),
        f64::from_bits(0x3c392d2ac40dac77),
    ],
    [
        f64::from_bits(0xbfa6e8557168cf8a),
        f64::from_bits(0x3c40b5cf2a57fc71),
    ],
    [
        f64::from_bits(0x3f895bb17ce427a4),
        f64::from_bits(0xbbd0f3a87df6aaa1),
    ],
    [
        f64::from_bits(0xbf6c5ac12d48f7e4),
        f64::from_bits(0x3c01d042aabfd8fb),
    ],
    [
        f64::from_bits(0x3f4ff816dad74e75),
        f64::from_bits(0x3bc24594bf64f76d),
    ],
    [
        f64::from_bits(0xbf3225f4a6a0c5b4),
        f64::from_bits(0xbbde1dd93bdafdda),
    ],
    [
        f64::from_bits(0x3f14ba3e5c007209),
        f64::from_bits(0x3b9a47689fcd5098),
    ],
    [
        f64::from_bits(0xbef7cb738040762c),
        f64::from_bits(0xbb9f100cc43c4b52),
    ],
    [
        f64::from_bits(0x3edb701200fca8b8),
        f64::from_bits(0x3b7d825b96ec4ceb),
    ],
    [
        f64::from_bits(0xbebfc33bee21552c),
        f64::from_bits(0x3b4096b388c1821d),
    ],
    [
        f64::from_bits(0x3ea272cc14fa7a53),
        f64::from_bits(0x3b3d22fcb28979f7),
    ],
    [
        f64::from_bits(0xbe857f61caeb3e48),
        f64::from_bits(0xbb003971054f5156),
    ],
    [
        f64::from_bits(0x3e691f102889ab68),
        f64::from_bits(0xbafb9949bdf154e2),
    ],
    [
        f64::from_bits(0xbe4d641c300c655e),
        f64::from_bits(0xbad47bf835c77fc5),
    ],
    [
        f64::from_bits(0x3e31319b4831f0e0),
        f64::from_bits(0xbab2688b8ef3d8f7),
    ],
    [
        f64::from_bits(0xbe14c01654897d29),
        f64::from_bits(0xba95c8d35e08cb51),
    ],
    [
        f64::from_bits(0x3dfa7ea9eb248f1e),
        f64::from_bits(0xba9346b798e3ef55),
    ],
    [
        f64::from_bits(0xbdd8afe498da415b),
        f64::from_bits(0x3a73f02ddea04afb),
    ],
];

const LGAMMA_ACC_X0_10: [f64; 3] = [
    f64::from_bits(0x401c0033fdedfe1f),
    f64::from_bits(0xbcb20bb7d2324678),
    f64::from_bits(0xb95f5536678d69d3),
];

const LGAMMA_ACC_C_10: [[f64; 2]; 24] = [
    [
        f64::from_bits(0x40b3b407aa387bd1),
        f64::from_bits(0x3d4da1e57343b1db),
    ],
    [
        f64::from_bits(0x40783e85daafbad6),
        f64::from_bits(0xbd1f37538d99bd29),
    ],
    [
        f64::from_bits(0x4043e552b5e3c226),
        f64::from_bits(0xbce07b1550d11cf4),
    ],
    [
        f64::from_bits(0x40125e42a45e905b),
        f64::from_bits(0x3c961a63054c47c5),
    ],
    [
        f64::from_bits(0x3fe216a3560743ee),
        f64::from_bits(0x3c6f5024887fdb9d),
    ],
    [
        f64::from_bits(0x3fb28e1c70ef5313),
        f64::from_bits(0x3c43234af5080105),
    ],
    [
        f64::from_bits(0x3f8393e2bc330081),
        f64::from_bits(0xbc145ef9ada326ba),
    ],
    [
        f64::from_bits(0x3f55164141f5ae6a),
        f64::from_bits(0xbbe802404a8fab89),
    ],
    [
        f64::from_bits(0x3f2712b3a86e1bdf),
        f64::from_bits(0xbbab2f4283e8796e),
    ],
    [
        f64::from_bits(0x3ef98fd36b906e08),
        f64::from_bits(0xbb9c84ad8d8481b7),
    ],
    [
        f64::from_bits(0x3ecc9ae6ef6262f6),
        f64::from_bits(0x3b6ff0ee333033bf),
    ],
    [
        f64::from_bits(0x3ea02382a958b63d),
        f64::from_bits(0xbb499f49c2b45652),
    ],
    [
        f64::from_bits(0x3e7256845ec9b6b1),
        f64::from_bits(0x3af822bd9a29ab32),
    ],
    [
        f64::from_bits(0x3e44f5ff36087486),
        f64::from_bits(0xbaddc262e92f86f3),
    ],
    [
        f64::from_bits(0x3e1814f9ce48442f),
        f64::from_bits(0x3ab544e6fd9ac4be),
    ],
    [
        f64::from_bits(0x3debca89f844ac84),
        f64::from_bits(0xba863ceb4f54c669),
    ],
    [
        f64::from_bits(0x3dc01946c9b7b26a),
        f64::from_bits(0xba59b4cb6a88b6fc),
    ],
    [
        f64::from_bits(0x3d92b75936f37bc3),
        f64::from_bits(0xba1d8809500cd622),
    ],
    [
        f64::from_bits(0x3d65d3d0a2486470),
        f64::from_bits(0xba0061c56c047c1d),
    ],
    [
        f64::from_bits(0x3d398297f709854f),
        f64::from_bits(0x3984cc4272ef2e04),
    ],
    [
        f64::from_bits(0x3d0dd4cc8e9d77b5),
        f64::from_bits(0x3990103a0db6c049),
    ],
    [
        f64::from_bits(0x3ce1ce46c495fa62),
        f64::from_bits(0x398622087cd95f2a),
    ],
    [
        f64::from_bits(0x3cb678ec65743609),
        f64::from_bits(0xb94850d52cb53a9f),
    ],
    [
        f64::from_bits(0x3c85fe290153bcaa),
        f64::from_bits(0x392a8b6414621298),
    ],
];

const LGAMMA_ACC_X0_11: [f64; 3] = [
    f64::from_bits(0x401ffff97f8159cf),
    f64::from_bits(0x3c8e54f415a91586),
    f64::from_bits(0x39253a5d106f9a3e),
];

const LGAMMA_ACC_C_11: [[f64; 2]; 27] = [
    [
        f64::from_bits(0xc0e3af76fe4c2fab),
        f64::from_bits(0xbd77cc92f0b999e0),
    ],
    [
        f64::from_bits(0x40b838e76caaf123),
        f64::from_bits(0x3d5292e15f502d1f),
    ],
    [
        f64::from_bits(0xc093de68b3256526),
        f64::from_bits(0x3d35456a490ea6e3),
    ],
    [
        f64::from_bits(0x407255c052530c71),
        f64::from_bits(0xbcf6700425996a73),
    ],
    [
        f64::from_bits(0xc0520c2a8418126a),
        f64::from_bits(0x3cd1d5e3d82023c1),
    ],
    [
        f64::from_bits(0x40328139342cef00),
        f64::from_bits(0x3cc0238cfdd7475d),
    ],
    [
        f64::from_bits(0xc01384066c322246),
        f64::from_bits(0x3cace82f225b48d1),
    ],
    [
        f64::from_bits(0x3ff502bc4dad47d3),
        f64::from_bits(0x3c5f45416fab08da),
    ],
    [
        f64::from_bits(0xbfd6faadfece0e30),
        f64::from_bits(0xbc76c5874c40cd48),
    ],
    [
        f64::from_bits(0x3fb9724323c89909),
        f64::from_bits(0xbc527460acd35293),
    ],
    [
        f64::from_bits(0xbf9c7684c96f2430),
        f64::from_bits(0x3c15d749ebbd28fe),
    ],
    [
        f64::from_bits(0x3f800d1f48743fb8),
        f64::from_bits(0x3c2ad60125333897),
    ],
    [
        f64::from_bits(0xbf623af6dd470e5c),
        f64::from_bits(0xbbf1d084246f0eda),
    ],
    [
        f64::from_bits(0x3f44d416958eec41),
        f64::from_bits(0xbbcdfeda2a7c037e),
    ],
    [
        f64::from_bits(0xbf27eb3eb405269f),
        f64::from_bits(0xbba95169511f695d),
    ],
    [
        f64::from_bits(0x3f0b972ec5a016f0),
        f64::from_bits(0x3ba0a9a114630056),
    ],
    [
        f64::from_bits(0xbeeff35af33eaac0),
        f64::from_bits(0xbb873922358a5b41),
    ],
    [
        f64::from_bits(0x3ed290662f9abea1),
        f64::from_bits(0xbb78a5267d544bb7),
    ],
    [
        f64::from_bits(0xbeb5a395633f9d69),
        f64::from_bits(0xbb3377120c488d2e),
    ],
    [
        f64::from_bits(0x3e994b2a8f6b566e),
        f64::from_bits(0x3b0685bd8f7318c7),
    ],
    [
        f64::from_bits(0xbe7da517a6c7ab89),
        f64::from_bits(0xbb1940f32a9a6413),
    ],
    [
        f64::from_bits(0x3e6168992efa74a9),
        f64::from_bits(0xbaebb9b7bf492795),
    ],
    [
        f64::from_bits(0xbe44665c665b014e),
        f64::from_bits(0x3ae91ef4857aca34),
    ],
    [
        f64::from_bits(0x3e27f49e136d0e91),
        f64::from_bits(0x3aa9d7da1584fe90),
    ],
    [
        f64::from_bits(0xbe0df48e4db635ef),
        f64::from_bits(0xbaaac6dd06d9a112),
    ],
    [
        f64::from_bits(0x3df3665697708dd0),
        f64::from_bits(0xba5e15f660c03e5b),
    ],
    [
        f64::from_bits(0xbdd015f65348c21a),
        f64::from_bits(0xba5e45b91e0267c0),
    ],
];

const LGAMMA_ACC_X0_12: [f64; 3] = [
    f64::from_bits(0x402000034028b3f9),
    f64::from_bits(0x3cbf60cb3cec1ced),
    f64::from_bits(0xb95ea26620d6b1ca),
];

const LGAMMA_ACC_C_12: [[f64; 2]; 27] = [
    [
        f64::from_bits(0x40e3b088fed67718),
        f64::from_bits(0xbd8505613ba29893),
    ],
    [
        f64::from_bits(0x40b83a3893550edc),
        f64::from_bits(0x3d5f52e3b240db30),
    ],
    [
        f64::from_bits(0x4093e0078db8ada4),
        f64::from_bits(0x3d3506573eecc6e3),
    ],
    [
        f64::from_bits(0x407257bec9464251),
        f64::from_bits(0x3ce8c4ea8394aa1e),
    ],
    [
        f64::from_bits(0x40520e9ea0755a47),
        f64::from_bits(0xbcf978e59e1a3d9f),
    ],
    [
        f64::from_bits(0x4032843e1313c83b),
        f64::from_bits(0xbc850610d6737717),
    ],
    [
        f64::from_bits(0x401387bd6a785478),
        f64::from_bits(0xbcb1f1fff3fa1b4b),
    ],
    [
        f64::from_bits(0x3ff5074e788de770),
        f64::from_bits(0x3c8aeeb83580ea25),
    ],
    [
        f64::from_bits(0x3fd7004dd990d7da),
        f64::from_bits(0x3c43f7f554ffc670),
    ],
    [
        f64::from_bits(0x3fb9792ed5f6dfb4),
        f64::from_bits(0xbc51f23ec3a02887),
    ],
    [
        f64::from_bits(0x3f9c7f08cdaef32f),
        f64::from_bits(0xbc11e18e4c767de0),
    ],
    [
        f64::from_bits(0x3f80125c81121f59),
        f64::from_bits(0xbc2a6731270f4d3f),
    ],
    [
        f64::from_bits(0x3f62416931f22d22),
        f64::from_bits(0xbbf2bc2867243e45),
    ],
    [
        f64::from_bits(0x3f44dc0543bea659),
        f64::from_bits(0x3bd635ca3a9b1a8d),
    ],
    [
        f64::from_bits(0x3f27f501645b2fde),
        f64::from_bits(0xbbce45b9db4c7680),
    ],
    [
        f64::from_bits(0x3f0ba331549d971d),
        f64::from_bits(0x3ba97379411be4a6),
    ],
    [
        f64::from_bits(0x3ef001110caa210f),
        f64::from_bits(0x3b968ce3f9f57bc4),
    ],
    [
        f64::from_bits(0x3ed2997db55cba74),
        f64::from_bits(0xbb50f0fa5ac41b0a),
    ],
    [
        f64::from_bits(0x3eb5aec550050044),
        f64::from_bits(0x3b57fdca933bee56),
    ],
    [
        f64::from_bits(0x3e9958ee8d695040),
        f64::from_bits(0x3b3c28e192bc3ba4),
    ],
    [
        f64::from_bits(0x3e7db608c6ccccc8),
        f64::from_bits(0xbaf71a440ef759be),
    ],
    [
        f64::from_bits(0x3e6173059c7e38cf),
        f64::from_bits(0x3b040d18009e343f),
    ],
    [
        f64::from_bits(0x3e44731fe682342c),
        f64::from_bits(0xbac535cdcf9e608b),
    ],
    [
        f64::from_bits(0x3e28043ec54aba15),
        f64::from_bits(0xbab62d6d36f1d19e),
    ],
    [
        f64::from_bits(0x3e0e08fc2d5d4a00),
        f64::from_bits(0xbaa6e4b1587c874e),
    ],
    [
        f64::from_bits(0x3df3742f9e468b63),
        f64::from_bits(0x3a74c5c5de4449bf),
    ],
    [
        f64::from_bits(0x3dd021df5d9d312d),
        f64::from_bits(0x3a5b2da076b6955f),
    ],
];

const LGAMMA_ACC_X0_13: [f64; 3] = [
    f64::from_bits(0x4021ffffa3884bd0),
    f64::from_bits(0x3caff90c9d2ae925),
    f64::from_bits(0xb9430c0efef78c04),
];

const LGAMMA_ACC_C_13: [[f64; 2]; 28] = [
    [
        f64::from_bits(0xc11625edfc63db2f),
        f64::from_bits(0x3dada7fc3ed68f99),
    ],
    [
        f64::from_bits(0x40eea8c150480a7a),
        f64::from_bits(0x3d8344e4cbf6d2a1),
    ],
    [
        f64::from_bits(0xc0cc4b30e4bc55c1),
        f64::from_bits(0xbd69ec40fdd810d7),
    ],
    [
        f64::from_bits(0x40ad5fe468dbbf03),
        f64::from_bits(0xbd480705f1856688),
    ],
    [
        f64::from_bits(0xc09043d21bc24dec),
        f64::from_bits(0xbd2b0db95c85f28b),
    ],
    [
        f64::from_bits(0x4072c334ae535e1d),
        f64::from_bits(0x3d1530b0b83f8e20),
    ],
    [
        f64::from_bits(0xc0564314b431cd64),
        f64::from_bits(0xbce5d7e4374b5c05),
    ],
    [
        f64::from_bits(0x403af6ed589b3a86),
        f64::from_bits(0xbcc622e836ea2287),
    ],
    [
        f64::from_bits(0xc02096e446edcfb4),
        f64::from_bits(0x3ccf849c504338f9),
    ],
    [
        f64::from_bits(0x4004aaf49e713c0c),
        f64::from_bits(0x3cabbab832e792a1),
    ],
    [
        f64::from_bits(0xbfea0246d9c1b57f),
        f64::from_bits(0xbc718f4c5d12b8b5),
    ],
    [
        f64::from_bits(0x3fd0806315c1a365),
        f64::from_bits(0x3c5074a0bafc599f),
    ],
    [
        f64::from_bits(0xbfb515dd6b891094),
        f64::from_bits(0xbc4bcca3edb6dc58),
    ],
    [
        f64::from_bits(0x3f9b1a5fe76de3d4),
        f64::from_bits(0xbc0f8d7b8e44c315),
    ],
    [
        f64::from_bits(0xbf8182229539b98c),
        f64::from_bits(0x3c15c2ae71e11bfe),
    ],
    [
        f64::from_bits(0x3f66b8b31630a98c),
        f64::from_bits(0xbbec5c500af742da),
    ],
    [
        f64::from_bits(0xbf4d9a3c6789b8e8),
        f64::from_bits(0x3be7160288934ba5),
    ],
    [
        f64::from_bits(0x3f3359c2ea889b3e),
        f64::from_bits(0x3bae44b62289b60b),
    ],
    [
        f64::from_bits(0xbf19606a9f683511),
        f64::from_bits(0xbbb3a3c7889c7af5),
    ],
    [
        f64::from_bits(0x3f00af84fb66238e),
        f64::from_bits(0xbba3d487a70215da),
    ],
    [
        f64::from_bits(0xbee5ff88107d908a),
        f64::from_bits(0x3b60c4172ceada86),
    ],
    [
        f64::from_bits(0x3ecd14070a537453),
        f64::from_bits(0x3b6964ddb3dba2ba),
    ],
    [
        f64::from_bits(0xbeb33eaa64e4654e),
        f64::from_bits(0x3b5e6deb19e51b4e),
    ],
    [
        f64::from_bits(0x3e9957ef375ea2ab),
        f64::from_bits(0x3b2065699585b36c),
    ],
    [
        f64::from_bits(0xbe80cd5c9d34a242),
        f64::from_bits(0xbb26f2ed041b79a1),
    ],
    [
        f64::from_bits(0x3e6836d31552e5d6),
        f64::from_bits(0xbb0789b28b32b093),
    ],
    [
        f64::from_bits(0xbe519aeee4b34ef6),
        f64::from_bits(0xbaec45d882aaae00),
    ],
    [
        f64::from_bits(0x3e2eb1912eaced35),
        f64::from_bits(0x3abe9842f288e554),
    ],
];

const LGAMMA_ACC_X0_14: [f64; 3] = [
    f64::from_bits(0x402200005c7768fb),
    f64::from_bits(0x3c9b5b610ffb70d4),
    f64::from_bits(0x393deb7ad09ec5ea),
];

const LGAMMA_ACC_C_14: [[f64; 2]; 28] = [
    [
        f64::from_bits(0x4116261203919440),
        f64::from_bits(0x3d97d5e8272ce41d),
    ],
    [
        f64::from_bits(0x40eea8f32fb7f586),
        f64::from_bits(0xbd8345b1cc20d4a5),
    ],
    [
        f64::from_bits(0x40cc4b75ee68e2ba),
        f64::from_bits(0xbd5812d7bcedccb5),
    ],
    [
        f64::from_bits(0x40ad6043fa1ffaa5),
        f64::from_bits(0xbd45a4eaff5f008a),
    ],
    [
        f64::from_bits(0x40904414411db7f4),
        f64::from_bits(0x3d3d742c54011402),
    ],
    [
        f64::from_bits(0x4072c3903ec9c90c),
        f64::from_bits(0x3d173da27a4b1b20),
    ],
    [
        f64::from_bits(0x40564393744bb9bd),
        f64::from_bits(0xbcf482bcb016f26f),
    ],
    [
        f64::from_bits(0x403af79ccdc71d33),
        f64::from_bits(0x3cd0a7a439456745),
    ],
    [
        f64::from_bits(0x4020975db7d71fc6),
        f64::from_bits(0x3ccbd477ffdd39f6),
    ],
    [
        f64::from_bits(0x4004ab9cba1e3478),
        f64::from_bits(0x3c91c929f5a50618),
    ],
    [
        f64::from_bits(0x3fea032f8f114635),
        f64::from_bits(0x3c892279b6d57eb2),
    ],
    [
        f64::from_bits(0x3fd0810426bfa5a1),
        f64::from_bits(0x3c70cd3cd3381e06),
    ],
    [
        f64::from_bits(0x3fb516bc616eaf53),
        f64::from_bits(0xbc50ebc405b4871a),
    ],
    [
        f64::from_bits(0x3f9b1b948b11a02d),
        f64::from_bits(0x3c38b0ffc5af8ce1),
    ],
    [
        f64::from_bits(0x3f8182f8343c9e61),
        f64::from_bits(0xbc2de95c6fadf8bb),
    ],
    [
        f64::from_bits(0x3f66b9dacc2e40d0),
        f64::from_bits(0x3c07d2afb82825d1),
    ],
    [
        f64::from_bits(0x3f4d9bd5c016a05a),
        f64::from_bits(0x3beafdd41ec7d841),
    ],
    [
        f64::from_bits(0x3f335ade3d85cd0e),
        f64::from_bits(0xbbd5080702a55ab7),
    ],
    [
        f64::from_bits(0x3f1961f2d2659b02),
        f64::from_bits(0x3bbf377d9fd7df3a),
    ],
    [
        f64::from_bits(0x3f00b0946f97b203),
        f64::from_bits(0x3ba9ecb24ffa7674),
    ],
    [
        f64::from_bits(0x3ee600ffd56f5733),
        f64::from_bits(0x3b82ffa0ff103bfd),
    ],
    [
        f64::from_bits(0x3ecd160f5b944255),
        f64::from_bits(0xbb46f04901227877),
    ],
    [
        f64::from_bits(0x3eb3401289d10771),
        f64::from_bits(0xbb3b2fca86443a23),
    ],
    [
        f64::from_bits(0x3e9959dec69f88ed),
        f64::from_bits(0x3b31945258ad2d80),
    ],
    [
        f64::from_bits(0x3e80ceb1d5fdcfdb),
        f64::from_bits(0xbb2e1bf43eb8d61f),
    ],
    [
        f64::from_bits(0x3e6838cdbbcfd467),
        f64::from_bits(0xbaf948ada050fa52),
    ],
    [
        f64::from_bits(0x3e519c71e9accac3),
        f64::from_bits(0x3ae74cb4459f6b90),
    ],
    [
        f64::from_bits(0x3e2eb4678e3fcf16),
        f64::from_bits(0xbac66454030f804f),
    ],
];

const LGAMMA_ACC_X0_15: [f64; 3] = [
    f64::from_bits(0x4023fffff6c0d7c0),
    f64::from_bits(0xbcc197cea8c42d7d),
    f64::from_bits(0xb967072c5a292198),
];

const LGAMMA_ACC_C_15: [[f64; 2]; 27] = [
    [
        f64::from_bits(0xc14baf7da5f3795d),
        f64::from_bits(0xbde16a79518c8367),
    ],
    [
        f64::from_bits(0x4117f3e8791fa0d2),
        f64::from_bits(0xbdb2aec811c9609e),
    ],
    [
        f64::from_bits(0xc0eba18befcaaa63),
        f64::from_bits(0xbd8d18c4e2604650),
    ],
    [
        f64::from_bits(0x40c1ede14765dc0c),
        f64::from_bits(0x3d613bc952384f0f),
    ],
    [
        f64::from_bits(0xc098d1a9ab5a5050),
        f64::from_bits(0xbd1904ee4f738385),
    ],
    [
        f64::from_bits(0x4071e4d8c35d22cc),
        f64::from_bits(0xbcf3b56f0fecec54),
    ],
    [
        f64::from_bits(0xc04a8a191db10900),
        f64::from_bits(0xbcc11eec467fa814),
    ],
    [
        f64::from_bits(0x4024174f65ff8680),
        f64::from_bits(0x3ccc939d4bfd1404),
    ],
    [
        f64::from_bits(0xbffee6d90f2332c7),
        f64::from_bits(0x3c83b18f44783de9),
    ],
    [
        f64::from_bits(0x3fd80fd3420fba09),
        f64::from_bits(0x3c655812b31052fb),
    ],
    [
        f64::from_bits(0xbfb2ecd481762eae),
        f64::from_bits(0xbc1d974aa28e0c70),
    ],
    [
        f64::from_bits(0x3f8e04a0b28db4c5),
        f64::from_bits(0xbc2f0d1f410ed5c5),
    ],
    [
        f64::from_bits(0xbf67f91af3f1200f),
        f64::from_bits(0x3c0da2dc770582dd),
    ],
    [
        f64::from_bits(0x3f4342652fcfe1b6),
        f64::from_bits(0xbbd7744fdcaebcfa),
    ],
    [
        f64::from_bits(0xbf1f1a88f796a348),
        f64::from_bits(0x3b5a8bb92c660ef4),
    ],
    [
        f64::from_bits(0x3ef93a68769f8a0f),
        f64::from_bits(0xbb9f24bb59f0e5da),
    ],
    [
        f64::from_bits(0xbed48af467f6d621),
        f64::from_bits(0x3b57948f005903ca),
    ],
    [
        f64::from_bits(0x3eb0c92181e40e32),
        f64::from_bits(0x3b565951755e7adf),
    ],
    [
        f64::from_bits(0xbe8b842459ce40bf),
        f64::from_bits(0x3b24503ca5bb9422),
    ],
    [
        f64::from_bits(0x3e669db729e3b047),
        f64::from_bits(0xbb0c0423483dea32),
    ],
    [
        f64::from_bits(0xbe42a3762889a1eb),
        f64::from_bits(0xbaec935e18265a70),
    ],
    [
        f64::from_bits(0x3e1ec8fbfde50a33),
        f64::from_bits(0x3a698aa5e946ac0f),
    ],
    [
        f64::from_bits(0xbdf95dcfae4547bc),
        f64::from_bits(0xba8c04aac4c4ba01),
    ],
    [
        f64::from_bits(0x3dd4f21457a90c36),
        f64::from_bits(0xba7e01d7e3277c24),
    ],
    [
        f64::from_bits(0xbdb26aacbe779418),
        f64::from_bits(0x3a3e6b067bac93ac),
    ],
    [
        f64::from_bits(0x3d90c603e5994262),
        f64::from_bits(0x3a2d6144d3f2ef17),
    ],
    [
        f64::from_bits(0xbd638f713b1343d3),
        f64::from_bits(0xba00c3eef336946b),
    ],
];

const LGAMMA_ACC_X0_16: [f64; 3] = [
    f64::from_bits(0x40240000093f2777),
    f64::from_bits(0x3cb927b45d95e154),
    f64::from_bits(0x3950780c21b6e452),
];

const LGAMMA_ACC_C_16: [[f64; 2]; 27] = [
    [
        f64::from_bits(0x414baf825a0c63b2),
        f64::from_bits(0xbdc20323f1015cdf),
    ],
    [
        f64::from_bits(0x4117f3ec8ae05f2e),
        f64::from_bits(0x3db2aec80d23ccb8),
    ],
    [
        f64::from_bits(0x40eba192fa62a5c8),
        f64::from_bits(0xbd825660ae99a81c),
    ],
    [
        f64::from_bits(0x40c1ede75ef431b0),
        f64::from_bits(0xbd6a691cbd9c9beb),
    ],
    [
        f64::from_bits(0x4098d1b435ece20f),
        f64::from_bits(0x3d05b052a017c6d3),
    ],
    [
        f64::from_bits(0x4071e4e1e218c99c),
        f64::from_bits(0x3d1775895e089534),
    ],
    [
        f64::from_bits(0x404a8a28e596ccce),
        f64::from_bits(0xbceaa094ee6e63a1),
    ],
    [
        f64::from_bits(0x4024175d0d35b3d4),
        f64::from_bits(0xbcc0cf2eb638feb3),
    ],
    [
        f64::from_bits(0x3ffee6f0af10b985),
        f64::from_bits(0xbc97d3173e84c277),
    ],
    [
        f64::from_bits(0x3fd80fe7b2913e68),
        f64::from_bits(0x3c758aa6aa5a2a1d),
    ],
    [
        f64::from_bits(0x3fb2ece6307c7cb0),
        f64::from_bits(0x3c5812bbc72b9bf4),
    ],
    [
        f64::from_bits(0x3f8e04bf4be02590),
        f64::from_bits(0x3c1b2fdc34bbffbe),
    ],
    [
        f64::from_bits(0x3f67f9356d1f8f5d),
        f64::from_bits(0x3c0200577f2d649b),
    ],
    [
        f64::from_bits(0x3f43427c173faa4e),
        f64::from_bits(0x3be963e715185901),
    ],
    [
        f64::from_bits(0x3f1f1ab0995dd7aa),
        f64::from_bits(0xbb9982b6a256e922),
    ],
    [
        f64::from_bits(0x3ef93a8ac07adde7),
        f64::from_bits(0x3b91938f0da933ad),
    ],
    [
        f64::from_bits(0x3ed48b1212550a17),
        f64::from_bits(0xbb773810e7394dec),
    ],
    [
        f64::from_bits(0x3eb0c93b2c55f9b3),
        f64::from_bits(0x3b4b11c8bd5d36ca),
    ],
    [
        f64::from_bits(0x3e8b8450c2eb24f2),
        f64::from_bits(0x3acb77d60943d51a),
    ],
    [
        f64::from_bits(0x3e669ddd961f0eb9),
        f64::from_bits(0x3af26555b649e2b9),
    ],
    [
        f64::from_bits(0x3e42a39767b036bd),
        f64::from_bits(0xbadff547cf250c2c),
    ],
    [
        f64::from_bits(0x3e1ec935873525f3),
        f64::from_bits(0xbab56d39c71ef428),
    ],
    [
        f64::from_bits(0x3df95e014a90402e),
        f64::from_bits(0x3a7ff041edd80219),
    ],
    [
        f64::from_bits(0x3dd4f23f078ebd25),
        f64::from_bits(0x3a622e3191f74383),
    ],
    [
        f64::from_bits(0x3db26ad387c7d522),
        f64::from_bits(0xba40581274dbf440),
    ],
    [
        f64::from_bits(0x3d90c628e2dc8b58),
        f64::from_bits(0x39e4473a911adb21),
    ],
    [
        f64::from_bits(0x3d638f9f8f1d99a8),
        f64::from_bits(0xba06893c733787d0),
    ],
];
#[inline(never)]
fn as_lgamma_database(x: f64, mut f: f64) -> f64 {
    let mut a: i32 = 0;
    let mut b: i32 = (LGAMMA_DB.len() as i32) - 1;
    let mut m = (a + b) / 2;
    while a <= b {
        let v = LGAMMA_DB[m as usize][0];
        if v < x {
            a = m + 1;
        } else if v == x {
            let row = LGAMMA_DB[m as usize];
            f = row[1] + row[2];
            break;
        } else {
            b = m - 1;
        }
        m = (a + b) / 2;
    }
    f
}

fn as_logd(x: f64, l: &mut f64) -> f64 {
    let mut t = (x).to_bits();
    let mut ex = (t >> 52) as i32;
    if ex == 0 {
        let k = t.leading_zeros() as i32;
        t <<= (k - 11) as u32;
        ex -= k - 12;
    }
    let e = ex - 0x3ff;
    t &= (!0u64) >> 12;
    let ed = e as f64;
    let i = (t >> (52 - 5)) as usize;
    let d = t & ((!0u64) >> 17);
    let c1_term = (LOG_B[i].c1 as i64).wrapping_mul((d >> 16) as i64) as u64;
    let j = t
        .wrapping_add((LOG_B[i].c0 as u64) << 33)
        .wrapping_add(c1_term)
        >> (52 - 10);
    t |= (0x3ffu64) << 52;
    let i1 = (j >> 5) as usize;
    let i2 = (j & 0x1f) as usize;
    let r = LOG_R1[i1] * LOG_R2[i2];
    let tf = asdouble(t);
    let o = r * tf;
    let dxl = fma_internal(r, tf, -o);
    let dxh = o - 1.0;
    let dx = fma_internal(r, tf, -1.0);
    let dx2 = dx * dx;
    let f = dx2 * ((LOG_C[0] + dx * LOG_C[1]) + dx2 * (LOG_C[2] + dx * LOG_C[3]));
    let lt = (LOG_L1[i1][1] + LOG_L2[i2][1]) + ed * f64::from_bits(0x3fe62e42fef80000);
    let lh = lt + dxh;
    let mut ll = (lt - lh) + dxh;
    ll += (LOG_L1[i1][0] + LOG_L2[i2][0]) + f64::from_bits(0x3db1cf79abc9e3b4) * ed + dxl;
    ll += f;
    *l = ll;
    lh
}

fn as_logd_accurate(x: f64, l: &mut f64, l2: &mut f64) -> f64 {
    let mut t = (x).to_bits();
    let mut ex = (t >> 52) as i32;
    if ex == 0 {
        let k = t.leading_zeros() as i32;
        t <<= (k - 11) as u32;
        ex -= k - 12;
    }
    let e = ex - 0x3ff;
    t &= (!0u64) >> 12;
    let ed = e as f64;
    let i = (t >> (52 - 5)) as usize;
    let d = t & ((!0u64) >> 17);
    let c1_term = (LOG_B[i].c1 as i64).wrapping_mul((d >> 16) as i64) as u64;
    let j = t
        .wrapping_add((LOG_B[i].c0 as u64) << 33)
        .wrapping_add(c1_term)
        >> (52 - 10);
    t |= (0x3ffu64) << 52;
    let i1 = (j >> 5) as usize;
    let i2 = (j & 0x1f) as usize;
    let r = LOG_R1[i1] * LOG_R2[i2];
    let tf = asdouble(t);
    let o = r * tf;
    let dxl = fma_internal(r, tf, -o);
    let mut dxh = o - 1.0;

    let mut dxl2 = 0.0;
    dxh = fasttwosum(dxh, dxl, &mut dxl2);
    let mut fl = dxh * (LOG_C_ACC[6][0] + dxh * (LOG_C_ACC[7][0] + dxh * LOG_C_ACC[8][0]));
    let mut fh = polydd2(dxh, dxl2, 6, &LOG_C_ACC, &mut fl);
    fh = muldd2(dxh, dxl2, fh, fl, &mut fl);

    let s2 = LOG_H1[i1][2] + LOG_H2[i2][2];
    let s1 = LOG_H1[i1][1] + LOG_H2[i2][1];
    let s0 = LOG_H1[i1][0] + LOG_H2[i2][0];

    let mut l0 = f64::from_bits(0x3fe62e42fefa3800) * ed;
    let mut l1 = f64::from_bits(0x3d2ef35793c76000) * ed;
    let mut l2v = f64::from_bits(0x3a8cc01f97b57a08) * ed;

    l0 += s2;
    l1 = sumdd(l1, l2v, s1, s0, &mut l2v);
    l1 = sumdd(l1, l2v, fh, fl, &mut l2v);

    l0 = fasttwosum(l0, l1, &mut l1);
    l1 = fasttwosum(l1, l2v, &mut l2v);

    *l = l1;
    *l2 = l2v;
    l0
}
#[inline(never)]
fn as_sinpipid(mut x: f64, l: &mut f64) -> f64 {
    x -= 0.5;
    let ax = x.abs();
    let sx = ax * 128.0;
    let ix = roundeven_finite(sx);
    let ky = ix as i32;
    let kx = 64 - ky;
    if kx < 2 {
        let z = 0.5 - ax;
        let z2 = z * z;
        let z2l = fma_internal(z, z, -z2);
        let mut fl = z2 * (SINPID_CL_NEAR[0] + z2 * (SINPID_CL_NEAR[1] + z2 * SINPID_CL_NEAR[2]));
        let mut e = 0.0;
        let mut fh = fasttwosum(SINPID_C_NEAR[0], fl, &mut fl);
        fl += SINPID_C_NEAR[1];
        fh = muldd2(z2, z2l, fh, fl, &mut fl);
        fh = mulddd2(z, fh, fl, &mut fl);
        fh = fasttwosum(z, fh, &mut e);
        fl += e;
        *l = fl;
        return fh;
    }

    let d = ix - sx;
    let d2 = d * d;

    let mut sh = STPI[kx as usize][1];
    let sl = STPI[kx as usize][0];
    let mut ch = STPI[ky as usize][1];
    let mut cl = STPI[ky as usize][0];

    let c0 = f64::from_bits(0xbbd692b66e3cf6e8);
    let s0 = f64::from_bits(0x3c31a624b88c9448);

    let p = d2 * (SINPID_C_COEFF[1] + d2 * (SINPID_C_COEFF[2] + d2 * SINPID_C_COEFF[3]));
    let q = d2 * (SINPID_S_COEFF[1] + d2 * (SINPID_S_COEFF[2] + d2 * SINPID_S_COEFF[3]));

    let mut ql = 0.0;
    let qh = fasttwosum(SINPID_S_COEFF[0], q, &mut ql);
    ql += s0;
    ch = muldd2(qh, ql, ch, cl, &mut cl);
    let mut tl = 0.0;
    let mut th = fasttwosum(SINPID_C_COEFF[0], p, &mut tl);
    tl += c0;
    th = mulddd2(d, th, tl, &mut tl);
    let mut pl = 0.0;
    let ph = muldd2(th, tl, sh, sl, &mut pl);
    ch = fastsum(ch, cl, ph, pl, &mut cl);
    ch = mulddd2(d, ch, cl, &mut cl);
    sh = fastsum(sh, sl, ch, cl, l);
    sh
}

#[inline(never)]
fn as_sinpipid_accurate(mut x: f64, l: &mut f64) -> f64 {
    x -= 0.5;
    x = x.abs();
    x *= 128.0;
    let ix = roundeven_finite(x);
    let d = ix - x;
    let ky = ix as i32;
    let kx = 64 - ky;

    let mut sh = STPI[kx as usize][1];
    let sl = STPI[kx as usize][0];
    let mut ch = STPI[ky as usize][1];
    let mut cl = STPI[ky as usize][0];

    let d2h = d * d;
    let d2l = fma_internal(d, d, -d2h);
    let mut pl = 0.0;
    let mut ph = polydd2(d2h, d2l, SINPID_ACC_C.len(), &SINPID_ACC_C, &mut pl);
    let mut ql = 0.0;
    let mut qh = polydd2(d2h, d2l, SINPID_ACC_S.len(), &SINPID_ACC_S, &mut ql);

    ph = mulddd2(d, ph, pl, &mut pl);
    ph = muldd2(sh, sl, ph, pl, &mut pl);
    qh = muldd2(ch, cl, qh, ql, &mut ql);

    ch = fastsum(qh, ql, ph, pl, &mut cl);
    ch = mulddd2(d, ch, cl, &mut cl);
    sh = fastsum(sh, sl, ch, cl, l);
    sh
}

const LGAMMA_ASYM_C_BIG: [[f64; 2]; 8] = [
    [
        f64::from_bits(0x3fdacfe390c97d69),
        f64::from_bits(0x3c734acf208a22c4),
    ],
    [
        f64::from_bits(0x3fb5555555555555),
        f64::from_bits(0x3c531799ffbcdddb),
    ],
    [
        f64::from_bits(0xbf66c16c16c165a9),
        f64::from_bits(0x3c01eefaee02f690),
    ],
    [
        f64::from_bits(0x3f4a01a019ada522),
        f64::from_bits(0xbbd4d52971deb155),
    ],
    [
        f64::from_bits(0xbf4381377e3a546d),
        f64::from_bits(0x3befd1b354a8db62),
    ],
    [
        f64::from_bits(0x3f4b9486dc1c9886),
        f64::from_bits(0xbbe2dac4b8cca031),
    ],
    [
        f64::from_bits(0xbf5f3ecd8799f337),
        f64::from_bits(0x3bfda5dd745e3963),
    ],
    [
        f64::from_bits(0x3f76d399e5618390),
        f64::from_bits(0x3c115e3000de141a),
    ],
];

const LGAMMA_ASYM_C_SMALL: [[f64; 2]; 13] = [
    [
        f64::from_bits(0x3fdacfe390c97d69),
        f64::from_bits(0x3c7f06a157d44d5b),
    ],
    [
        f64::from_bits(0x3fb5555555555541),
        f64::from_bits(0x3c59d5fc10df4161),
    ],
    [
        f64::from_bits(0xbf66c16c16bfb733),
        f64::from_bits(0xbbf557d8fba9e97a),
    ],
    [
        f64::from_bits(0x3f4a01a01651819c),
        f64::from_bits(0xbbedd3c0f402122a),
    ],
    [
        f64::from_bits(0xbf438136b229bfb4),
        f64::from_bits(0xbbc879990edddc5f),
    ],
    [
        f64::from_bits(0x3f4b94c0472d00a0),
        f64::from_bits(0x3be215a15f7d9289),
    ],
    [
        f64::from_bits(0xbf5f619a122c3918),
        f64::from_bits(0x3bf13405abdba76d),
    ],
    [
        f64::from_bits(0x3f79edef47081644),
        f64::from_bits(0x3c1d9d833b12b9b0),
    ],
    [
        f64::from_bits(0xbf9bfc20185bf7cc),
        f64::from_bits(0xbc38aa555605e3b1),
    ],
    [
        f64::from_bits(0x3fc0e832a9372330),
        f64::from_bits(0x3c5871cbde1ab342),
    ],
    [
        f64::from_bits(0xbfe2beb46518ed4a),
        f64::from_bits(0xbc50298c44c99cee),
    ],
    [
        f64::from_bits(0x3ffe5717107e0999),
        f64::from_bits(0x3c75bdfe7ac38f81),
    ],
    [
        f64::from_bits(0xc0090c04fbd840a6),
        f64::from_bits(0xbc95d2fbfe47e148),
    ],
];

#[inline(never)]
fn as_lgamma_asym(xh: f64, xl: &mut f64) -> f64 {
    let zh = 1.0 / xh;
    let dz = *xl * zh;
    let zl = (fma_internal(zh, -xh, 1.0) - dz) * zh;
    let mut ll = 0.0;
    let mut lh = as_logd(xh, &mut ll);
    ll += dz;
    lh = muldd2(xh - 0.5, *xl, lh - 1.0, ll, &mut ll);
    let mut z2l = 0.0;
    let z2h = muldd2(zh, zl, zh, zl, &mut z2l);
    let x2 = z2h * z2h;
    let mut fh;
    let mut fl;
    if xh > 11.5 {
        let c = &LGAMMA_ASYM_C_BIG;
        lh = fastsum(lh, ll, c[0][0], c[0][1], &mut ll);
        let q0 = c[2][0] + z2h * c[3][0];
        let q2 = c[4][0] + z2h * c[5][0];
        let q4 = c[6][0] + z2h * c[7][0];
        fl = z2h * (q0 + x2 * (q2 + x2 * q4));
        fh = polydd2(z2h, z2l, 1, &c[1..], &mut fl);
    } else {
        let c = &LGAMMA_ASYM_C_SMALL;
        lh = fastsum(lh, ll, c[0][0], c[0][1], &mut ll);
        let x4 = x2 * x2;
        let q0 = c[3][0] + z2h * c[4][0];
        let q2 = c[5][0] + z2h * c[6][0];
        let q4 = c[7][0] + z2h * c[8][0];
        let q6 = c[9][0] + z2h * c[10][0];
        let q8 = c[11][0] + z2h * c[12][0];
        let q4 = q4 + x2 * (q6 + x2 * q8);
        let q0 = q0 + x2 * q2 + x4 * q4;
        fl = z2h * q0;
        fh = polydd2(z2h, z2l, 2, &c[1..], &mut fl);
    }
    fh = muldd2(zh, zl, fh, fl, &mut fl);
    fastsum(lh, ll, fh, fl, xl)
}

#[inline(never)]
fn as_lgamma_asym_accurate(xh: f64, xl: &mut f64, e: &mut f64) -> f64 {
    let mut l2 = 0.0;
    let mut l1 = 0.0;
    let mut l0 = as_logd_accurate(xh, &mut l1, &mut l2);
    let mut l0x: f64;
    let mut l1x: f64;
    let mut l2x: f64;

    if xh < f64::from_bits(0x47b0000000000000) {
        // 0x1p120
        let zh = 1.0 / xh;
        let dz = *xl * zh;
        let zl = (fma_internal(zh, -xh, 1.0) - dz) * zh;
        if *xl != 0.0 {
            let mut dl2 = 0.0;
            let dl1 = mulddd2(*xl, zh, fma_internal(zh, -xh, 1.0) * zh, &mut dl2);
            dl2 -= dl1 * dl1 * 0.5;
            l1 = sumdd(l1, l2, dl1, dl2, &mut l2);
        }

        let (wh, wl) = if ((xh).to_bits() >> 52) > (0x3ff + 51) {
            (xh, *xl - 0.5)
        } else {
            (xh - 0.5, *xl)
        };
        l0 -= 1.0;

        l0x = l0 * wh;
        let l0xl = fma_internal(l0, wh, -l0x);
        l1x = l1 * wh;
        let l1xl = fma_internal(l1, wh, -l1x);
        l2x = l2 * wh;

        l1x = sumdd(l1x, l2x, l0xl, l1xl, &mut l2x);
        l1x = sumdd(l1x, l2x, l0 * wl, l1 * wl, &mut l2x);

        let mut z2l = 0.0;
        let z2h = muldd2(zh, zl, zh, zl, &mut z2l);
        let mut fh: f64;
        let mut fl: f64;

        if xh >= 48.0 {
            l1x = sumdd(
                l1x,
                l2x,
                LGAMMA_ASYM_ACC_C1[0][0],
                LGAMMA_ASYM_ACC_C1[0][1],
                &mut l2x,
            );
            fl = 0.0;
            fh = polydd2(
                z2h,
                z2l,
                LGAMMA_ASYM_ACC_C1.len() - 1,
                &LGAMMA_ASYM_ACC_C1[1..],
                &mut fl,
            );
        } else if xh >= 14.5 {
            l1x = sumdd(
                l1x,
                l2x,
                LGAMMA_ASYM_ACC_C2[0][0],
                LGAMMA_ASYM_ACC_C2[0][1],
                &mut l2x,
            );
            fl = 0.0;
            fh = polydd2(
                z2h,
                z2l,
                LGAMMA_ASYM_ACC_C2.len() - 1,
                &LGAMMA_ASYM_ACC_C2[1..],
                &mut fl,
            );
        } else {
            l1x = sumdd(
                l1x,
                l2x,
                LGAMMA_ASYM_ACC_C3[0][0],
                LGAMMA_ASYM_ACC_C3[0][1],
                &mut l2x,
            );
            fl = 0.0;
            fh = polydd2(
                z2h,
                z2l,
                LGAMMA_ASYM_ACC_C3.len() - 1,
                &LGAMMA_ASYM_ACC_C3[1..],
                &mut fl,
            );
        }
        fh = muldd2(zh, zl, fh, fl, &mut fl);
        l1x = sumdd(l1x, l2x, fh, fl, &mut l2x);
        l0x = fasttwosum(l0x, l1x, &mut l1x);
        l1x = fasttwosum(l1x, l2x, &mut l2x);
    } else {
        let wl = *xl - 0.5;
        l0 -= 1.0;
        l0x = l0 * xh;
        let l0xl = fma_internal(l0, xh, -l0x);
        l1x = l1 * xh;
        let l1xl = fma_internal(l1, xh, -l1x);
        l2x = l2 * xh;
        l1x = sumdd(l1x, l2x, l0xl, l1xl, &mut l2x);
        l1x = sumdd(l1x, l2x, l0 * wl, l1 * wl, &mut l2x);
    }

    *xl = l1x;
    *e = l2x;
    l0x
}

#[inline(always)]
fn copysign(x: f64, y: f64) -> f64 {
    let xb = x.to_bits();
    let yb = y.to_bits();
    f64::from_bits((xb & 0x7fff_ffff_ffff_ffff) | (yb & 0x8000_0000_0000_0000))
}

#[inline(always)]
fn exp_dd_parts(h: f64, l: f64) -> (f64, f64) {
    let eh = exp(h);
    if !eh.is_finite() || l == 0.0 {
        return (eh, 0.0);
    }
    let em1 = expm1(l);
    let ch = 1.0 + em1;
    let cl = em1 - (ch - 1.0);
    let mut pl = 0.0;
    let ph = muldd2(eh, 0.0, ch, cl, &mut pl);
    (ph, pl)
}

#[inline(always)]
fn exp_dd(h: f64, l: f64) -> f64 {
    let (ph, pl) = exp_dd_parts(h, l);
    ph + pl
}

// ---- tgamma fast-path coefficients (core-math) ----

const TGAMMA_SMALL_CC: [[f64; 2]; 28] = [
    [
        f64::from_bits(0xbfe2788cfc6fb619),
        f64::from_bits(0x3c56cb90701fbfba),
    ],
    [
        f64::from_bits(0x3fefa658c23b1578),
        f64::from_bits(0x3c8dd92b465a81dd),
    ],
    [
        f64::from_bits(0xbfed0a118f324b63),
        f64::from_bits(0x3c53a4f48373c073),
    ],
    [
        f64::from_bits(0x3fef6a51055096b5),
        f64::from_bits(0x3c7fabe4f73da654),
    ],
    [
        f64::from_bits(0xbfef6c80ec38b67b),
        f64::from_bits(0x3c8c9fc797e7567e),
    ],
    [
        f64::from_bits(0x3fefc7e0a6eb310b),
        f64::from_bits(0xbc5042340fb21e2f),
    ],
    [
        f64::from_bits(0xbfefdf3f157b7a39),
        f64::from_bits(0xbc86fd12a61ae3a9),
    ],
    [
        f64::from_bits(0x3feff07b5a17ff6c),
        f64::from_bits(0xbc89b8f7ad70e4bc),
    ],
    [
        f64::from_bits(0xbfeff803d68a0bd4),
        f64::from_bits(0x3c4ed4f2b2dcc156),
    ],
    [
        f64::from_bits(0x3feffc0841d585a3),
        f64::from_bits(0xbc87089ffbe760f0),
    ],
    [
        f64::from_bits(0xbfeffe018c484f48),
        f64::from_bits(0x3c8f8a632e2ff912),
    ],
    [
        f64::from_bits(0x3fefff00b768f1c4),
        f64::from_bits(0x3c750de70bb4e28b),
    ],
    [
        f64::from_bits(0xbfefff8035584df4),
        f64::from_bits(0xbc8aae8f6d8b868d),
    ],
    [
        f64::from_bits(0x3fefffc012f95019),
        f64::from_bits(0x3c8d26498825213d),
    ],
    [
        f64::from_bits(0xbfefffe0062afaf7),
        f64::from_bits(0x3c8ef9359641ed4b),
    ],
    [
        f64::from_bits(0x3feffff002146fec),
        f64::from_bits(0x3c81ce6548eaa3c0),
    ],
    [
        f64::from_bits(0xbfeffff800af4ca8),
        f64::from_bits(0x3c767121c2223cb7),
    ],
    [
        f64::from_bits(0x3feffffc0037a50d),
        f64::from_bits(0xbc67c2694e8ff6af),
    ],
    [
        f64::from_bits(0xbfeffffe00676ecc),
        f64::from_bits(0xbc82563151dfa334),
    ],
    [
        f64::from_bits(0x3fefffff00ac53bf),
        f64::from_bits(0x3c83effd975c6fec),
    ],
    [
        f64::from_bits(0xbfefffff72b5b5d1),
        f64::from_bits(0x3c8fb25bd33b12e8),
    ],
    [
        f64::from_bits(0x3fefffffa828c9f3),
        f64::from_bits(0x3c85fde8ff2fb275),
    ],
    [
        f64::from_bits(0xbff00000bc9c338d),
        f64::from_bits(0xbc77ced37913c915),
    ],
    [
        f64::from_bits(0x3ff000014750b54c),
        f64::from_bits(0x3c9980564dfb043c),
    ],
    [
        f64::from_bits(0xbfefffdabf3f1dfb),
        f64::from_bits(0xbc8e34acb3149e71),
    ],
    [
        f64::from_bits(0x3fefffc7cf7439c3),
        f64::from_bits(0xbc82599c887c744a),
    ],
    [
        f64::from_bits(0xbff001452e7ff0b7),
        f64::from_bits(0xbc7acc4a235f0de4),
    ],
    [
        f64::from_bits(0x3ff001c6c5b18192),
        f64::from_bits(0x3c935874aa96ce14),
    ],
];

const TGAMMA_SMALL_C: [f64; 8] = [
    f64::from_bits(0xbfefdf5126e6a83b),
    f64::from_bits(0x3fefd56c1b531709),
    f64::from_bits(0xbff0956741759e58),
    f64::from_bits(0x3ff0b6167d27c5c4),
    f64::from_bits(0xbfe8e69e55b6dca0),
    f64::from_bits(0x3fe7e124660c827d),
    f64::from_bits(0xbffc797b42bd6745),
    f64::from_bits(0x3ffd688bb8db046c),
];

const TGAMMA_CC: [[f64; 2]; 18] = [
    [
        f64::from_bits(0x400a96390899a074),
        f64::from_bits(0xbc56e95430fab070),
    ],
    [
        f64::from_bits(0x400d545472146024),
        f64::from_bits(0x3c7c07f9774e12b3),
    ],
    [
        f64::from_bits(0x400491ad1cb98836),
        f64::from_bits(0x3ca51e26c4cfd792),
    ],
    [
        f64::from_bits(0x3ff4a0b6a8230929),
        f64::from_bits(0x3c9c1c6993b10594),
    ],
    [
        f64::from_bits(0x3fe0e5d232b95859),
        f64::from_bits(0x3c7d4248748dd78b),
    ],
    [
        f64::from_bits(0x3fc71d1672129fee),
        f64::from_bits(0x3c43b47c61245ee6),
    ],
    [
        f64::from_bits(0x3fabd2afde7e4816),
        f64::from_bits(0xbc325466b734902d),
    ],
    [
        f64::from_bits(0x3f8d8376e1031a16),
        f64::from_bits(0x3c22cd76af7fbb20),
    ],
    [
        f64::from_bits(0x3f6c9e94992c88c1),
        f64::from_bits(0x3bf5d7be78c93d16),
    ],
    [
        f64::from_bits(0x3f490ba7276a0c19),
        f64::from_bits(0xbbd6cad258076bb3),
    ],
    [
        f64::from_bits(0x3f249cfed9d63c8b),
        f64::from_bits(0x3b50a8ada0cff18d),
    ],
    [
        f64::from_bits(0x3efec018849c245b),
        f64::from_bits(0x3b9cea7c4e5e9d4f),
    ],
    [
        f64::from_bits(0x3ed65e5a18d31c17),
        f64::from_bits(0x3b712fc2f27069ec),
    ],
    [
        f64::from_bits(0x3eaca1890add8727),
        f64::from_bits(0x3b469c0fe53eb0fa),
    ],
    [
        f64::from_bits(0x3e8378b3b91f9033),
        f64::from_bits(0xbb1d62590f524392),
    ],
    [
        f64::from_bits(0x3e5432cdb3640fca),
        f64::from_bits(0xbae33987f0b3b6b6),
    ],
    [
        f64::from_bits(0x3e2f239fc9cf2155),
        f64::from_bits(0xbac8a95d04bfb2e4),
    ],
    [
        f64::from_bits(0x3dee3ea4e1366932),
        f64::from_bits(0xba25c950f5465458),
    ],
];

const TGAMMA_C: [f64; 18] = [
    f64::from_bits(0x400a96390899a074),
    f64::from_bits(0x400d545472146024),
    f64::from_bits(0x400491ad1cb98836),
    f64::from_bits(0x3ff4a0b6a8230929),
    f64::from_bits(0x3fe0e5d232b95859),
    f64::from_bits(0x3fc71d1672129fee),
    f64::from_bits(0x3fabd2afde7e4816),
    f64::from_bits(0x3f8d8376e1031a16),
    f64::from_bits(0x3f6c9e94992c88c1),
    f64::from_bits(0x3f490ba7276a0c19),
    f64::from_bits(0x3f249cfed9d63c8b),
    f64::from_bits(0x3efec018849c245b),
    f64::from_bits(0x3ed65e5a18d31c17),
    f64::from_bits(0x3eaca1890add8727),
    f64::from_bits(0x3e8378b3b91f9033),
    f64::from_bits(0x3e5432cdb3640fca),
    f64::from_bits(0x3e2f239fc9cf2155),
    f64::from_bits(0x3dee3ea4e1366932),
];

const TGAMMA_EPS: f64 = f64::from_bits(0x3b92e3b40a0e9b4f);
const TGAMMA_EPS_BIG: f64 = f64::from_bits(0x3ba2e3b40a0e9b4f);
const TGAMMA_EPS_BIG2: f64 = f64::from_bits(0x3b66aad80c11872c);

#[inline(always)]
fn tgamma_fast(x: f64) -> Option<f64> {
    let ax = x.abs();
    if ax < 0.25 {
        let x2 = x * x;
        let x4 = x2 * x2;
        let mut c0 = TGAMMA_SMALL_C[0]
            + x * TGAMMA_SMALL_C[1]
            + x2 * (TGAMMA_SMALL_C[2] + x * TGAMMA_SMALL_C[3]);
        let c4 = TGAMMA_SMALL_C[4]
            + x * TGAMMA_SMALL_C[5]
            + x2 * (TGAMMA_SMALL_C[6] + x * TGAMMA_SMALL_C[7]);
        c0 += x4 * c4;

        let mut cl = x * c0;
        let ch = polyddd(x, TGAMMA_SMALL_CC.len(), &TGAMMA_SMALL_CC, &mut cl);
        let mut fh = 1.0 / x;
        let dh = fma_internal(fh, -x, 1.0);
        let mut fl = dh * fh;
        let mut fll = fma_internal(fl, -x, dh) * fh;
        fl = sumdd(fl, fll, ch, cl, &mut fll);
        fl = twosum(fl, fll, &mut fll);
        fh = fasttwosum(fh, fl, &mut fl);
        fl = fasttwosum(fl, fll, &mut fll);
        let mut et = 0.0;
        fasttwosum(fh, 2.0 * fl, &mut et);
        if et == 0.0 {
            let bump = f64::from_bits(0x3cb0000000000000);
            if copysign(1.0, fl) * copysign(1.0, fll) > 0.0 {
                fl *= 1.0 + bump;
            } else {
                fl *= 1.0 - bump;
            }
        }
        return Some(fh + fl);
    }

    if x >= f64::from_bits(0x406573fae561f648) {
        return Some(f64::INFINITY);
    }

    let fx = floor(x);
    if fx == x {
        if x == 0.0 {
            return Some(f64::INFINITY);
        }
        if x < 0.0 {
            return Some(f64::NAN);
        }
        let n = x as i64;
        if n <= 20 {
            let mut acc: u128 = 1;
            for i in 2..n {
                acc *= i as u128;
            }
            return Some(acc as f64);
        }
        let mut hi = 1.0;
        let mut lo = 0.0;
        let mut x0 = 1.0;
        for _ in 1..n {
            let mut l = 0.0;
            hi = mulddd2(x0, hi, lo, &mut l);
            lo = l;
            x0 += 1.0;
        }
        return Some(hi + lo);
    }

    if x <= -184.0 {
        return None;
    }

    if x < -3.0 {
        let mut ll = 0.0;
        let lh = as_lgamma_asym(fasttwosum(-x, 1.0, &mut ll), &mut ll);
        let (yh, yl) = exp_dd_parts(lh, ll);
        let mut sl = 0.0;
        let sh = sinpi_parts(x, &mut sl);
        let mut pl = 0.0;
        let ph = muldd2(sh, sl, yh, yl, &mut pl);
        let rh = 1.0 / ph;
        let rl = (fma_internal(rh, -ph, 1.0) - pl * rh) * rh;
        let eps =
            rh * (f64::from_bits(0x3bbeb2049057bc61) - x * f64::from_bits(0x3b661019f74442b7));
        let ub = rh + (rl + eps);
        let lb = rh + (rl - eps);
        if ub == lb {
            return Some(ub);
        }
        return None;
    }

    if x > 4.0 {
        let mut ll = 0.0;
        let lh = as_lgamma_asym(x, &mut ll);
        let (yh, yl) = exp_dd_parts(lh, ll);
        let eps = yh * (TGAMMA_EPS_BIG + x * TGAMMA_EPS_BIG2);
        let ub = yh + (yl + eps);
        let lb = yh + (yl - eps);
        if ub == lb {
            return Some(ub);
        }
        return None;
    }

    let z = x;
    let m = z - f64::from_bits(0x400c000000000000);
    let i = roundeven_finite(m);
    let d = z - (i + f64::from_bits(0x400c000000000000));
    let d2 = d * d;
    let d4 = d2 * d2;
    let mut fl = d
        * ((TGAMMA_C[10] + d * TGAMMA_C[11])
            + d2 * (TGAMMA_C[12] + d * TGAMMA_C[13])
            + d4 * ((TGAMMA_C[14] + d * TGAMMA_C[15]) + d2 * (TGAMMA_C[16] + d * TGAMMA_C[17])));
    let mut fh = polyddd(d, 10, &TGAMMA_CC, &mut fl);
    let jm = i.abs() as i32;
    let mut wh = 1.0;
    let mut wl = 0.0;
    let mut xph = z;
    let mut xpl = 0.0;
    if jm != 0 {
        wh = xph;
        for _ in 1..jm {
            let mut l = 0.0;
            if xph.abs() > 1.0 {
                xph = fasttwosum(xph, 1.0, &mut l);
            } else {
                xph = fasttwosum(1.0, xph, &mut l);
            }
            xpl += l;
            wh = muldd2(xph, xpl, wh, wl, &mut wl);
        }
    }
    let rh = 1.0 / wh;
    let rl = (fma_internal(rh, -wh, 1.0) - wl * rh) * rh;
    fh = muldd2(rh, rl, fh, fl, &mut fl);
    let eps = fh * TGAMMA_EPS;
    let ub = fh + (fl + eps);
    let lb = fh + (fl - eps);
    if ub == lb {
        return Some(ub);
    }
    None
}

#[inline(always)]
fn sinpi_parts(x: f64, sl: &mut f64) -> f64 {
    let k = floor(x);
    let phi = x - k;
    let sh = as_sinpipid_accurate(phi, sl);
    if ((k as i64) & 1) != 0 {
        *sl = -*sl;
        return -sh;
    }
    sh
}

#[inline(always)]
fn tgamma_half(x: f64) -> f64 {
    let n = rint(x - 0.5) as i32;
    let mut y = f64::from_bits(0x3ffc5bf891b4ef6a); // sqrt(pi)
    if n >= 0 {
        for i in 0..n {
            y *= i as f64 + 0.5;
        }
    } else {
        for k in 1..=(-n) {
            y /= 0.5 - k as f64;
        }
    }
    y
}

#[inline(always)]
fn tgamma_quarter(x: f64, frac: f64) -> f64 {
    let (hi, lo) = tgamma_quarter_dd(x, frac);
    hi + lo
}

#[inline(always)]
fn tgamma_quarter_dd(x: f64, frac: f64) -> (f64, f64) {
    let base = if frac == 0.25 { 0.25 } else { 0.75 };
    let (mut hi, mut lo) = if frac == 0.25 {
        (
            f64::from_bits(0x400d013fc47eeeea), // gamma(1/4)
            f64::from_bits(0x3c9e6ce29429451b),
        )
    } else {
        (
            f64::from_bits(0x3ff39b4e8b50f62c), // gamma(3/4)
            f64::from_bits(0x3c43d7a9256698c6),
        )
    };
    let n = rint(x - base) as i32;
    if n >= 0 {
        for i in 0..n {
            let mut l = 0.0;
            hi = mulddd2(base + i as f64, hi, lo, &mut l);
            lo = l;
        }
    } else {
        for k in 1..=(-n) {
            let d = base - k as f64;
            let inv = 1.0 / d;
            let inv_err = fma_internal(inv, d, -1.0) * inv;
            let mut l = 0.0;
            hi = muldd2(hi, lo, inv, inv_err, &mut l);
            lo = l;
        }
    }
    (hi, lo)
}

#[inline(always)]
fn eval_lgamma_acc_block(
    sx: f64,
    x0: &[f64; 3],
    c: &[[f64; 2]],
    sc: f64,
    k: usize,
    fh: &mut f64,
    fl: &mut f64,
) {
    let mut zl = 0.0;
    let zh = fasttwosum(x0[0] + sx, x0[1], &mut zl);
    zl += x0[2];
    let sh = zh * sc;
    let sl = zl * sc;
    let n = c.len();
    *fl = sh * polyd(sh, k, &c[n - k..]);
    *fh = polydd2(sh, sl, n - k, c, fl);
    *fh = muldd2(zh, zl, *fh, *fl, fl);
}

#[inline(never)]
#[allow(unused_assignments)]
fn as_lgamma_accurate_dd(mut x: f64) -> (f64, f64) {
    let sx = x;
    let mut fh = 0.0;
    let mut fl = 0.0;
    let mut fll = 0.0;
    x = x.abs();

    if x < f64::from_bits(0x3ab0000000000000) {
        // 0x1p-100
        let mut lll = 0.0;
        let mut ll = 0.0;
        let lh = as_logd_accurate(x, &mut ll, &mut lll);
        fh = -lh;
        fl = -ll;
        fll = -lll;
        fh = fasttwosum(fh, fl, &mut fl);
        fl = fasttwosum(fl, fll, &mut fll);
        let mut e = 0.0;
        fasttwosum(fh, 2.0 * fl, &mut e);
        if e == 0.0
            && (1.0 + f64::from_bits(0x3c90000000000000))
                == (1.0 - f64::from_bits(0x3c90000000000000))
        {
            fl *= 1.0 + copysign(f64::from_bits(0x3cb0000000000000), fl) * copysign(1.0, fll);
        }
    } else if x < f64::from_bits(0x3fd0000000000000) {
        // 0x1p-2
        fh = polydddfst(sx, LGAMMA_C0.len(), &LGAMMA_C0, &mut fl);
        fh = mulddd2(sx, fh, fl, &mut fl);
        let mut lll = 0.0;
        let mut ll = 0.0;
        let lh = as_logd_accurate(x, &mut ll, &mut lll);
        fh = sumdd(fh, fl, -ll, -lll, &mut fl);
        fh = twosum(-lh, fh, &mut fll);
        fl = twosum(fll, fl, &mut fll);
        let mut e = 0.0;
        fasttwosum(fh, 2.0 * fl, &mut e);
        if e == 0.0
            && (1.0 + f64::from_bits(0x3c90000000000000))
                == (1.0 - f64::from_bits(0x3c90000000000000))
        {
            fl *= 1.0 + copysign(f64::from_bits(0x3cb0000000000000), fl) * copysign(1.0, fll);
        }
    } else {
        if (x - 0.5).abs() < f64::from_bits(0x3fd0000000000000) {
            // 0x1p-2
            fh = polydddfst(x - 0.5, LGAMMA_B.len(), &LGAMMA_B, &mut fl);
            if sx > 0.0 {
                let mut lll = 0.0;
                let mut ll = 0.0;
                let lh = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, 0.0, -ll, -lll, &mut fll);
                let mut lh_err = 0.0;
                fh = twosum(fh, -lh, &mut lh_err);
                fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            }
        } else if (x - 2.5).abs() < f64::from_bits(0x3fd0000000000000) {
            fh = polydddfst(x - 2.5, LGAMMA_B.len(), &LGAMMA_B, &mut fl);
            let mut lll = 0.0;
            let mut ll = 0.0;
            let lh = as_logd_accurate(x - 1.0, &mut ll, &mut lll);
            fl = sumdd(fl, 0.0, ll, lll, &mut fll);
            let mut lh_err = 0.0;
            fh = twosum(fh, lh, &mut lh_err);
            fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            if sx < 0.0 {
                let lh2 = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, fll, ll, lll, &mut fll);
                let mut lh_err2 = 0.0;
                fh = twosum(fh, lh2, &mut lh_err2);
                fl = sumdd(fl, fll, lh_err2, 0.0, &mut fll);
            }
        } else if (x - 3.5).abs() < f64::from_bits(0x3fd0000000000000) {
            let mut l2ll = 0.0;
            let mut l2l = 0.0;
            let l2h = as_logd_accurate(x - 2.0, &mut l2l, &mut l2ll);
            let mut l1ll = 0.0;
            let mut l1l = 0.0;
            let l1h = as_logd_accurate(x - 1.0, &mut l1l, &mut l1ll);
            l1l = sumdd(l1l, l1ll, l2l, l2ll, &mut l1ll);
            let mut l2h_err = 0.0;
            let l1h2 = fasttwosum(l1h, l2h, &mut l2h_err);
            l1l = sumdd(l1l, l1ll, l2h_err, 0.0, &mut l1ll);
            fh = polydddfst(x - 3.5, LGAMMA_B.len(), &LGAMMA_B, &mut fl);
            fl = sumdd(fl, 0.0, l1l, l1ll, &mut fll);
            let mut l1h_err = 0.0;
            fh = twosum(fh, l1h2, &mut l1h_err);
            fl = sumdd(fl, fll, l1h_err, 0.0, &mut fll);
            if sx < 0.0 {
                let mut lll = 0.0;
                let mut ll = 0.0;
                let lh = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, fll, ll, lll, &mut fll);
                let mut lh_err = 0.0;
                fh = twosum(fh, lh, &mut lh_err);
                fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            }
        } else if (x - 1.0).abs() < f64::from_bits(0x3fd0000000000000) {
            fh = polydddfst(x - 1.0, LGAMMA_C0.len(), &LGAMMA_C0, &mut fl);
            fh = mulddd2(x - 1.0, fh, fl, &mut fl);
            if sx < 0.0 {
                let mut lll = 0.0;
                let mut ll = 0.0;
                let lh = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, 0.0, ll, lll, &mut fll);
                let mut lh_err = 0.0;
                fh = twosum(fh, lh, &mut lh_err);
                fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            }
        } else if (x - 1.5).abs() < f64::from_bits(0x3fd0000000000000) {
            fh = polydddfst(x - 1.5, LGAMMA_B.len(), &LGAMMA_B, &mut fl);
            if sx < 0.0 {
                let mut lll = 0.0;
                let mut ll = 0.0;
                let lh = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, 0.0, ll, lll, &mut fll);
                let mut lh_err = 0.0;
                fh = twosum(fh, lh, &mut lh_err);
                fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            }
        } else if (x - 2.0).abs() < f64::from_bits(0x3fd0000000000000) {
            let mut lll = 0.0;
            let mut ll = 0.0;
            let lh = as_logd_accurate(x - 1.0, &mut ll, &mut lll);
            fh = polydddfst(x - 2.0, LGAMMA_C0.len(), &LGAMMA_C0, &mut fl);
            fh = mulddd2(x - 2.0, fh, fl, &mut fl);
            fl = sumdd(fl, 0.0, ll, lll, &mut fll);
            let mut lh_err = 0.0;
            fh = twosum(fh, lh, &mut lh_err);
            fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            if sx < 0.0 {
                let lh2 = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, fll, ll, lll, &mut fll);
                let mut lh_err2 = 0.0;
                fh = twosum(fh, lh2, &mut lh_err2);
                fl = sumdd(fl, fll, lh_err2, 0.0, &mut fll);
            }
        } else if (x - 3.0).abs() < f64::from_bits(0x3fd0000000000000) {
            let mut l2ll = 0.0;
            let mut l2l = 0.0;
            let l2h = as_logd_accurate(x - 2.0, &mut l2l, &mut l2ll);
            let mut l1ll = 0.0;
            let mut l1l = 0.0;
            let l1h = as_logd_accurate(x - 1.0, &mut l1l, &mut l1ll);
            l1l = sumdd(l1l, l1ll, l2l, l2ll, &mut l1ll);
            let mut l2h_err = 0.0;
            let l1h2 = fasttwosum(l1h, l2h, &mut l2h_err);
            l1l = sumdd(l1l, l1ll, l2h_err, 0.0, &mut l1ll);
            fh = polydddfst(x - 3.0, LGAMMA_C0.len(), &LGAMMA_C0, &mut fl);
            fh = mulddd2(x - 3.0, fh, fl, &mut fl);
            fl = sumdd(fl, 0.0, l1l, l1ll, &mut fll);
            let mut l1h_err = 0.0;
            fh = twosum(fh, l1h2, &mut l1h_err);
            fl = sumdd(fl, fll, l1h_err, 0.0, &mut fll);
            if sx < 0.0 {
                let mut lll = 0.0;
                let mut ll = 0.0;
                let lh = as_logd_accurate(x, &mut ll, &mut lll);
                fl = sumdd(fl, fll, ll, lll, &mut fll);
                let mut lh_err = 0.0;
                fh = twosum(fh, lh, &mut lh_err);
                fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            }
        } else {
            if sx < 0.0 {
                let mut tmp = 0.0;
                x = fasttwosum(x, 1.0, &mut tmp);
                fl = tmp;
            }
            fh = as_lgamma_asym_accurate(x, &mut fl, &mut fll);
        }

        if sx < 0.0 {
            let phi = if sx < -0.5 { sx - floor(sx) } else { -sx };
            let mut sl = 0.0;
            let sh = as_sinpipid_accurate(phi, &mut sl);
            let mut lll = 0.0;
            let mut ll = 0.0;
            let lh = as_logd_accurate(sh, &mut ll, &mut lll);
            ll += sl / sh + lll;
            fl = sumdd(fl, fll, ll, 0.0, &mut fll);
            let mut lh_err = 0.0;
            fh = twosum(fh, lh, &mut lh_err);
            fl = sumdd(fl, fll, lh_err, 0.0, &mut fll);
            fh = -fh;
            fl = -fl;
            fll = -fll;
        }
        fh = fasttwosum(fh, fl, &mut fl);
        fl = fasttwosum(fl, fll, &mut fll);
        fh = fasttwosum(fh, fl, &mut fl);
        fl = fasttwosum(fl, fll, &mut fll);
        let mut e = 0.0;
        fasttwosum(fh, 2.0 * fl, &mut e);
        if e == 0.0 {
            let dfl = 1.0
                + copysign(f64::from_bits(0x3e50000000000000), fl)
                    * copysign(f64::from_bits(0x3e50000000000000), fll);
            fl *= dfl;
        }
    }

    const SX_BND: f64 = f64::from_bits(0xc004e147ae147ae1);

    if fh.abs() < f64::from_bits(0x3fd8000000000000) {
        if fh.abs() < f64::from_bits(0x3fb7400000000000) && sx > SX_BND && sx < -2.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_0,
                &LGAMMA_ACC_C_0,
                f64::from_bits(0x4020000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+3
        } else if fh.abs() < f64::from_bits(0x3fb1680000000000) && sx > -3.0 && sx < SX_BND {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_1,
                &LGAMMA_ACC_C_1,
                f64::from_bits(0x4030000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+4
        } else if fh.abs() < f64::from_bits(0x3fb2d40000000000) && sx > -3.5 && sx < -3.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_2,
                &LGAMMA_ACC_C_2,
                f64::from_bits(0x4050000000000000),
                7,
                &mut fh,
                &mut fl,
            ); // 0x1p+6
        } else if fh.abs() < f64::from_bits(0x3faef00000000000) && sx > -4.0 && sx < -3.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_3,
                &LGAMMA_ACC_C_3,
                f64::from_bits(0x4070000000000000),
                7,
                &mut fh,
                &mut fl,
            ); // 0x1p+8
        } else if fh.abs() < f64::from_bits(0x3fc0000000000000) && sx > -4.5 && sx < -4.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_4,
                &LGAMMA_ACC_C_4,
                f64::from_bits(0x4060000000000000),
                7,
                &mut fh,
                &mut fl,
            ); // 0x1p+7
        } else if fh.abs() < f64::from_bits(0x3fb0800000000000) && sx > -5.0 && sx < -4.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_5,
                &LGAMMA_ACC_C_5,
                f64::from_bits(0x4090000000000000),
                7,
                &mut fh,
                &mut fl,
            ); // 0x1p+10
        } else if fh.abs() < f64::from_bits(0x3fb1000000000000) && sx > -5.5 && sx < -5.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_6,
                &LGAMMA_ACC_C_6,
                f64::from_bits(0x4090000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+10
        } else if fh.abs() < f64::from_bits(0x3fc1300000000000) && sx > -6.0 && sx < -5.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_7,
                &LGAMMA_ACC_C_7,
                f64::from_bits(0x40b0000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+12
        } else if fh.abs() < f64::from_bits(0x3fc1500000000000) && sx > -6.5 && sx < -6.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_8,
                &LGAMMA_ACC_C_8,
                f64::from_bits(0x40b0000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+12
        } else if fh.abs() < f64::from_bits(0x3fc3400000000000) && sx > -7.0 && sx < -6.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_9,
                &LGAMMA_ACC_C_9,
                f64::from_bits(0x40d0000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+14
        } else if fh.abs() < f64::from_bits(0x3fc1450000000000) && sx > -7.5 && sx < -7.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_10,
                &LGAMMA_ACC_C_10,
                f64::from_bits(0x40e0000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+15
        } else if fh.abs() < f64::from_bits(0x3fc76b0000000000) && sx > -8.0 && sx < -7.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_11,
                &LGAMMA_ACC_C_11,
                f64::from_bits(0x4100000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+17
        } else if fh.abs() < f64::from_bits(0x3fc76c0000000000)
            && sx > f64::from_bits(0xc0210f5c28f5c28f)
            && sx < -8.0
        {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_12,
                &LGAMMA_ACC_C_12,
                f64::from_bits(0x4100000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+17
        } else if fh.abs() < f64::from_bits(0x3fc9900000000000) && sx > -9.0 && sx < -8.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_13,
                &LGAMMA_ACC_C_13,
                f64::from_bits(0x4130000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+20
        } else if fh.abs() < f64::from_bits(0x3fc9900000000000) && sx > -9.5 && sx < -9.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_14,
                &LGAMMA_ACC_C_14,
                f64::from_bits(0x4130000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+20
        } else if fh.abs() < f64::from_bits(0x3fc76b0000000000) && sx > -10.0 && sx < -9.5 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_15,
                &LGAMMA_ACC_C_15,
                f64::from_bits(0x4170000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+24
        } else if fh.abs() < f64::from_bits(0x3fc76b0000000000) && sx > -10.5 && sx < -10.0 {
            eval_lgamma_acc_block(
                sx,
                &LGAMMA_ACC_X0_16,
                &LGAMMA_ACC_C_16,
                f64::from_bits(0x4170000000000000),
                6,
                &mut fh,
                &mut fl,
            ); // 0x1p+24
        }
    }

    let ft = (fl).to_bits().wrapping_add(2) & ((!0u64) >> 12);
    if ft <= 2 {
        return (as_lgamma_database(sx, fh + fl), 0.0);
    }
    (fh, fl)
}
#[inline(never)]
fn as_lgamma_accurate(x: f64) -> f64 {
    let (fh, fl) = as_lgamma_accurate_dd(x);
    fh + fl
}
#[inline(never)]
#[allow(unused_assignments)]
fn ieee754_lgamma_r(x: f64, signgamp: &mut i32) -> f64 {
    let t = (x).to_bits();
    let nx = t << 1;
    if nx >= 0xfeaea9b24f16a34c {
        // |x| >= 0x1.006df1bfac84ep+1015
        *signgamp = 1;
        if t == 0x7f5754d9278b51a6 {
            return f64::from_bits(0x7feffffffffffffe) - f64::from_bits(0x7c80000000000000); // 0x1.ffffffffffffep+1023 - 0x1p+969
        }
        if t == 0x7f5754d9278b51a7 {
            return f64::from_bits(0x7fefffffffffffff) - f64::from_bits(0x7c80000000000000); // 0x1.fffffffffffffp+1023 - 0x1p+969
        }
        if nx >= ((0x7ffu64) << 53) {
            if nx == ((0x7ffu64) << 53) {
                return x.abs();
            }
            return x + x;
        }
        if (t >> 63) != 0 {
            return f64::INFINITY;
        }
        return f64::from_bits(0x7fe0000000000000) * f64::from_bits(0x7fe0000000000000);
    }

    let fx = floor(x);
    if fx == x {
        if x <= 0.0 {
            *signgamp = 1 - 2 * ((t >> 63) as i32);
            return f64::INFINITY;
        }
        if x == 1.0 || x == 2.0 {
            *signgamp = 1;
            return 0.0;
        }
    }

    let mut au = (nx >> 38) as u32;
    let mut fh = 0.0;
    let mut fl = 0.0;
    let mut eps = 0.0;

    if au < UBRD[0] {
        // |x| < 0.5
        *signgamp = 1 - 2 * ((t >> 63) as i32);
        let mut ll = 0.0;
        let lh = as_logd(x.abs(), &mut ll);
        if au < 0x1da0000 {
            // |x| < 0x1p-75
            fh = -lh;
            fl = -ll;
            eps = 1.5e-22;
        } else if au < 0x1fd0000 {
            // |x| < 0.03125
            let z = x;
            let z2 = z * z;
            let z4 = z2 * z2;
            let q0 = LGAMMA_SMALL_Q[0] + z * LGAMMA_SMALL_Q[1];
            let q2 = LGAMMA_SMALL_Q[2] + z * LGAMMA_SMALL_Q[3];
            let q4 = LGAMMA_SMALL_Q[4] + z * LGAMMA_SMALL_Q[5];
            let q6 = LGAMMA_SMALL_Q[6] + z * LGAMMA_SMALL_Q[7];
            fl = z * ((q0 + z2 * q2) + z4 * (q4 + z2 * q6));
            fh = polydddfst(z, LGAMMA_SMALL_C0.len(), &LGAMMA_SMALL_C0, &mut fl);
            fh = mulddd2(x, fh, fl, &mut fl);
            fh = sumdd(-lh, -ll, fh, fl, &mut fl);
            eps = 1.5e-22;
        } else {
            // -0.5 < x < -0.03125 || 0.03125 < x < 0.5
            let mut xl = 0.0;
            let tf = fasttwosum(1.0, asdouble(t), &mut xl);
            au = ((tf).to_bits() >> 37) as u32;
            let ou = au - UBRD[0];
            let mut j =
                ((0x157ced865u64 - (ou as u64) * 0x150d) * (ou as u64) + 0x128000000000) >> 45;
            if au < UBRD[j as usize] {
                j -= 1;
            }
            let z = (tf - LGAMMA_OFFS[j as usize]) + xl;
            let z2 = z * z;
            let z4 = z2 * z2;
            let q = &LGAMMA_CL[j as usize];
            let q0 = q[0] + z * q[1];
            let q2 = q[2] + z * q[3];
            let q4 = q[4] + z * q[5];
            let q6 = q[6] + z * q[7];
            fl = z * ((q0 + z2 * q2) + z4 * (q4 + z2 * q6));
            fh = polydddfst(z, 5, &LGAMMA_CH[j as usize], &mut fl);
            if j == 4 {
                let z = -x;
                fh = mulddd2(z, fh, fl, &mut fl);
            }
            eps = fh.abs() * 8.3e-20;
            fh = sumdd(-lh, -ll, fh, fl, &mut fl);
            eps += lh.abs() * 5e-22;
        }
    } else {
        let ax = x.abs();
        if au >= UBRD[19] {
            // |x| >= 8.29541
            let mut ll = 0.0;
            let mut lh = as_logd(ax, &mut ll);
            lh -= 1.0;
            if au >= 0x2198000 {
                if au >= 0x3fabaa6 {
                    let mut tmp = 0.0;
                    lh = fasttwosum(lh, ll, &mut tmp);
                    ll = tmp;
                }
                let hlh = lh * 0.5;
                lh = mulddd2(ax, lh, ll, &mut ll);
                ll -= hlh;
            } else {
                lh = mulddd2(ax - 0.5, lh, ll, &mut ll);
            }

            lh = fastsum(lh, ll, LGAMMA_ASYM_C[0][0], LGAMMA_ASYM_C[0][1], &mut ll);
            let mut fh1;
            let mut fl1;
            if ax < f64::from_bits(0x47b0000000000000) {
                // 0x1p100
                let zh = 1.0 / ax;
                let zl = fma_internal(zh, -ax, 1.0) * zh;
                let z2h = zh * zh;
                let z4h = z2h * z2h;
                let q0 = LGAMMA_ASYM_Q[0] + z2h * LGAMMA_ASYM_Q[1];
                let q2 = LGAMMA_ASYM_Q[2] + z2h * LGAMMA_ASYM_Q[3];
                let q4 = LGAMMA_ASYM_Q[4];
                fl1 = z2h * (q0 + z4h * (q2 + z4h * q4));
                fh1 = fasttwosum(LGAMMA_ASYM_C[1][0], fl1, &mut fl1);
                fl1 += LGAMMA_ASYM_C[1][1];
                fh1 = muldd2(fh1, fl1, zh, zl, &mut fl1);
            } else {
                fh1 = 0.0;
                fl1 = 0.0;
            }
            fh = fastsum(lh, ll, fh1, fl1, &mut fl);
            eps = fh.abs() * 4.5e-20;
        } else {
            // x in [0.5, 8.29541]
            let ou = au - UBRD[0];
            let mut j =
                ((0x157ced865u64 - (ou as u64) * 0x150d) * (ou as u64) + 0x128000000000) >> 45;
            if au < UBRD[j as usize] {
                j -= 1;
            }
            let z = ax - LGAMMA_OFFS[j as usize];
            let z2 = z * z;
            let z4 = z2 * z2;
            let q = &LGAMMA_CL[j as usize];
            let q0 = q[0] + z * q[1];
            let q2 = q[2] + z * q[3];
            let q4 = q[4] + z * q[5];
            let q6 = q[6] + z * q[7];
            fl = z * ((q0 + z2 * q2) + z4 * (q4 + z2 * q6));
            fh = polydddfst(z, 5, &LGAMMA_CH[j as usize], &mut fl);
            if j == 4 {
                let z = 1.0 - ax;
                fh = mulddd2(z, fh, fl, &mut fl);
            }
            if j == 10 {
                let z = ax - 2.0;
                fh = mulddd2(z, fh, fl, &mut fl);
            }
            eps = fh.abs() * 8.3e-20 + 1e-24;
        }

        if (t >> 63) != 0 {
            let mut sl = 0.0;
            let sh = as_sinpipid(x - floor(x), &mut sl);
            let sh = mulddd2(-x, sh, sl, &mut sl);
            let mut ll = 0.0;
            let lh = as_logd(sh, &mut ll);
            ll += sl / sh;
            fh = -sumdd(fh, fl, lh, ll, &mut fl);
            fl = -fl;
            eps += lh.abs() * 4e-22;
            let k = fx as i64;
            *signgamp = 1 - 2 * ((k & 1) as i32);
        } else {
            *signgamp = 1;
        }
    }

    let ub = fh + (fl + eps);
    let lb = fh + (fl - eps);
    if ub != lb {
        return as_lgamma_accurate(x);
    }
    ub
}
#[inline(always)]
pub fn lgamma(x: f64) -> f64 {
    let mut sign = 1;
    ieee754_lgamma_r(x, &mut sign)
}

#[inline(always)]
fn tgamma_slow(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x == 0.0 {
        return if x.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            f64::NAN
        };
    }
    if x <= 0.0 && x == floor(x) {
        return f64::NAN;
    }
    let ax = x.abs();
    let frac = x - floor(x);
    if ax < 64.0 && frac == 0.5 {
        return tgamma_half(x);
    }
    if x > 0.0 && ax < 64.0 && (frac == 0.25 || frac == 0.75) {
        return tgamma_quarter(x, frac);
    }
    if x < 0.0 {
        let frac = x - floor(x);
        if frac == 0.25 || frac == 0.5 || frac == 0.75 {
            let mut sl = 0.0;
            let sh = sinpi_parts(x, &mut sl);
            let (yh, yl) = tgamma_pos_dd(1.0 - x);
            let mut pl = 0.0;
            let ph = muldd2(sh, sl, yh, yl, &mut pl);
            let mut r = 1.0 / ph;
            let mut e = fma_internal(ph, r, -1.0) + pl * r;
            r -= r * e;
            e = fma_internal(ph, r, -1.0) + pl * r;
            r -= r * e;
            return r;
        }
        let mut sl = 0.0;
        let sh = sinpi_parts(x, &mut sl);
        let sign = if sh.is_sign_negative() { -1.0 } else { 1.0 };
        let (lh, ll) = as_lgamma_accurate_dd(x);
        let (yh, yl) = exp_dd_parts(lh, ll);
        let y = yh + yl;
        return if sign < 0.0 { -y } else { y };
    }
    tgamma_pos(x)
}

#[inline(always)]
pub fn tgamma(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x == 0.5 {
        return f64::from_bits(0x3ffc5bf891b4ef6a);
    }
    if x == 0.0 {
        return if x.is_sign_negative() {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            f64::INFINITY
        } else {
            f64::NAN
        };
    }
    if x <= 0.0 && x == floor(x) {
        return f64::NAN;
    }
    if let Some(v) = tgamma_fast(x) {
        return v;
    }
    tgamma_slow(x)
}

#[inline(never)]
fn tgamma_pos(x: f64) -> f64 {
    let (yh, yl) = tgamma_pos_dd(x);
    yh + yl
}

#[inline(never)]
fn tgamma_pos_dd(x: f64) -> (f64, f64) {
    if x < 64.0 {
        let frac = x - floor(x);
        if frac == 0.5 {
            return (tgamma_half(x), 0.0);
        }
        if frac == 0.25 || frac == 0.75 {
            return tgamma_quarter_dd(x, frac);
        }
    }
    if x == floor(x) {
        let n = x as i32;
        if n <= 20 {
            let mut acc: u128 = 1;
            for i in 2..n {
                acc *= i as u128;
            }
            return (acc as f64, 0.0);
        }
        if n <= 170 {
            let mut hi = 1.0;
            let mut lo = 0.0;
            for i in 2..n {
                let mut l = 0.0;
                hi = mulddd2(i as f64, hi, lo, &mut l);
                lo = l;
            }
            let y = hi + lo;
            if y.is_finite() {
                return (hi, lo);
            }
        }
    }
    let (lh, ll) = as_lgamma_accurate_dd(x);
    exp_dd_parts(lh, ll)
}
