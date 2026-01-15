use super::{f64_from_bits, f64_to_bits};

// ========= glibc-derived log tables (N=128) =========

const LOG_TABLE_BITS: u32 = 7;
const N: u64 = 1u64 << LOG_TABLE_BITS;
const OFF: u64 = 0x3fe6000000000000;

const LOG_INVC_U64: [u64; 128] = [
    0x3ff734f0c3e0de9fu64,
    0x3ff713786a2ce91fu64,
    0x3ff6f26008fab5a0u64,
    0x3ff6d1a61f138c7du64,
    0x3ff6b1490bc5b4d1u64,
    0x3ff69147332f0cbau64,
    0x3ff6719f18224223u64,
    0x3ff6524f99a51ed9u64,
    0x3ff63356aa8f24c4u64,
    0x3ff614b36b9ddc14u64,
    0x3ff5f66452c65c4cu64,
    0x3ff5d867b5912c4fu64,
    0x3ff5babccb5b90deu64,
    0x3ff59d61f2d91a78u64,
    0x3ff5805612465687u64,
    0x3ff56397cee76bd3u64,
    0x3ff54725e2a77f93u64,
    0x3ff52aff42064583u64,
    0x3ff50f22dbb2bddfu64,
    0x3ff4f38f4734ded7u64,
    0x3ff4d843cfde2840u64,
    0x3ff4bd3ec078a3c8u64,
    0x3ff4a27fc3e0258au64,
    0x3ff4880524d48434u64,
    0x3ff46dce1b192d0bu64,
    0x3ff453d9d3391854u64,
    0x3ff43a2744b4845au64,
    0x3ff420b54115f8fbu64,
    0x3ff40782da3ef4b1u64,
    0x3ff3ee8f5d57fe8fu64,
    0x3ff3d5d9a00b4ce9u64,
    0x3ff3bd60c010c12bu64,
    0x3ff3a5242b75dab8u64,
    0x3ff38d22cd9fd002u64,
    0x3ff3755bc5847a1cu64,
    0x3ff35dce49ad36e2u64,
    0x3ff34679984dd440u64,
    0x3ff32f5cceffcb24u64,
    0x3ff3187775a10d49u64,
    0x3ff301c8373e3990u64,
    0x3ff2eb4ebb95f841u64,
    0x3ff2d50a0219a9d1u64,
    0x3ff2bef9a8b7fd2au64,
    0x3ff2a91c7a0c1babu64,
    0x3ff293726014b530u64,
    0x3ff27dfa5757a1f5u64,
    0x3ff268b39b1d3bbfu64,
    0x3ff2539d838ff5bdu64,
    0x3ff23eb7aac9083bu64,
    0x3ff22a012ba940b6u64,
    0x3ff2157996cc4132u64,
    0x3ff201201dd2fc9bu64,
    0x3ff1ecf4494d480bu64,
    0x3ff1d8f5528f6569u64,
    0x3ff1c52311577e7cu64,
    0x3ff1b17c74cb26e9u64,
    0x3ff19e010c2c1ab6u64,
    0x3ff18ab07bb670bdu64,
    0x3ff1778a25efbcb6u64,
    0x3ff1648d354c31dau64,
    0x3ff151b990275fddu64,
    0x3ff13f0ea432d24cu64,
    0x3ff12c8b7210f9dau64,
    0x3ff11a3028ecb531u64,
    0x3ff107fbda8434afu64,
    0x3ff0f5ee0f4e6bb3u64,
    0x3ff0e4065d2a9fceu64,
    0x3ff0d244632ca521u64,
    0x3ff0c0a77ce2981au64,
    0x3ff0af2f83c636d1u64,
    0x3ff09ddb98a01339u64,
    0x3ff08cabaf52e7dfu64,
    0x3ff07b9f2f4e28fbu64,
    0x3ff06ab58c358f19u64,
    0x3ff059eea5ecf92cu64,
    0x3ff04949cdd12c90u64,
    0x3ff038c6c6f0ada9u64,
    0x3ff02865137932a9u64,
    0x3ff0182427ea7348u64,
    0x3ff008040614b195u64,
    0x3fefe01ff726fa1au64,
    0x3fefa11cc261ea74u64,
    0x3fef6310b081992eu64,
    0x3fef25f63ceeadcdu64,
    0x3feee9c8039113e7u64,
    0x3feeae8078cbb1abu64,
    0x3fee741aa29d0c9bu64,
    0x3fee3a91830a99b5u64,
    0x3fee01e009609a56u64,
    0x3fedca01e577bb98u64,
    0x3fed92f20b7c9103u64,
    0x3fed5cac66fb5cceu64,
    0x3fed272caa5ede9du64,
    0x3fecf26e3e6b2ccdu64,
    0x3fecbe6da2a77902u64,
    0x3fec8b266d37086du64,
    0x3fec5894bd5d5804u64,
    0x3fec26b533bb9f8cu64,
    0x3febf583eeece73fu64,
    0x3febc4fd75db96c1u64,
    0x3feb951e0c864a28u64,
    0x3feb65e2c5ef3e2cu64,
    0x3feb374867c9888bu64,
    0x3feb094b211d304au64,
    0x3feadbe885f2ef7eu64,
    0x3feaaf1d31603da2u64,
    0x3fea82e63fd358a7u64,
    0x3fea5740ef09738bu64,
    0x3fea2c2a90ab4b27u64,
    0x3fea01a01393f2d1u64,
    0x3fe9d79f24db3c1bu64,
    0x3fe9ae2505c7b190u64,
    0x3fe9852ef297ce2fu64,
    0x3fe95cbaeea44b75u64,
    0x3fe934c69de74838u64,
    0x3fe90d4f2f6752e6u64,
    0x3fe8e6528effd79du64,
    0x3fe8bfce9fcc007cu64,
    0x3fe899c0dabec30eu64,
    0x3fe87427aa2317fbu64,
    0x3fe84f00acb39a08u64,
    0x3fe82a49e8653e55u64,
    0x3fe8060195f40260u64,
    0x3fe7e22563e0a329u64,
    0x3fe7beb377dcb5adu64,
    0x3fe79baa679725c2u64,
    0x3fe77907f2170657u64,
    0x3fe756cadbd6130cu64,
];

const LOG_LOGC_U64: [u64; 128] = [
    0xbfd7cc7f79e69000u64,
    0xbfd76feec20d0000u64,
    0xbfd713e31351e000u64,
    0xbfd6b85b38287800u64,
    0xbfd65d5590807800u64,
    0xbfd602d076180000u64,
    0xbfd5a8ca86909000u64,
    0xbfd54f4356035000u64,
    0xbfd4f637c36b4000u64,
    0xbfd49da7fda85000u64,
    0xbfd445923989a800u64,
    0xbfd3edf439b0b800u64,
    0xbfd396ce448f7000u64,
    0xbfd3401e17bda000u64,
    0xbfd2e9e2ef468000u64,
    0xbfd2941b3830e000u64,
    0xbfd23ec58cda8800u64,
    0xbfd1e9e129279000u64,
    0xbfd1956d2b48f800u64,
    0xbfd141679ab9f800u64,
    0xbfd0edd094ef9800u64,
    0xbfd09aa518db1000u64,
    0xbfd047e65263b800u64,
    0xbfcfeb224586f000u64,
    0xbfcf474a7517b000u64,
    0xbfcea4443d103000u64,
    0xbfce020d44e9b000u64,
    0xbfcd60a22977f000u64,
    0xbfccc00104959000u64,
    0xbfcc202956891000u64,
    0xbfcb81178d811000u64,
    0xbfcae2c9ccd3d000u64,
    0xbfca45402e129000u64,
    0xbfc9a877681df000u64,
    0xbfc90c6d69483000u64,
    0xbfc87120a645c000u64,
    0xbfc7d68fb4143000u64,
    0xbfc73cb83c627000u64,
    0xbfc6a39a9b376000u64,
    0xbfc60b3154b7a000u64,
    0xbfc5737d76243000u64,
    0xbfc4dc7b8fc23000u64,
    0xbfc4462c51d20000u64,
    0xbfc3b08abc830000u64,
    0xbfc31b996b490000u64,
    0xbfc2875490a44000u64,
    0xbfc1f3b9f879a000u64,
    0xbfc160c8252ca000u64,
    0xbfc0ce7f57f72000u64,
    0xbfc03cdc49fea000u64,
    0xbfbf57bdbc4b8000u64,
    0xbfbe370896404000u64,
    0xbfbd17983ef94000u64,
    0xbfbbf9674ed8a000u64,
    0xbfbadc79202f6000u64,
    0xbfb9c0c3e7288000u64,
    0xbfb8a646b372c000u64,
    0xbfb78d01b3ac0000u64,
    0xbfb674f145380000u64,
    0xbfb55e0e6d878000u64,
    0xbfb4485cdea1e000u64,
    0xbfb333d94d6aa000u64,
    0xbfb22079f8c56000u64,
    0xbfb10e4698622000u64,
    0xbfaffa6c6ad20000u64,
    0xbfadda8d4a774000u64,
    0xbfabbcece4850000u64,
    0xbfa9a1894012c000u64,
    0xbfa788583302c000u64,
    0xbfa5715e67d68000u64,
    0xbfa35c8a49658000u64,
    0xbfa149e364154000u64,
    0xbf9e72c082eb8000u64,
    0xbf9a55f152528000u64,
    0xbf963d62cf818000u64,
    0xbf9228fb8caa0000u64,
    0xbf8c317b20f90000u64,
    0xbf8419355daa0000u64,
    0xbf781203c2ec0000u64,
    0xbf60040979240000u64,
    0x3f6feff384900000u64,
    0x3f87dc41353d0000u64,
    0x3f93cea3c4c28000u64,
    0x3f9b9fc114890000u64,
    0x3fa1b0d8ce110000u64,
    0x3fa58a5bd001c000u64,
    0x3fa95c8340d88000u64,
    0x3fad276aef578000u64,
    0x3fb07598e598c000u64,
    0x3fb253f5e30d2000u64,
    0x3fb42edd8b380000u64,
    0x3fb606598757c000u64,
    0x3fb7da76356a0000u64,
    0x3fb9ab434e1c6000u64,
    0x3fbb78c7bb0d6000u64,
    0x3fbd431332e72000u64,
    0x3fbf0a3171de6000u64,
    0x3fc067152b914000u64,
    0x3fc147858292b000u64,
    0x3fc2266ecdca3000u64,
    0x3fc303d7a6c55000u64,
    0x3fc3dfc33c331000u64,
    0x3fc4ba366b7a8000u64,
    0x3fc5933928d1f000u64,
    0x3fc66acd2418f000u64,
    0x3fc740f8ec669000u64,
    0x3fc815c0f51af000u64,
    0x3fc8e92954f68000u64,
    0x3fc9bb3602f84000u64,
    0x3fca8bed1c2c0000u64,
    0x3fcb5b515c01d000u64,
    0x3fcc2967ccbcc000u64,
    0x3fccf635d5486000u64,
    0x3fcdc1bd3446c000u64,
    0x3fce8c01b8cfe000u64,
    0x3fcf5509c0179000u64,
    0x3fd00e6c121fb800u64,
    0x3fd071b80e93d000u64,
    0x3fd0d46b9e867000u64,
    0x3fd13687334bd000u64,
    0x3fd1980d67234800u64,
    0x3fd1f8ffe0cc8000u64,
    0x3fd2595fd7636800u64,
    0x3fd2b9300914a800u64,
    0x3fd3187210436000u64,
    0x3fd377266dec1800u64,
    0x3fd3d54ffbaf3000u64,
    0x3fd432eee32fe000u64,
];

// ========= ln(x) =========

const LN2_HI: f64 = f64::from_bits(0x3fe62e42fefa3800);
const LN2_LO: f64 = f64::from_bits(0x3d2ef35793c76730);

const LOG_A0: f64 = f64::from_bits(0xbfe0000000000001);
const LOG_A1: f64 = f64::from_bits(0x3fd555555551305b);
const LOG_A2: f64 = f64::from_bits(0xbfcfffffffeb4590);
const LOG_A3: f64 = f64::from_bits(0x3fc999b324f10111);
const LOG_A4: f64 = f64::from_bits(0xbfc55575e506c89f);

const LOG_B0: f64 = f64::from_bits(0xbfe0000000000000);
const LOG_B1: f64 = f64::from_bits(0x3fd5555555555577);
const LOG_B2: f64 = f64::from_bits(0xbfcffffffffffdcb);
const LOG_B3: f64 = f64::from_bits(0x3fc999999995dd0c);
const LOG_B4: f64 = f64::from_bits(0xbfc55555556745a7);
const LOG_B5: f64 = f64::from_bits(0x3fc24924a344de30);
const LOG_B6: f64 = f64::from_bits(0xbfbfffffa4423d65);
const LOG_B7: f64 = f64::from_bits(0x3fbc7184282ad6ca);
const LOG_B8: f64 = f64::from_bits(0xbfb999eb43b068ff);
const LOG_B9: f64 = f64::from_bits(0x3fb78182f7afd085);
const LOG_B10: f64 = f64::from_bits(0xbfb5521375d145cd);

#[inline(always)]
pub fn ln(x: f64) -> f64 {
    let mut ix = f64_to_bits(x);
    let top = (ix >> 48) as u32;

    const LO: u64 = 0x3fee000000000000; // 1 - 2^-4
    const HI: u64 = 0x3ff1090000000000; // 1 + 0x1.09p-4

    if ix.wrapping_sub(LO) < (HI - LO) {
        if ix == f64_to_bits(1.0) {
            return 0.0;
        }
        let r = x - 1.0;
        let r2 = r * r;
        let r3 = r * r2;
        let y = r3
            * (LOG_B1
                + r * LOG_B2
                + r2 * LOG_B3
                + r3 * (LOG_B4
                    + r * LOG_B5
                    + r2 * LOG_B6
                    + r3 * (LOG_B7 + r * LOG_B8 + r2 * LOG_B9 + r3 * LOG_B10)));
        let w = r * f64::from_bits(0x4190000000000000); // 2^27
        let rhi = r + w - w;
        let rlo = r - rhi;
        let w = rhi * rhi * LOG_B0;
        let hi = r + w;
        let mut lo = r - hi + w;
        lo += LOG_B0 * rlo * (rhi + r);
        return y + lo + hi;
    }

    if top.wrapping_sub(0x0010) >= 0x7ff0 - 0x0010 {
        if (ix << 1) == 0 {
            return f64::NEG_INFINITY;
        }
        if ix == 0x7ff0_0000_0000_0000 {
            return f64::INFINITY;
        }
        if (top & 0x8000) != 0 || (top & 0x7ff0) == 0x7ff0 {
            return f64::NAN;
        }
        ix = f64_to_bits(x * f64::from_bits(0x4330_0000_0000_0000));
        ix = ix.wrapping_sub(52u64 << 52);
    }

    let tmp = ix.wrapping_sub(OFF);
    let i = ((tmp >> (52 - LOG_TABLE_BITS)) & (N - 1)) as usize;
    let k = ((tmp as i64) >> 52) as i32;
    let iz = ix.wrapping_sub(tmp & (0xfff_u64 << 52));
    let invc = f64_from_bits(LOG_INVC_U64[i]);
    let logc = f64_from_bits(LOG_LOGC_U64[i]);
    let z = f64_from_bits(iz);

    let r = z.mul_add(invc, -1.0);
    let kd = k as f64;

    let w = kd * LN2_HI + logc;
    let hi = w + r;
    let lo = w - hi + r + kd * LN2_LO;

    let r2 = r * r;
    lo + r2 * LOG_A0 + r * r2 * (LOG_A1 + r * LOG_A2 + r2 * (LOG_A3 + r * LOG_A4)) + hi
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ln_near_one() {
        let values = [
            0.9375, // LO boundary
            0.94, 0.99, 1.0, 1.01, 1.06, 1.0625, // Near HI boundary
        ];
        for &x in &values {
            let actual = ln(x);
            let expected = x.ln();
            let diff = (actual - expected).abs();
            assert!(
                diff < 1e-15,
                "ln({x}) failed: got {actual}, expected {expected}"
            );
        }
    }

    #[test]
    fn test_ln_table_boundaries() {
        // Test values that cross table indices
        for i in 0..127 {
            let tmp = (i as u64) << (52 - LOG_TABLE_BITS);
            let ix = tmp + OFF;
            let x = f64_from_bits(ix);
            let actual = ln(x);
            let expected = x.ln();
            let diff = (actual - expected).abs();
            assert!(diff < 1e-15, "ln({x}) at index {i} failed");

            // Just above boundary
            let x_plus = f64_from_bits(ix + 1);
            let actual_plus = ln(x_plus);
            let expected_plus = x_plus.ln();
            assert!((actual_plus - expected_plus).abs() < 1e-15);
        }
    }
}
