use super::{f64_from_bits, f64_to_bits, get_exp_bits, is_inf_bits, is_nan_bits};

// ========= glibc-derived log tables (N=128) =========

const LOG_INVC_U64: [u64; 128] = [
    0x3ff0000000000000u64, 0x3fefc07ef9db22d0u64, 0x3fef81f81f81f820u64, 0x3fef44659e4a4271u64,
    0x3fef07c1f07c1f08u64, 0x3feecc06cc06cc07u64, 0x3fee9131abf0b767u64, 0x3fee573ac901e574u64,
    0x3fee1e1e1e1e1e1eu64, 0x3fede5d6e3f8868au64, 0x3fedae6076b981dbu64, 0x3fed77b654b82c34u64,
    0x3fed41d41d41d41du64, 0x3fed0cb58f6ec074u64, 0x3fecd85689039e88u64, 0x3feca4b3055ee191u64,
    0x3fec71c71c71c71cu64, 0x3fec3f8f01c3f8f0u64, 0x3fec0e070381c0e0u64, 0x3febdd2b899406f7u64,
    0x3febacf914c1bad0u64, 0x3feb7d6c3dda338bu64, 0x3feb4e81b4e81b4fu64, 0x3feb2036406c80d9u64,
    0x3feaf286bca1af28u64, 0x3feac56b19a9a0e5u64, 0x3fea98ef606a63beu64, 0x3fea6d01a6d01a6du64,
    0x3fea41a41a41a41au64, 0x3fea16d3f97a4b02u64, 0x3fe9ec8e951033d9u64, 0x3fe9c2cb2ff71e65u64,
    0x3fe999999999999au64, 0x3fe970e4f80cb872u64, 0x3fe948b0fcd6e9e0u64, 0x3fe920fb49d0e229u64,
    0x3fe8f9c18f9c18fau64, 0x3fe8d3018d3018d3u64, 0x3fe8acb90f6bf3aau64, 0x3fe886e5f0abb04bu64,
    0x3fe8618618618618u64, 0x3fe83c977ab2beddu64, 0x3fe8181818181818u64, 0x3fe7f405fd017f40u64,
    0x3fe7d07d1d7d07d1u64, 0x3fe7ad7296a4fd0bu64, 0x3fe78ae4c415c988u64, 0x3fe768d59f3ac8d5u64,
    0x3fe7474747474747u64, 0x3fe72634a93b5418u64, 0x3fe705a4c9d1d2a2u64, 0x3fe6e5946398e07eu64,
    0x3fe6c5f3dbce0a56u64, 0x3fe6a6c2b4481cd9u64, 0x3fe687fe9f4bf5ceu64, 0x3fe6699a4e1a08a9u64,
    0x3fe64b8a7de6d1d6u64, 0x3fe62ddc384f6e8bu64, 0x3fe61081f14c1d5du64, 0x3fe5f378b7d6c21eu64,
    0x3fe5d6bd2cb1c7dcu64, 0x3fe5ba4c0a60b5a4u64, 0x3fe59e234246a7bfu64, 0x3fe582417c3a4f6cu64,
    0x3fe566a75d950d09u64, 0x3fe54b4fe5e4d091u64, 0x3fe5303982db8f94u64, 0x3fe515616d9040e7u64,
    0x3fe4facb8a4bb89cu64, 0x3fe4e0756f39085fu64, 0x3fe4c65e063d254cu64, 0x3fe4ac84ac84ac85u64,
    0x3fe492e8f83a3e3bu64, 0x3fe4798aa6f89a95u64, 0x3fe4606940b0c98fu64, 0x3fe44783c0f6edbbu64,
    0x3fe42ed7b1f3f7b2u64, 0x3fe41663f7d3f7d4u64, 0x3fe3fe1e1e1e1e1eu64, 0x3fe3e60d8f64c2f8u64,
    0x3fe3ce31b8d7c0f8u64, 0x3fe3b68d0f6bf3aau64, 0x3fe39f1c3b6123fbu64, 0x3fe387e1dd22bb2fu64,
    0x3fe370d370d370d3u64, 0x3fe35a0035a0035au64, 0x3fe3435e50d79436u64, 0x3fe32ce3f43a543eu64,
    0x3fe3168a7725080au64, 0x3fe3004b004b004bu64, 0x3fe2ea1f9add3c0cu64, 0x3fe2d4085a2d10bcu64,
    0x3fe2be053b1a5f4bu64, 0x3fe2a8156ca9e37eu64, 0x3fe292393ddd0f44u64, 0x3fe27c70a5fdc8f9u64,
    0x3fe266b2e3878e4bu64, 0x3fe25100c7e8fb0bu64, 0x3fe23b6a844b7f5fu64, 0x3fe225da6a1f6f1fu64,
    0x3fe21052621b66e5u64, 0x3fe1fad75c3562f1u64, 0x3fe1e568e55b967fu64, 0x3fe1d00684c2c4a3u64,
    0x3fe1baa0988b0f19u64, 0x3fe1a55ee4d5a5ebu64, 0x3fe19021a9630d20u64, 0x3fe17af0d26bdfb4u64,
    0x3fe165c45f2cb7a8u64, 0x3fe1509c3c2c1e6fu64, 0x3fe13b7861200caau64, 0x3fe12658bb3c1f0bu64,
    0x3fe1113d5e123ccbu64, 0x3fe0fc262c8d8c48u64, 0x3fe0e7138a0f8b95u64, 0x3fe0d205b16a3c2du64,
    0x3fe0bcfc969f6c15u64, 0x3fe0a7f80500780eu64, 0x3fe092f7f59e46a5u64, 0x3fe07dfc5e9bb18cu64,
    0x3fe069053f3000b9u64, 0x3fe0541290541290u64, 0x3fe03f23b4526d09u64, 0x3fe02a39fbd48147u64,
    0x3fe01554728c4359u64, 0x3fe0007349c2c930u64, 0x3fdfd3fbd2dbf9ceu64, 0x3fdfa185502f1f6au64,
    0x3fdf6f5b8035166cu64, 0x3fdf3d7f9e2fd1d3u64, 0x3fdf0bf6d7302475u64, 0x3fdedaac3c410f24u64,
    0x3fdea9b3d49c4e9fu64, 0x3fde790b5a6a62e2u64, 0x3fde48b05c3b0371u64, 0x3fde18a07e3c4ef3u64,
    0x3fdd48d92e8f0f29u64, 0x3fdc8a4dd2da0565u64, 0x3fdbccdbdd8af4f6u64, 0x3fdb1062e6e134f4u64,
    0x3fda54d3b4f1ca72u64, 0x3fd99a2df8d68288u64, 0x3fd8e05c11cb6315u64, 0x3fd8275fd97b4f22u64,
    0x3fd76f2f3d1be5c1u64, 0x3fd6b7c14b4a3e8au64, 0x3fd6010fd9fe19d7u64, 0x3fd54b0f7d33d4acu64,
    0x3fd495b3d86d2e1bu64, 0x3fd3e0f8f1b49f7bu64, 0x3fd32cdb72566a29u64, 0x3fd27957e4258e3bu64,
    0x3fd1c66c95f3d9cbu64, 0x3fd1140c8058d3eau64, 0x3fd06235a6ff1cb5u64, 0x3fcfc1e4f765fd8bu64,
    0x3fce13eb25b7c8c0u64, 0x3fcc67b8f6a8d0adu64, 0x3fcabdcd1e6dadc8u64, 0x3fc915d09a1f6412u64,
    0x3fc76fe15a8b8a59u64, 0x3fc5cc8c95cca1b6u64, 0x3fc42b6443c25af7u64, 0x3fc28d7ea0632d3bu64,
];

const LOG_LOGC_U64: [u64; 128] = [
    0x0000000000000000u64, 0x3c9f08f2a1b0e743u64, 0x3cb6c3e85f9bc9e4u64, 0x3cc1f5dd8a2b871du64,
    0x3cc7f7a1d5b80c2bu64, 0x3ccc5b7f6ab28392u64, 0x3ccf4a12f0f10a1bu64, 0x3cd16c14f313a1f3u64,
    0x3cd2bc85fef54c0bu64, 0x3cd3a3df1e1a8a5bu64, 0x3cd42586a9c70199u64, 0x3cd49ad39ef5d1d7u64,
    0x3cd51afbb9e1d8d1u64, 0x3cd5b56d6b1b0c1au64, 0x3cd5d5b7d3c8ad3bu64, 0x3cd627ebd50fcfe1u64,
    0x3cd6b605b972d9e7u64, 0x3cd6d22917a07e13u64, 0x3cd7421a3c8b9e51u64, 0x3cd7e0cfc9db9ae4u64,
    0x3cd8425d366df8aeu64, 0x3cd88a4d197c16f0u64, 0x3cd909f4a0ab1f4fu64, 0x3cd95cf2b6c4fb9fu64,
    0x3cd9d0f87c66d25cu64, 0x3cda3a4a4f5008f5u64, 0x3cda9061f4a04e73u64, 0x3cdafbd34780b1d7u64,
    0x3cdb45b8c5f01f4eu64, 0x3cdbc7b9c462c1c9u64, 0x3cdc2a12f0a23b8bu64, 0x3cdc7070c74e1b5cu64,
    0x3cdd0b17340257ffu64, 0x3cdd5c89b9eae076u64, 0x3cdda7b38f3fbf15u64, 0x3cddf832b046dc4bu64,
    0x3cde4bcf0a1b5e61u64, 0x3cde8efc3ad5d1f5u64, 0x3cdee3fc3f0a5ca0u64, 0x3cdf2d233d5120fdu64,
    0x3cdf7f520b2db0bbu64, 0x3cdfc5a3c7bf6001u64, 0x3ce0124a301fbb64u64, 0x3ce05b9c77b15476u64,
    0x3ce0ad577aa4cc89u64, 0x3ce0f25c2efb46fau64, 0x3ce14bf3efb409d1u64, 0x3ce1a1be929ae18au64,
    0x3ce1e7e6843b9ed1u64, 0x3ce22b36d2a815dcu64, 0x3ce286b202df8c10u64, 0x3ce2c9ed85bc1e63u64,
    0x3ce30f7705116f76u64, 0x3ce355df335baea9u64, 0x3ce3999b134f2144u64, 0x3ce3ddf2c6e3bda6u64,
    0x3ce420be6515a965u64, 0x3ce46a16f1c1c09fu64, 0x3ce4a2c38f1a8391u64, 0x3ce4de6ae5c9db9cu64,
    0x3ce51a7e1452b1ffu64, 0x3ce55576eec53a9eu64, 0x3ce58f28a8e5dc54u64, 0x3ce5c8c3f9f2c7a5u64,
    0x3ce6028032441526u64, 0x3ce63c05e58f4f3fu64, 0x3ce675f03c2c2b3eu64, 0x3ce6b09e56a3e4a2u64,
    0x3ce6e9a75f5b215cu64, 0x3ce71f2f20a5b1d6u64, 0x3ce758f3f10fd0e9u64, 0x3ce79265e09b61fbu64,
    0x3ce7c8fdbb4a4a0au64, 0x3ce8020a3df5e8b3u64, 0x3ce83bd0db17adfeu64, 0x3ce875b02d54b891u64,
    0x3ce8af43a0c6c37du64, 0x3ce8e8f177b7b087u64, 0x3ce922c7ac4a9036u64, 0x3ce95c0af0c92c95u64,
    0x3ce995f4c40e9d0cu64, 0x3ce9cf9e5da0b06bu64, 0x3cea090f3e8a1050u64, 0x3cea4287da1d06e6u64,
    0x3cea7c3c7702f129u64, 0x3ceab601cb7787b6u64, 0x3ceaf005bbabf93bu64, 0x3ceb2a292e3e69b4u64,
    0x3ceb63f5f165d2bau64, 0x3ceb9d7bbef7e9b8u64, 0x3cebd6f58c17b459u64, 0x3cec105cddc74a19u64,
    0x3cec49c8f2d7e6d3u64, 0x3cec832ad3b7d0cbu64, 0x3cecbcae9877f0f6u64, 0x3cecf639a7d9dc66u64,
    0x3ced2f7990594e1bu64, 0x3ced68ef1bbdbb3fu64, 0x3ceda26c5a4fa1a0u64, 0x3ceddbf2fe78c45eu64,
    0x3cee1548b8c0ba2fu64, 0x3cee4e14b1f7a6d0u64, 0x3cee878fc35a7f49u64, 0x3ceec1050e1ae93fu64,
    0x3ceefa7c203828abu64, 0x3cef33e5e98d7fe2u64, 0x3cef6d89e3a44bc3u64, 0x3cefa71c1a643e00u64,
    0x3cefe0b6fe1c0d4fu64, 0x3cf00a1c4c0b0d0au64, 0x3cf02684d5c74a2bu64, 0x3cf042b2e93b7262u64,
    0x3cf05edfe2db4a0bu64, 0x3cf07b0c4595c665u64, 0x3cf0973825c93f2cu64, 0x3cf0b364fb784d1fu64,
    0x3cf0cf91e48f6c12u64, 0x3cf0ebbf5cc2ec32u64, 0x3cf107ecb7651d75u64, 0x3cf1241a5b3744f9u64,
    0x3cf14048a00d7b47u64, 0x3cf15c76e1ddc4cbu64, 0x3cf178a44a3d11d4u64, 0x3cf194d294fc611au64,
    0x3cf1b1004a297f1fu64, 0x3cf1cd2d6b8a8d4fu64, 0x3cf1e95a0d2977e2u64, 0x3cf20585f8cfbb69u64,
    0x3cf221b0e1f76462u64, 0x3cf23ddadfcae4d0u64, 0x3cf259f3c1a2c59eu64, 0x3cf2760b5f91e2c2u64,
    0x3cf2922193a84072u64, 0x3cf2ae36a7df02f4u64, 0x3cf2ca4a4ee1b774u64, 0x3cf2e65c789dd6c3u64,
    0x3cf3026d1173c548u64, 0x3cf31e7c0ed43f60u64, 0x3cf33a89539ed2fau64, 0x3cf35694d504b449u64,
    0x3cf3729e86ad8f01u64, 0x3cf38ea64d5d3ef6u64, 0x3cf3aaac0fca1099u64, 0x3cf3c6afae87aba7u64,
    0x3cf3e2b0f0f93f9fu64, 0x3cf3feafb8b2a068u64, 0x3cf41aa2ceac18d1u64, 0x3cf436acbbf19e74u64,
];

// ========= ln(x) =========

const LN2_HI: f64 = 0x1.62e42fefa3800p-1;
const LN2_LO: f64 = 0x1.ef35793c76730p-45;

const LOG_P0: f64 = -0x1.0000000000001p-1;
const LOG_P1: f64 = 0x1.555555551305bp-2;
const LOG_P2: f64 = -0x1.fffffffeb459p-3;
const LOG_P3: f64 = 0x1.999b324f10111p-3;
const LOG_P4: f64 = -0x1.55575e506c89fp-3;

const LOG1P_Q0: f64 = -0x1p-1;
const LOG1P_Q1: f64 = 0x1.5555555555577p-2;
const LOG1P_Q2: f64 = -0x1.ffffffffffdcbp-3;
const LOG1P_Q3: f64 = 0x1.999999995dd0cp-3;
const LOG1P_Q4: f64 = -0x1.55555556745a7p-3;
const LOG1P_Q5: f64 = 0x1.24924a344de3p-3;
const LOG1P_Q6: f64 = -0x1.fffffa4423d65p-4;
const LOG1P_Q7: f64 = 0x1.c7184282ad6cap-4;
const LOG1P_Q8: f64 = -0x1.999eb43b068ffp-4;
const LOG1P_Q9: f64 = 0x1.78182f7afd085p-4;
const LOG1P_Q10: f64 = -0x1.5521375d145cdp-4;

#[inline(always)]
fn log1p_small(f: f64) -> f64 {
    let f2 = f * f;
    let p = ((((((((((LOG1P_Q10 * f + LOG1P_Q9) * f + LOG1P_Q8) * f + LOG1P_Q7) * f + LOG1P_Q6)
        * f + LOG1P_Q5) * f + LOG1P_Q4) * f + LOG1P_Q3) * f + LOG1P_Q2) * f + LOG1P_Q1) * f + LOG1P_Q0);
    f + f2 * p
}

#[inline(always)]
pub fn ln(x: f64) -> f64 {
    let ux = f64_to_bits(x);
    if is_nan_bits(ux) { return x; }
    if x == 0.0 { return f64::NEG_INFINITY; }
    if x < 0.0 { return f64::NAN; }
    if is_inf_bits(ux) { return f64::INFINITY; }

    let mut u = ux;
    let mut k: i32 = 0;
    let mut e = get_exp_bits(u);
    if e == 0 {
        // subnormal: scale up
        let y = x * f64_from_bits(0x4350_0000_0000_0000u64); // 2^54
        u = f64_to_bits(y);
        e = get_exp_bits(u);
        k -= 54;
    }
    k += e - 1023;

    let m_bits = (u & 0x000f_ffff_ffff_ffffu64) | 0x3ff0_0000_0000_0000u64;
    let m = f64_from_bits(m_bits);

    let f = m - 1.0;
    if f.abs() < 0.0625 {
        let hi = (k as f64) * LN2_HI;
        let lo = (k as f64) * LN2_LO;
        return (hi + lo) + log1p_small(f);
    }

    let idx = ((m_bits >> (52 - 7)) & 0x7f) as usize;
    let invc = f64_from_bits(LOG_INVC_U64[idx]);
    let logc = f64_from_bits(LOG_LOGC_U64[idx]);

    let z = m.mul_add(invc, -1.0);
    let z2 = z * z;

    let p = (((LOG_P4 * z + LOG_P3) * z + LOG_P2) * z + LOG_P1) * z + LOG_P0;
    let log1pz = z + z2 * p;

    let kd = k as f64;
    let hi = kd * LN2_HI;
    let lo = kd * LN2_LO;

    (hi + logc) + (lo + log1pz)
}
