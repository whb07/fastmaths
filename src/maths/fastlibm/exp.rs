use super::{f64_from_bits, f64_to_bits, is_inf_bits, is_nan_bits};

// ========= glibc-derived exp table (N=128) =========

const EXP_TAB_U64: [u64; 256] = [
    0x0000000000000000u64,
    0x3ff0000000000000u64,
    0x3c9b3b4f1a88bf6eu64,
    0x3feff63da9fb3335u64,
    0xbc7160139cd8dc5du64,
    0x3fefec9a3e778061u64,
    0xbc905e7a108766d1u64,
    0x3fefe315e86e7f85u64,
    0x3c8cd2523567f613u64,
    0x3fefd9b0d3158574u64,
    0xbc8bce8023f98efau64,
    0x3fefd06b29ddf6deu64,
    0x3c60f74e61e6c861u64,
    0x3fefc74518759bc8u64,
    0x3c90a3e45b33d399u64,
    0x3fefbe3ecac6f383u64,
    0x3c979aa65d837b6du64,
    0x3fefb5586cf9890fu64,
    0x3c8eb51a92fdeffcu64,
    0x3fefac922b7247f7u64,
    0x3c3ebe3d702f9cd1u64,
    0x3fefa3ec32d3d1a2u64,
    0xbc6a033489906e0bu64,
    0x3fef9b66affed31bu64,
    0xbc9556522a2fbd0eu64,
    0x3fef9301d0125b51u64,
    0xbc5080ef8c4eea55u64,
    0x3fef8abdc06c31ccu64,
    0xbc91c923b9d5f416u64,
    0x3fef829aaea92de0u64,
    0x3c80d3e3e95c55afu64,
    0x3fef7a98c8a58e51u64,
    0xbc801b15eaa59348u64,
    0x3fef72b83c7d517bu64,
    0xbc8f1ff055de323du64,
    0x3fef6af9388c8deau64,
    0x3c8b898c3f1353bfu64,
    0x3fef635beb6fcb75u64,
    0xbc96d99c7611eb26u64,
    0x3fef5be084045cd4u64,
    0x3c9aecf73e3a2f60u64,
    0x3fef54873168b9aau64,
    0xbc8fe782cb86389du64,
    0x3fef4d5022fcd91du64,
    0x3c8a6f4144a6c38du64,
    0x3fef463b88628cd6u64,
    0x3c807a05b0e4047du64,
    0x3fef3f49917ddc96u64,
    0x3c968efde3a8a894u64,
    0x3fef387a6e756238u64,
    0x3c875e18f274487du64,
    0x3fef31ce4fb2a63fu64,
    0x3c80472b981fe7f2u64,
    0x3fef2b4565e27cddu64,
    0xbc96b87b3f71085eu64,
    0x3fef24dfe1f56381u64,
    0x3c82f7e16d09ab31u64,
    0x3fef1e9df51fdee1u64,
    0xbc3d219b1a6fbffau64,
    0x3fef187fd0dad990u64,
    0x3c8b3782720c0ab4u64,
    0x3fef1285a6e4030bu64,
    0x3c6e149289cecb8fu64,
    0x3fef0cafa93e2f56u64,
    0x3c834d754db0abb6u64,
    0x3fef06fe0a31b715u64,
    0x3c864201e2ac744cu64,
    0x3fef0170fc4cd831u64,
    0x3c8fdd395dd3f84au64,
    0x3feefc08b26416ffu64,
    0xbc86a3803b8e5b04u64,
    0x3feef6c55f929ff1u64,
    0xbc924aedcc4b5068u64,
    0x3feef1a7373aa9cbu64,
    0xbc9907f81b512d8eu64,
    0x3feeecae6d05d866u64,
    0xbc71d1e83e9436d2u64,
    0x3feee7db34e59ff7u64,
    0xbc991919b3ce1b15u64,
    0x3feee32dc313a8e5u64,
    0x3c859f48a72a4c6du64,
    0x3feedea64c123422u64,
    0xbc9312607a28698au64,
    0x3feeda4504ac801cu64,
    0xbc58a78f4817895bu64,
    0x3feed60a21f72e2au64,
    0xbc7c2c9b67499a1bu64,
    0x3feed1f5d950a897u64,
    0x3c4363ed60c2ac11u64,
    0x3feece086061892du64,
    0x3c9666093b0664efu64,
    0x3feeca41ed1d0057u64,
    0x3c6ecce1daa10379u64,
    0x3feec6a2b5c13cd0u64,
    0x3c93ff8e3f0f1230u64,
    0x3feec32af0d7d3deu64,
    0x3c7690cebb7aafb0u64,
    0x3feebfdad5362a27u64,
    0x3c931dbdeb54e077u64,
    0x3feebcb299fddd0du64,
    0xbc8f94340071a38eu64,
    0x3feeb9b2769d2ca7u64,
    0xbc87deccdc93a349u64,
    0x3feeb6daa2cf6642u64,
    0xbc78dec6bd0f385fu64,
    0x3feeb42b569d4f82u64,
    0xbc861246ec7b5cf6u64,
    0x3feeb1a4ca5d920fu64,
    0x3c93350518fdd78eu64,
    0x3feeaf4736b527dau64,
    0x3c7b98b72f8a9b05u64,
    0x3feead12d497c7fdu64,
    0x3c9063e1e21c5409u64,
    0x3feeab07dd485429u64,
    0x3c34c7855019c6eau64,
    0x3feea9268a5946b7u64,
    0x3c9432e62b64c035u64,
    0x3feea76f15ad2148u64,
    0xbc8ce44a6199769fu64,
    0x3feea5e1b976dc09u64,
    0xbc8c33c53bef4da8u64,
    0x3feea47eb03a5585u64,
    0xbc845378892be9aeu64,
    0x3feea34634ccc320u64,
    0xbc93cedd78565858u64,
    0x3feea23882552225u64,
    0x3c5710aa807e1964u64,
    0x3feea155d44ca973u64,
    0xbc93b3efbf5e2228u64,
    0x3feea09e667f3bcdu64,
    0xbc6a12ad8734b982u64,
    0x3feea012750bdabfu64,
    0xbc6367efb86da9eeu64,
    0x3fee9fb23c651a2fu64,
    0xbc80dc3d54e08851u64,
    0x3fee9f7df9519484u64,
    0xbc781f647e5a3ecfu64,
    0x3fee9f75e8ec5f74u64,
    0xbc86ee4ac08b7db0u64,
    0x3fee9f9a48a58174u64,
    0xbc8619321e55e68au64,
    0x3fee9feb564267c9u64,
    0x3c909ccb5e09d4d3u64,
    0x3feea0694fde5d3fu64,
    0xbc7b32dcb94da51du64,
    0x3feea11473eb0187u64,
    0x3c94ecfd5467c06bu64,
    0x3feea1ed0130c132u64,
    0x3c65ebe1abd66c55u64,
    0x3feea2f336cf4e62u64,
    0xbc88a1c52fb3cf42u64,
    0x3feea427543e1a12u64,
    0xbc9369b6f13b3734u64,
    0x3feea589994cce13u64,
    0xbc805e843a19ff1eu64,
    0x3feea71a4623c7adu64,
    0xbc94d450d872576eu64,
    0x3feea8d99b4492edu64,
    0x3c90ad675b0e8a00u64,
    0x3feeaac7d98a6699u64,
    0x3c8db72fc1f0eab4u64,
    0x3feeace5422aa0dbu64,
    0xbc65b6609cc5e7ffu64,
    0x3feeaf3216b5448cu64,
    0x3c7bf68359f35f44u64,
    0x3feeb1ae99157736u64,
    0xbc93091fa71e3d83u64,
    0x3feeb45b0b91ffc6u64,
    0xbc5da9b88b6c1e29u64,
    0x3feeb737b0cdc5e5u64,
    0xbc6c23f97c90b959u64,
    0x3feeba44cbc8520fu64,
    0xbc92434322f4f9aau64,
    0x3feebd829fde4e50u64,
    0xbc85ca6cd7668e4bu64,
    0x3feec0f170ca07bau64,
    0x3c71affc2b91ce27u64,
    0x3feec49182a3f090u64,
    0x3c6dd235e10a73bbu64,
    0x3feec86319e32323u64,
    0xbc87c50422622263u64,
    0x3feecc667b5de565u64,
    0x3c8b1c86e3e231d5u64,
    0x3feed09bec4a2d33u64,
    0xbc91bbd1d3bcbb15u64,
    0x3feed503b23e255du64,
    0x3c90cc319cee31d2u64,
    0x3feed99e1330b358u64,
    0x3c8469846e735ab3u64,
    0x3feede6b5579fdbfu64,
    0xbc82dfcd978e9db4u64,
    0x3feee36bbfd3f37au64,
    0x3c8c1a7792cb3387u64,
    0x3feee89f995ad3adu64,
    0xbc907b8f4ad1d9fau64,
    0x3feeee07298db666u64,
    0xbc55c3d956dcaebau64,
    0x3feef3a2b84f15fbu64,
    0xbc90a40e3da6f640u64,
    0x3feef9728de5593au64,
    0xbc68d6f438ad9334u64,
    0x3feeff76f2fb5e47u64,
    0xbc91eee26b588a35u64,
    0x3fef05b030a1064au64,
    0x3c74ffd70a5fddcdu64,
    0x3fef0c1e904bc1d2u64,
    0xbc91bdfbfa9298acu64,
    0x3fef12c25bd71e09u64,
    0x3c736eae30af0cb3u64,
    0x3fef199bdd85529cu64,
    0x3c8ee3325c9ffd94u64,
    0x3fef20ab5fffd07au64,
    0x3c84e08fd10959acu64,
    0x3fef27f12e57d14bu64,
    0x3c63cdaf384e1a67u64,
    0x3fef2f6d9406e7b5u64,
    0x3c676b2c6c921968u64,
    0x3fef3720dcef9069u64,
    0xbc808a1883ccb5d2u64,
    0x3fef3f0b555dc3fau64,
    0xbc8fad5d3ffffa6fu64,
    0x3fef472d4a07897cu64,
    0xbc900dae3875a949u64,
    0x3fef4f87080d89f2u64,
    0x3c74a385a63d07a7u64,
    0x3fef5818dcfba487u64,
    0xbc82919e2040220fu64,
    0x3fef60e316c98398u64,
    0x3c8e5a50d5c192acu64,
    0x3fef69e603db3285u64,
    0x3c843a59ac016b4bu64,
    0x3fef7321f301b460u64,
    0xbc82d52107b43e1fu64,
    0x3fef7c97337b9b5fu64,
    0xbc892ab93b470dc9u64,
    0x3fef864614f5a129u64,
    0x3c74b604603a88d3u64,
    0x3fef902ee78b3ff6u64,
    0x3c83c5ec519d7271u64,
    0x3fef9a51fbc74c83u64,
    0xbc8ff7128fd391f0u64,
    0x3fefa4afa2a490dau64,
    0xbc8dae98e223747du64,
    0x3fefaf482d8e67f1u64,
    0x3c8ec3bc41aa2008u64,
    0x3fefba1bee615a27u64,
    0x3c842b94c3a9eb32u64,
    0x3fefc52b376bba97u64,
    0x3c8a64a931d185eeu64,
    0x3fefd0765b6e4540u64,
    0xbc8e37bae43be3edu64,
    0x3fefdbfdad9cbe14u64,
    0x3c77893b4d91cd9du64,
    0x3fefe7c1819e90d8u64,
    0x3c5305c14160cc89u64,
    0x3feff3c22b8f71f1u64,
];

// ========= exp(x) =========

const EXP_TABLE_BITS: u32 = 7;
const N: u64 = 1u64 << EXP_TABLE_BITS;
const INV_LN2_N: f64 = f64::from_bits(0x3ff71547652b82fe) * 128.0;
const NEG_LN2_HI_N: f64 = f64::from_bits(0xbf762e42fefa0000);
const NEG_LN2_LO_N: f64 = f64::from_bits(0xbd0cf79abc9e3b3a);
const SHIFT: f64 = f64::from_bits(0x4338000000000000);

const EXP_C2: f64 = f64::from_bits(0x3fdffffffffffdbd);
const EXP_C3: f64 = f64::from_bits(0x3fc555555555543c);
const EXP_C4: f64 = f64::from_bits(0x3fa55555cf172b91);
const EXP_C5: f64 = f64::from_bits(0x3f81111167a4d017);

const EXP_HI: f64 = 709.782712893384;
const EXP_LO: f64 = -745.1332191019411;

#[inline(always)]
fn specialcase(tmp: f64, sbits: u64, k: i64) -> f64 {
    if k > 0 {
        let sbits = sbits.wrapping_sub(1009u64 << 52);
        let scale = f64_from_bits(sbits);
        return f64::from_bits(0x7f00_0000_0000_0000) * (scale + scale * tmp);
    }

    let sbits = sbits.wrapping_add(1022u64 << 52);
    let scale = f64_from_bits(sbits);
    let mut y = scale + scale * tmp;
    if y < 1.0 {
        let lo = scale - y + scale * tmp;
        let hi = 1.0 + y;
        let lo = 1.0 - hi + y + lo;
        y = (hi + lo) - 1.0;
        if y == 0.0 {
            y = 0.0;
        }
    }
    f64::from_bits(0x0010_0000_0000_0000) * y
}

#[inline(always)]
pub fn exp(x: f64) -> f64 {
    let ux = f64_to_bits(x);
    if is_nan_bits(ux) {
        return x;
    }
    if is_inf_bits(ux) {
        return if x.is_sign_negative() {
            0.0
        } else {
            f64::INFINITY
        };
    }
    if x > EXP_HI {
        return f64::INFINITY;
    }
    if x < EXP_LO {
        return 0.0;
    }

    let z = x * INV_LN2_N;
    let kd = z + SHIFT;
    let ki = f64_to_bits(kd);
    let kd = kd - SHIFT;
    let k = kd as i64;
    let r = x + kd * NEG_LN2_HI_N + kd * NEG_LN2_LO_N;

    let idx = (ki % N) as usize * 2;
    let top = ki << (52 - EXP_TABLE_BITS);
    let tail = f64_from_bits(EXP_TAB_U64[idx]);
    let sbits = EXP_TAB_U64[idx + 1].wrapping_add(top);
    let scale = f64_from_bits(sbits);

    let r2 = r * r;
    let tmp = tail + r + r2 * (EXP_C2 + r * EXP_C3) + r2 * r2 * (EXP_C4 + r * EXP_C5);
    if k <= -(1023 * N as i64) || k >= (1024 * N as i64) {
        return specialcase(tmp, sbits, k);
    }
    scale + scale * tmp
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_exp_near_limits() {
        let values = [EXP_HI, EXP_HI - 1e-10, EXP_LO, EXP_LO + 1e-10];
        for &x in &values {
            let actual = exp(x);
            let expected = x.exp();
            if expected.is_infinite() {
                assert!(actual.is_infinite());
            } else {
                let diff = (actual - expected).abs();
                let rel_diff = diff / expected;
                assert!(
                    rel_diff < 1e-15,
                    "exp({x}) failed: got {actual}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn test_exp_subnormals() {
        // Values that result in subnormal outputs
        // e^x < 2^-1022 => x < -708.396...
        for i in 0..100 {
            let x = -708.4 - (i as f64) * 0.1;
            if x < EXP_LO {
                break;
            }
            let actual = exp(x);
            let expected = x.exp();
            let diff = (actual - expected).abs();
            // For subnormals, absolute difference is more appropriate as relative diff can be large
            assert!(
                diff < 1e-300,
                "exp({x}) subnormal failed: got {actual}, expected {expected}"
            );
        }
    }
}
