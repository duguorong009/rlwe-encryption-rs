use rug::{Float, Integer};

use crate::{sampling::Sampling, util::ZZX};

#[derive(Debug, Clone)]
pub struct EncryptionScheme {
    /* Ring parameters */
    p: i32,
    q: i32,
    f: ZZX,

    /* Knuth-Yao discrete Gaussian sampler parameters */
    tail_cut: f32,
    sigma: Float,
    center: Float,

    gauss: Sampling,
}

impl EncryptionScheme {
    fn set_F() {
        todo!("impl `void SetF()` func");
    }

    fn poly_sampling(a: &ZZX) {
        todo!("impl `void PolySampling(ZZX& a)` func");
    }

    fn _Mod(a: &ZZX) {
        todo!("impl `void Mod(ZZX& a)` func");
    }

    fn _mod(i: Integer, n: Integer) {
        todo!("impl `inline ZZ mod(ZZ& i, IZZ& n)` func");
    }
}

impl EncryptionScheme {
    pub fn new(p: i32, q: i32, preicsion: u64, tail_cut: f32, sigma: Float, center: Float) -> Self {
        // RR::SetPrecision(to_long(precision));

        let mut f = ZZX::new();
        f.set_length(p as usize + 1);
        f.set_coeff(p as usize, Some(1));
        f.set_coeff(1, Some(-1));
        f.set_coeff(0, Some(-1));

        let gauss = Sampling::new(preicsion, tail_cut, sigma.clone(), center.clone());

        Self {
            p,
            q,
            f,
            tail_cut,
            sigma,
            center,
            gauss,
        }
    }

    pub fn new_with_instance(orig: &EncryptionScheme) -> Self {
        orig.clone()
    }

    fn key_generation(a: &ZZX, r2: &ZZX, p1: &ZZX) {
        todo!("impl `void KeyGeneration(ZZX& a, ZZX& r2, ZZX& p1)` func");
    }

    fn encode(aprime: &ZZX, a: Vec<i32>) {
        todo!("impl `void Encode(ZZX& aprime, const int32_t a[])` func");
    }

    fn decode(a: Vec<i32>, aprime: &ZZX) {
        todo!("impl `void Decode(int32_t a[], const ZZX& aprime)` func");
    }

    fn encryption(c1: &ZZX, c2: &ZZX, a: &ZZX, p1: &ZZX, m: &ZZX) {
        todo!("impl `void Encryption(ZZX& c1, ZZX& c2, const ZZX& a, const ZZX& p1, const ZZX& m)` func");
    }

    fn decrtyption(m: &ZZX, c1: &ZZX, c2: &ZZX, r2: &ZZX) {
        todo!("impl `void Decryption(ZZX& m, const ZZX& c1, const ZZX& c2, const ZZX& r2)` func");
    }
}