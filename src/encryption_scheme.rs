use rug::{Float, Integer};

use crate::{sampling::Sampling, util::ZZX};

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
        todo!("impl `EncryptionScheme(const int& p, const int& q, long precision, float tailcut, RR sigma, RR center);` func");
    }

    pub fn new_with_instance(orig: &EncryptionScheme) -> Self {
        todo!("impl `EncryptionScheme(const EncryptionScheme& orig);` func");
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