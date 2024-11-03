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

