use rug::{Float, Integer};

use crate::{
    sampling::Sampling,
    util::{add, mulmod, ZZX},
};

#[derive(Debug, Clone)]
pub struct EncryptionScheme {
    /* Ring parameters */
    p: i32,
    q: i32,
    f: ZZX,

    /* Knuth-Yao discrete Gaussian sampler parameters */
    tailcut: f32,
    sigma: Float,
    center: Float,

    gauss: Sampling,
}

impl EncryptionScheme {
    fn poly_sampling(&self, a: &mut ZZX) {
        // TODO: check if match c code
        // int bound = ((int)tailcut)*to_int(sigma);
        // int center = to_int(center);
        let bound = (self.tailcut * self.sigma.clone().to_f32()).round() as i32;
        let center = self.center.to_f32().round() as i32;

        a.set_length(self.p as usize);
        for i in 0..self.p as usize {
            let mut sample = self.gauss.knuth_yao();
            while (sample >= (center + bound)) || (sample <= (center - bound)) {
                sample = self.gauss.knuth_yao();
            }
            a.set_coeff(i, Some(sample));
        }
    }

    fn _mod(&self, a: &mut ZZX) {
        for i in 0..self.p as usize {
            a.set_coeff(i, Some(_mod(a.coeff(i).clone(), self.q.into())));
        }
    }
}

impl EncryptionScheme {
    pub fn new(p: i32, q: i32, preicsion: u64, tailcut: f32, sigma: Float, center: Float) -> Self {
        // RR::SetPrecision(to_long(precision));

        let mut f = ZZX::new();
        f.set_length(p as usize + 1);
        f.set_coeff(p as usize, Some(1));
        f.set_coeff(1, Some(-1));
        f.set_coeff(0, Some(-1));

        // f[6] = Integer::from(1);
        // f[5] = Integer::from(1);
        // f[4] = Integer::from(-5);
        // f[3] = Integer::from(-4);
        // f[2] = Integer::from(6);
        // f[1] = Integer::from(3);
        // f[0] = Integer::from(-1);

        let gauss = Sampling::new(preicsion, tailcut, sigma.clone(), center.clone());

        Self {
            p,
            q,
            f,
            tailcut,
            sigma,
            center,
            gauss,
        }
    }

    pub fn new_with_instance(orig: &EncryptionScheme) -> Self {
        orig.clone()
    }

    fn key_generation(&self, a: &ZZX, r2: &mut ZZX, p1: &mut ZZX) {
        let mut c: ZZX = ZZX::new();
        let mut r1: ZZX = ZZX::new();

        c.set_length(self.p as usize);
        r1.set_length(self.p as usize);

        self.poly_sampling(&mut r1);
        self.poly_sampling(r2);

        c = mulmod(a, r2, &self.f);
        *p1 = r1 - c;

        self._mod(p1);
    }

    fn encode(&self, aprime: &mut ZZX, a: Vec<i32>) {
        aprime.set_length(self.p as usize);

        let bound = (self.q - 1) / 2;
        for i in 0..self.p as usize {
            aprime[i] = Integer::from(a[i] * bound);
        }
    }

    fn decode(&self, a: &mut Vec<i32>, aprime: &ZZX) {
        let lbound = Integer::from((self.q - 1) / 4);
        let ubound = Integer::from(3 * lbound.clone());

        for i in 0..self.p as usize {
            if aprime[i] >= lbound && aprime[i] < ubound {
                a[i] = 1;
            } else {
                a[i] = 0;
            }
        }
    }

    fn encryption(&self, c1: &mut ZZX, c2: &mut ZZX, a: &ZZX, p1: &ZZX, m: &ZZX) {
        c1.set_length(self.p as usize);
        c2.set_length(self.p as usize);

        let mut add = ZZX::new();
        add.set_length(self.p as usize);

        let mut mult = ZZX::new();
        mult.set_length(self.p as usize);

        let mut e1 = ZZX::new();
        let mut e2 = ZZX::new();
        let mut e3 = ZZX::new();

        self.poly_sampling(&mut e1);
        self.poly_sampling(&mut e2);
        self.poly_sampling(&mut e3);

        add = e3 + m;
        mult = mulmod(p1, &e1, &self.f);
        *c2 = &mult + &add;
        mult = mulmod(a, &e1, &self.f);
        *c1 = mult + e2;

        self._mod(c1);
        self._mod(c2);
    }

    fn decrtyption(&self, m: &mut ZZX, c1: &ZZX, c2: &ZZX, r2: &ZZX) {
        let mut mult = ZZX::new();
        mult.set_length(self.p as usize);

        m.set_length(self.p as usize);

        mult = mulmod(c1, r2, &self.f);
        *m = mult - c2;

        self._mod(m);
    }
}

fn _mod(i: Integer, n: Integer) -> Integer {
    (i % n.clone() + n.clone()) % n
}
