use rug::Integer;

/// Custom clone of NTL::ZZX 
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ZZX {
    coeffs: Vec<Integer>,
}

impl ZZX {
    /************************
     * Constructors, Assignment
     ************************/
    /// init 0
    pub fn new() -> Self {
        ZZX { coeffs: vec![] }
    }

    /// init with single coefficient
    pub fn new_with_val<T: Into<Integer>>(n: T) -> Self {
        ZZX { coeffs: vec![n.into()] }
    }

    /// init with vector of coefficients
    pub fn new_with_vec<T: Into<Integer>>(coeffs: Vec<T>) -> Self {
        let coeffs = coeffs.into_iter().map(Into::into).collect();
        ZZX { coeffs }
    }

    /// intial value 0, but space is pre-allocated for n coefficients
    pub fn new_with_size(n: usize) -> Self {
        let coeffs = Vec::with_capacity(n);
        ZZX { coeffs }
    }

    /// Strip leading zeros
    pub fn normalize(&mut self) {
        let mut i = self.coeffs.len() - 1;
        while i > 0 && self.coeffs[i] == 0 {
            self.coeffs.pop();
            i -= 1;
        }
    }

    pub fn set_length(&mut self, len: usize) {
        self.coeffs.resize(len, Integer::from(0));
    }

    /*************************
     * Some utility functions
     *************************/
    pub fn deg(&self) -> i64 {
        self.coeffs.len() as i64 - 1
    }

    pub fn coeff(&self, i: usize) -> Integer {
        if i >= self.coeffs.len() {
            Integer::from(0)
        } else {
            self.coeffs[i].clone()
        }
    }

    pub fn get_coeff(&self, i: usize) -> Integer {
        if i >= self.coeffs.len() {
            Integer::from(0)
        } else {
            self.coeffs[i].clone()
        }
    }

    pub fn set_coeff<T>(&mut self, i: usize, n: Option<T>) where T: Into<Integer> {
        self.coeffs[i] = n.map(Into::into).unwrap_or(Integer::from(1));
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn lead_coeff(&self) -> Integer {
        if self.is_zero() {
            Integer::from(0)
        } else {
            self.coeffs[self.deg() as usize].clone()   
        }
    }

    pub fn max_size(&self) -> u32 {
        let mut res = 0;
        for i in 0..self.coeffs.len() {
            let t = self.coeffs[i].significant_bits();
            if t > res {
                res = t;
            }
        }
        res
    }

    pub fn max_bits(&self) -> u32 {
        let mut m = 0;
        for i in 0..self.deg() as usize {
            m = m.max(self.coeffs[i].significant_bits());
        }
        m
    }
}

pub fn mulmod(x: &mut ZZX, a: &ZZX, b: &ZZX, f: &ZZX) {    
    if a.deg() >= f.deg() || b.deg() >= f.deg() || f.deg() == 0 || f.lead_coeff() != *Integer::ONE {
        panic!("MulMod: bad args");
    }

    let mut t = ZZX::new();
    mul(&mut t, a, b);
    rem(x, &t, f);
}

pub fn mul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    if a.is_zero() || b.is_zero() {
        c.set_length(0);
        return;
    }
    if a == b {
        sqr(c, a);
        return;
    }

    let maxa: u32 = a.max_size(); // MaxSize(a)
    let maxb: u32 = b.max_size(); // MaxSize(b)

    let k = maxa.min(maxb);
    let s = a.deg().min(b.deg()) + 1;
   
    // // FIXME: I should have a way of setting all these crossovers
    // // automatically

    // if s == 1 || (k == 1 && s < 40) || (k == 2 && s < 20) || (k == 3 && s < 10) {
    //     PlainMul(c, a, b);
    // } else if s < 80 || (k < 30 && s < 150) {
    //     KarMul(c, a, b);
    // } else if choose_ss(a.deg(), a.max_bits(), b.deg(), b.max_bits()){
    //     SSMul(c, a, b);
    // } else {
    //     HomMul(c, a, b);
    // }

    // TODO: uncomment the above.
    PlainMul(c, a, b);
}

pub fn rem(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    todo!("impl `void Rem(ZZX& r, const ZZX& a, const ZZX& b)` func");
}

fn sqr(c: &mut ZZX, a: &ZZX) {
    if a.is_zero() {
        c.set_length(0);
        return;
    }

    let maxa = a.max_size();
    let k = maxa;
    let s = a.deg() + 1;
    
    // if s == 1 || (k == 1 && s < 50) || (k == 2 && s < 25) || (k == 3 && s < 25) || (k == 4 && s < 10) {
    //     PlainSqr(c, a);
    // } else if s < 80 || (k < 30 && s < 150) {
    //     KarSqr(c, a);
    // } else if choose_ss(a.deg(), a.max_bits(), a.deg(), a.max_bits()) {
    //     SSSqr(c, a);
    // } else {
    //     HomSqr(c, a);
    // }

    // TODO: uncomment the above.
    PlainSqr(c, a);
}

fn PlainMul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    if a == b {
        PlainSqr(c, a);
        return;
    }

    let da = a.deg();
    let db = b.deg();

    if da < 0 || db < 0 {
        c.set_length(0);
        return;
    }

    let d = da + db;
    let ap: &Vec<Integer> = &a.coeffs;
    let bp: &Vec<Integer> = &b.coeffs;

    c.set_length(d as usize + 1);

    for i in 0..=d {
        let j_min = 0.max(i - db);
        let j_max = da.min(i);
        let mut accum = Integer::from(0);
        for j in j_min..=j_max {
            accum += ap[j as usize].clone() * bp[i as usize - j as usize].clone();
        }
        c.set_coeff(i as usize, Some(accum));
    }

    c.normalize();
}

fn KarMul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    // if a.is_zero() || b.is_zero() {
    //     c.set_length(0);
    //     return;
    // }

    // if a == b {
    //     KarSqr(c, a);
    //     return;
    // }

    // let sa = a.coeffs.len();
    // let sb = b.coeffs.len();

    // let ap = &a.coeffs;
    // let bp = &b.coeffs;

    // c.set_length(sa + sb - 1);

    todo!("impl `void KarMul(ZZX& c, const ZZX& a, const ZZX& b)` func");
}

fn SSMul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    todo!("impl `void SSMul(ZZX& c, const ZZX& a, const ZZX& b)` func");
}

fn HomMul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    todo!("impl `void HomMul(ZZX& c, const ZZX& a, const ZZX& b)` func");
}

fn choose_ss(a_deg: i64, a_max_bits: u32, b_deg: i64, b_max_bits: u32) -> bool {
    todo!("impl `bool choose_ss(int a_deg, int a_max_bits, int b_deg, int b_max_bits)` func");
}

fn PlainSqr(c: &mut ZZX, a: &ZZX) {
    let da = a.deg();
    if da < 0 {
        c.set_length(0);
        return;
    }

    let d = 2 * a.deg();
    let ap = &a.coeffs;

    c.set_length(d as usize + 1);

    for i in 0..=d {
        let j_min = 0.max(i - da);
        let j_max = i.min(da);
        let m = j_max - j_min + 1;
        let m2 = m >> 1;
        let j_max = j_min + m2 - 1;
        
        let mut accum = Integer::from(0);
        let mut t = Integer::from(0);
        for j in j_min..=j_max {
            t = ap[j as usize].clone() * ap[i as usize - j as usize].clone();
            accum += t;
        }
        accum *= 2;
        if m & 1 == 1 {
            t = ap[j_max as usize + 1].clone().square();
            accum += t;
        }
        c.set_coeff(i as usize, Some(accum));
    }
    c.normalize();
}

fn KarSqr(c: &mut ZZX, a: &ZZX) {
    todo!("impl `void KarSqr(ZZX& c, const ZZX& a)` func");
}

fn SSSqr(c: &mut ZZX, a: &ZZX) {
    todo!("impl `void SSSqr(ZZX& c, const ZZX& a)` func");
}

fn HomSqr(c: &mut ZZX, a: &ZZX) {
    todo!("impl `void HomSqr(ZZX& c, const ZZX& a)` func");
}
