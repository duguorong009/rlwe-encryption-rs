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
        ZZX {
            coeffs: vec![n.into()],
        }
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
    /// degree of polynomial, note that zero polynomial has degree -1
    pub fn deg(&self) -> i64 {
        self.coeffs.len() as i64 - 1
    }
    /// zero if i not in range
    pub fn coeff(&self, i: usize) -> Integer {
        if i >= self.coeffs.len() {
            Integer::from(0)
        } else {
            self.coeffs[i].clone()
        }
    }

    /// return coeff[i], or 0 if i not in range
    pub fn get_coeff(&self, i: usize) -> Integer {
        if i >= self.coeffs.len() {
            Integer::from(0)
        } else {
            self.coeffs[i].clone()
        }
    }

    /// zero if poly is zero
    pub fn lead_coeff(&self) -> Integer {
        if self.is_zero() {
            Integer::from(0)
        } else {
            self.coeffs[self.deg() as usize].clone()
        }
    }

    /// zero if poly is zero
    pub fn const_term(&self) -> Integer {
        if self.is_zero() {
            Integer::from(0)
        } else {
            self.coeffs[0].clone()
        }
    }

    // set coeff[i] = n or 1
    pub fn set_coeff<T>(&mut self, i: usize, n: Option<T>)
    where
        T: Into<Integer>,
    {
        let m = self.deg();
        if i > m as usize && self.is_zero() {
            return;
        }

        if i > m as usize {
            self.set_length(i + 1);
        }
        self.coeffs[i] = n.map(Into::into).unwrap_or(Integer::from(1));

        self.normalize();
    }

    /// set to the monomial X
    pub fn set_x() -> Self {
        let mut c = ZZX::new();
        c.set_length(2);
        c.set_coeff(1, Some(1));
        c
    }

    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn is_x(&self) -> bool {
        self.deg() == 1 && self.coeffs[0] == 0 && self.coeffs[1] == 1
    }

    pub fn clear(&mut self) {
        self.coeffs.clear();
    }

    pub fn set(&mut self) {
        self.coeffs = vec![Integer::from(1)];
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

pub fn rem(r: &mut ZZX, a: &ZZX, b: &ZZX) {
    let da = a.deg();
    let db = b.deg();

    if db < 0 {
        panic!("rem: division by zero");
    }

    if da < db {
        r.coeffs = a.coeffs.clone();
    } else if db == 0 {
        const_rem(r, a, b.const_term());
    } else if b.lead_coeff() == 1 {
        pseudo_rem(r, a, b);
    } else if b.lead_coeff() == -1 {
        let b1 = b.neg();
        pseudo_rem(r, a, &b1);
    } else if divide(a, b) {
        r.coeffs = vec![Integer::from(0)];
    } else {
        let mut r1 = ZZX::new();
        pseudo_rem(&mut r1, a, b);
        power(m, b.lead_coeff(), da - db + 1);
        if !divide(r, r1, m) {
            panic!("rem: remainder not defined over ZZ");
        }
    }
}

fn const_rem(r: &mut ZZX, a: &ZZX, b: Integer) {
    if b == 0 {
        panic!("const_rem: division by zero");
    }
    r.coeffs = vec![Integer::from(0)];
}

fn pseudo_rem(r: &mut ZZX, a: &ZZX, b: &ZZX) {
    plain_pseudo_rem(r, a, b);
}

fn plain_pseudo_rem(r: &mut ZZX, a: &ZZX, b: &ZZX) {
    let mut q = ZZX::new();
    plain_pseudo_div_rem(&mut q, r, a, b);
}

fn plain_pseudo_div_rem(q: &mut ZZX, r: &mut ZZX, a: &ZZX, b: &ZZX) {
    let da = a.deg();
    let db = b.deg();

    if db < 0 {
        panic!("pseudo_div_rem: division by zero");
    }

    if da < db {
        r.coeffs = a.coeffs.clone();
        q.clear();
        return;
    }

    let bp = &b.coeffs;
    let lc = bp[db as usize].clone();
    let lc_is_one = lc == 1;

    let mut xp = a.coeffs.clone();
    let dq = da - db;

    q.set_length(dq as usize + 1);
    let mut qp = q.coeffs.clone();

    if !lc_is_one {
        let mut t = lc.clone();
        for i in (0..=dq - 1).rev() {
            xp[i as usize] = xp[i as usize].clone() * t.clone();
            if i > 0 {
                t = t * lc.clone();
            }
        }
    }

    for i in (0..=dq).rev() {
        let t = xp[(i + db) as usize].clone();
        qp[i as usize] = t.clone();
        for j in (0..db).rev() {
            let s = t.clone() * bp[j as usize].clone();
            if !lc_is_one {
                xp[i as usize + j as usize] = xp[i as usize + j as usize].clone() * lc.clone();
            }
            xp[i as usize + j as usize] = xp[i as usize + j as usize].clone() - s.clone();
        }
    }

    if !lc_is_one {
        let mut t = lc.clone();
        for i in 1..=dq {
            qp[i as usize] = qp[i as usize].clone() * t.clone();
            if i < dq {
                t = t.clone() * lc.clone();
            }
        }
    }

    r.set_length(db as usize);
    for i in 0..db {
        r.coeffs[i as usize] = xp[i as usize].clone();
    }
    r.normalize();
}

fn divide(q: &mut ZZX, a: &ZZX, b: &ZZX) -> bool {
    let da = a.deg();
    let db = b.deg();

    if db <= 8 || da - db <= 8 {
        plain_divide(q, a, b)
    } else {
        hom_divide(q, a, b)
    }
}

fn plain_divide(qq: &mut ZZX, aa: &ZZX, bb: &ZZX) -> bool {
    if bb.is_zero() {
        if aa.is_zero() {
            qq.clear();
            return true; // 1
        } else {
            return false; // 0
        }
    }

    if bb.deg() == 0 {
        divide(qq, aa, bb.const_term())
    }


    let da = aa.deg();
    let db = bb.deg();

    if da < db {
        return false;  // 0
    }

    let mut ca = Integer::from(0);
    let mut cb = Integer::from(0);
    content(&mut ca, aa);
    content(&mut cb, bb);

    if !divide(cq, ca, cb) {
        return false; // 0
    }

    let mut a = ZZX::new();
    let mut b = ZZX::new();
    divide(&mut a, aa, bb);
    divide(&mut b, bb, cb);

    if !divide(a.lead_coeff(), b.lead_coeff()) {
        return false; // 0
    }

    if !divide(a.const_term(), b.const_term()) {
        return false; // 0
    }

    let coeff_bnd = a.max_bits() + ((da + 1).num_bits() + 1)/2  + (da - db);

    let mut bp = b.coeffs.clone();
    let lc = bp[db as usize].clone();

    let lc_is_one = lc == 1;

    let mut xp = a.coeffs.clone();

    let dq = da - db;
    let mut q = ZZX::new();
    q.set_length(dq as usize + 1);
    
    for i in (0..=dq).rev() {
        if !lc_is_one {
            if !divide(t, xp[i + db], lc) {
                return 0;
            }
        } else {
            t = xp[i + db].clone();
        }

        if t.num_bits() > coeff_bnd {
            return 0;
        }

        q.coeffs[i as usize] = t.clone();
        for j in (0..db).rev() {
            s = t.clone() * bp[j as usize].clone();
            xp[i as usize + j as usize] = xp[i as usize + j as usize].clone() - s.clone();
        }
    }

    for i in 0..db {
        if !xp[i].is_zero() {
            return 0;
        }
    }

    mul(qq, q, cq);
    return 1;
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

/// x = a % X^m
fn _trunc(x: &mut ZZX, a: &ZZX, m: usize) {
    let n = m.min(a.coeffs.len());
    x.set_length(n);

    for i in 0..n {
        x.coeffs[i] = a.coeffs[i].clone();
    }
    x.normalize();
}

pub fn trunc(a: &ZZX, m: usize) -> ZZX {
    let mut x = ZZX::new();
    _trunc(&mut x, a, m);
    x
}

/// x = a / X^n
fn _right_shift(x: &mut ZZX, a: &ZZX, n: i64) {
    if a.is_zero() {
        x.clear();
        return;
    }

    if n < 0 {
        // handle case n < -MAX_LONG
        _left_shift(x, a, -n);
        return;
    }

    let da = a.deg();
    if da < n {
        x.clear();
        return;
    }

    let n = n as usize;
    let da = da as usize;

    x.set_length(a.coeffs.len() - n);

    for i in 0..da - n {
        x.coeffs[i] = a.coeffs[i + n].clone();
    }

    x.normalize();
}

pub fn right_shift(a: &ZZX, n: i64) -> ZZX {
    let mut x = ZZX::new();
    _right_shift(&mut x, a, n);
    x
}

/// x = a * X^n
fn _left_shift(x: &mut ZZX, a: &ZZX, n: i64) {
    if a.is_zero() {
        x.clear();
        return;
    }

    if n < 0 {
        // handle case n < -MAX_LONG
        _right_shift(x, a, -n);
        return;
    }

    let n = n as usize;
    let m = a.coeffs.len();
    x.set_length(m + n);
    for i in (0..m).rev() {
        x.coeffs[i + n] = a.coeffs[i].clone();
    }
    for i in 0..n {
        x.coeffs[i] = Integer::from(0);
    }
    x.normalize();
}

pub fn left_shift(a: &ZZX, n: i64) -> ZZX {
    let mut x = ZZX::new();
    _left_shift(&mut x, a, n);
    x
}

/* TODO: implement operators >>, <<, >>=, <<= with above shift funcs */

/// x = derivative of a
fn _diff(x: &mut ZZX, a: &ZZX) {
    let n = a.deg();
    if n <= 0 {
        x.clear();
        return;
    }

    x.set_length(n as usize);

    for i in 0..n as usize {
        x.coeffs[i] = a.coeffs[i + 1].clone() * (i + 1) as i64;
    }
    x.set_length(n as usize);

    x.normalize();
}

pub fn diff(a: &ZZX) -> ZZX {
    let mut x = ZZX::new();
    _diff(&mut x, a);
    x
}

// /// x = a^{-1} % X^m
// fn _inv_trunc(x: &mut ZZX, a: &ZZX, m: i64) {
//     if m < 0 {
//         panic!("inv_trunc: m < 0");
//     }

//     if m == 0 {
//         x.clear();
//         return;
//     }

//     // TODO: add ntl overflow check `NTL_OVERFLOW(e, 1, 0)`

//     newton_inv_trunc(x, a, m);
// }

// pub fn inv_trunc(a: &ZZX, m: i64) -> ZZX {
//     let mut x = ZZX::new();
//     _inv_trunc(&mut x, a, m);
//     x
// }

// fn newton_inv_trunc(c: &mut ZZX, a: &ZZX, m: i64) {
//     let x = if a.const_term() == 1 {
//         Integer::from(1)
//     } else if a.const_term() == -1 {
//         Integer::from(-1)
//     } else {
//         panic!("inv_trunc: non-invertible constant term");
//     };

//     if m == 1 {
//         conv(c, x);
//         return;
//     }

//     let mut _E = Vec::new();
//     _E.push(m);
//     let mut m = m;
//     while m > 1 {
//         m = (m + 1) / 2;
//         _E.push(m);
//     }

//     let _L = _E.len();

//     let mut g = ZZX::new();
//     g.set_length(_E[0] as usize);

//     let mut g0 = ZZX::new();
//     g0.set_length(_E[0] as usize);

//     let mut g1 = ZZX::new();
//     g1.set_length(((3 * _E[0] + 1) / 2) as usize);

//     let mut g2 = ZZX::new();
//     g2.set_length(_E[0] as usize);

//     conv(g, x);

//     for i in (1.._L).rev() {
//         // lift from _E[i] to _E[i - 1]
//         let k = _E[i];
//         let l = _E[i - 1] - _E[i];

//         _trunc(&mut g0, a, (k + l) as usize);
//         mul(&mut g1, &g0, &g);
//         _right_shift(&mut g1, &g1, k);
//         _trunc(&mut g, &g1, l as usize);

//         mul(&mut g2, &g1, &g);
//         _trunc(&mut g2, &g2, l as usize);
//         _left_shift(&mut g2, &g2, k);

//         sub(&mut g, &g, &g2);
//     }

//     c = g;

// }
