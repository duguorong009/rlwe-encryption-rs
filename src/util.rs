use std::ops::{
    Add, AddAssign, Div, Index, IndexMut, Mul, MulAssign, Neg, Rem, Shl, ShlAssign, Shr, ShrAssign,
    Sub, SubAssign,
};

use rug::{
    ops::{NegAssign, Pow},
    Complete, Integer,
};

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
        let n = n.into();
        if n == 0 {
            ZZX::new()
        } else {
            ZZX {
                coeffs: vec![n.into()],
            }
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

    pub fn set_max_length(&mut self, len: usize) {
        if self.coeffs.len() < len {
            self.coeffs.reserve_exact(len - self.coeffs.len());
        }
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
        c.set_coeff::<i64>(1, None);
        c
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
            res = res.max(self.coeffs[i].significant_bits());
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

    /*************************
     * Comparison functions
     *************************/
    pub fn is_zero(&self) -> bool {
        self.coeffs.is_empty()
    }

    pub fn is_one(&self) -> bool {
        self.coeffs.len() == 1 && self.coeffs[0] == 1
    }
}

impl Index<usize> for ZZX {
    type Output = Integer;

    fn index(&self, index: usize) -> &Self::Output {
        if index >= self.coeffs.len() {
            panic!("index out of range");
        } else {
            &self.coeffs[index]
        }
    }
}

impl IndexMut<usize> for ZZX {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        if index >= self.coeffs.len() {
            self.set_length(index + 1);
        }

        &mut self.coeffs[index]
    }
}

// Modular arithmetic -- f must be monic, and other args
// must have degree less than that of f
fn _mulmod(x: &mut ZZX, a: &ZZX, b: &ZZX, f: &ZZX) {
    if a.deg() >= f.deg() || b.deg() >= f.deg() || f.deg() == 0 || f.lead_coeff() != 1 {
        panic!("MulMod: bad args");
    }

    let mut t = ZZX::new();
    mul(&mut t, a, b);
    rem(x, &t, f);
}

pub fn mulmod(a: &ZZX, b: &ZZX, f: &ZZX) -> ZZX {
    let mut x = ZZX::new();
    _mulmod(&mut x, a, b, f);
    x
}

fn _sqrmod(x: &mut ZZX, a: &ZZX, f: &ZZX) {
    if a.deg() >= f.deg() || f.deg() == 0 || f.lead_coeff() != 1 {
        panic!("SqrMod: bad args");
    }
    let mut t = ZZX::new();
    _sqr(&mut t, a);
    rem(x, &t, f);
}

pub fn sqrmod(a: &ZZX, f: &ZZX) -> ZZX {
    let mut x = ZZX::new();
    _sqrmod(&mut x, a, f);
    x
}

fn _mulbyxmod_aux(h: &mut ZZX, a: &ZZX, f: &ZZX) {
    let n = f.deg();
    let m = a.deg();

    if m >= n || n == 0 || f.lead_coeff() != 1 {
        panic!("MulByXMod: bad args");
    }

    if m < 0 {
        h.clear();
        return;
    }

    if m < n - 1 {
        h.set_length(m as usize + 2);
        for i in (1..m).rev() {
            h.coeffs[i as usize] = a.coeffs[i as usize - 1].clone();
            h.coeffs[0] = Integer::new();
        }
    } else {
        h.set_length(n as usize);
        let z = -a.coeffs[n as usize - 1].clone();
        for i in (1..n).rev() {
            let t = z.clone() * f.coeffs[i as usize].clone();
            h.coeffs[i as usize] = a.coeffs[i as usize - 1].clone() * t;
        }
        h.coeffs[0] = z * f.coeffs[0].clone();
        h.normalize();
    }
}

fn _mulbyxmod(h: &mut ZZX, a: &ZZX, f: &ZZX) {
    _mulbyxmod_aux(h, a, f);
}

pub fn mulbyxmod(a: &ZZX, f: &ZZX) -> ZZX {
    let mut x = ZZX::new();
    _mulbyxmod(&mut x, a, f);
    x
}

/// x = a * b
pub fn mul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    if a.is_zero() || b.is_zero() {
        c.set_length(0);
        return;
    }
    if a == b {
        _sqr(c, a);
        return;
    }

    // let maxa: u32 = a.max_size();
    // let maxb: u32 = b.max_size();

    // let k = maxa.min(maxb);
    // let s = a.deg().min(b.deg()) + 1;

    // if s == 1 || (k == 1 && s < 40) || (k == 2 && s < 20) || (k == 3 && s < 10) {
    //     plain_mul(c, a, b);
    // } else if s < 80 || (k < 30 && s < 150) {
    //     kar_mul(c, a, b);
    // } else if choose_ss(a.deg(), a.max_bits(), b.deg(), b.max_bits()) {
    //     ss_mul(c, a, b);
    // } else {
    //     hom_mul(c, a, b);
    // }

    // TODO: uncomment the above.
    plain_mul(c, a, b);
}

impl Mul for ZZX {
    type Output = ZZX;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut res = ZZX::new();
        mul(&mut res, &self, &rhs);
        res
    }
}

impl MulAssign for ZZX {
    fn mul_assign(&mut self, rhs: Self) {
        let mut output = ZZX::new();
        mul(&mut output, &self, &rhs);
        *self = output;
    }
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
        const_rem(r, a, &b.const_term());
    } else if b.lead_coeff() == 1 {
        pseudo_rem(r, a, b);
    } else if b.lead_coeff() == -1 {
        let mut b1 = ZZX::new();
        negate(&mut b1, b);
        pseudo_rem(r, a, &b1);
    } else if divide_(a, b) {
        r.coeffs = vec![Integer::from(0)];
    } else {
        let mut r1 = ZZX::new();
        pseudo_rem(&mut r1, a, b);
        let m = Integer::from(b.lead_coeff().pow(da as u32 - db as u32 + 1));
        if !divide_with_integer(r, &r1, &m) {
            panic!("rem: remainder not defined over ZZ");
        }
    }
}

fn const_rem(r: &mut ZZX, a: &ZZX, b: &Integer) {
    if b == &0 {
        panic!("const_rem: division by zero");
    }
    r.coeffs = vec![Integer::from(0)];
}

fn plain_pseudo_rem(r: &mut ZZX, a: &ZZX, b: &ZZX) {
    let mut q = ZZX::new();
    plain_pseudo_div_rem(&mut q, r, a, b);
}

fn plain_pseudo_div(q: &mut ZZX, a: &ZZX, b: &ZZX) {
    let mut r = ZZX::new();
    plain_pseudo_div_rem(q, &mut r, a, b);
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
                t = Integer::from(&t * &lc);
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

fn pseudo_div(q: &mut ZZX, a: &ZZX, b: &ZZX) {
    plain_pseudo_div(q, a, b);
}

fn pseudo_rem(r: &mut ZZX, a: &ZZX, b: &ZZX) {
    plain_pseudo_rem(r, a, b);
}

fn pseudo_div_rem(q: &mut ZZX, r: &mut ZZX, a: &ZZX, b: &ZZX) {
    plain_pseudo_div_rem(q, r, a, b);
}

fn div_rem(q: &mut ZZX, r: &mut ZZX, a: &ZZX, b: &ZZX) {
    let da = a.deg();
    let db = b.deg();

    if db < 0 {
        panic!("div_rem: division by zero");
    }

    if da < db {
        r.coeffs = a.coeffs.clone();
        q.clear();
    } else if db == 0 {
        const_div_rem(q, r, a, &b.const_term());
    } else if b.lead_coeff() == 1 {
        pseudo_div_rem(q, r, a, b);
    } else if b.lead_coeff() == -1 {
        let mut b1 = ZZX::new();
        negate(&mut b1, b);
        pseudo_div_rem(q, r, a, &b1);
        negate(q, &q.clone());
    } else if divide(q, a, b) {
        r.clear();
    } else {
        let mut q1 = ZZX::new();
        let mut r1 = ZZX::new();
        let mut m = Integer::new();
        pseudo_div_rem(&mut q1, &mut r1, a, b);
        m = b.lead_coeff().pow(da as u32 - db as u32 + 1);
        if !divide_with_integer(q, &mut q1, &m) {
            panic!("div_rem: quotient not defined over ZZ");
        }
        if !divide_with_integer(r, &mut r1, &m) {
            panic!("div_rem: remainder not defined over ZZ");
        }
    }
}

fn const_div_rem(q: &mut ZZX, r: &mut ZZX, a: &ZZX, b: &Integer) {
    if b == &0 {
        panic!("const_div_rem: division by zero");
    }
    if !divide_with_integer(q, a, b) {
        panic!("const_div_rem: quotient not defined over ZZ");
    }
    r.clear();
}

fn div(q: &mut ZZX, a: &ZZX, b: &ZZX) {
    let da = a.deg();
    let db = b.deg();

    if db < 0 {
        panic!("div: division by zero");
    }

    if da < db {
        q.clear();
    } else if db == 0 {
        div_with_integer(q, a, &b.const_term());
    } else if b.lead_coeff() == 1 {
        pseudo_div(q, a, b);
    } else if b.lead_coeff() == -1 {
        let mut b1 = ZZX::new();
        negate(&mut b1, b);
        pseudo_div(q, a, &b1);
        negate(q, &q.clone());
    } else if divide(q, a, b) {
        // nothing to do
    } else {
        let mut q1 = ZZX::new();
        let mut m = Integer::new();
        pseudo_div(&mut q1, a, b);
        m = b.lead_coeff().pow(da as u32 - db as u32 + 1);
        if !divide_with_integer(q, &mut q1, &m) {
            panic!("div: quotient not defined over ZZ");
        }
    }
}

fn div_with_integer(q: &mut ZZX, a: &ZZX, b: &Integer) {
    if b == &0 {
        panic!("div: division by zero");
    }
    if !divide_with_integer(q, a, b) {
        panic!("div: quotient not defined over ZZ");
    }
}

impl Div for ZZX {
    type Output = ZZX;

    fn div(self, rhs: Self) -> Self::Output {
        let mut res = ZZX::new();
        div(&mut res, &self, &rhs);
        res
    }
}

impl Div<Integer> for ZZX {
    type Output = ZZX;

    fn div(self, rhs: Integer) -> Self::Output {
        let mut res = ZZX::new();
        div_with_integer(&mut res, &self, &rhs);
        res
    }
}

impl Rem for ZZX {
    type Output = ZZX;

    fn rem(self, rhs: Self) -> Self::Output {
        let mut res = ZZX::new();
        rem(&mut res, &self, &rhs);
        res
    }
}

impl Rem<Integer> for ZZX {
    type Output = ZZX;

    fn rem(self, rhs: Integer) -> Self::Output {
        let mut res = ZZX::new();
        const_rem(&mut res, &self, &rhs);
        res
    }
}

// if b | a, sets q = a / b and returns true; otherwise returns false
fn divide(q: &mut ZZX, a: &ZZX, b: &ZZX) -> bool {
    let da = a.deg();
    let db = b.deg();

    // if db <= 8 || da - db <= 8 {
    //     plain_divide(q, a, b)
    // } else {
    //     hom_divide(q, a, b)
    // }

    // TODO: uncomment the above.
    plain_divide(q, a, b)
}

fn divide_(a: &ZZX, b: &ZZX) -> bool {
    let da = a.deg();
    let db = b.deg();

    if db <= 8 || da - db <= 8 {
        plain_divide_(a, b)
    } else {
        hom_divide_(a, b)
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
        return divide_with_integer(qq, aa, &bb.const_term());
    }

    let da = aa.deg();
    let db = bb.deg();

    if da < db {
        return false; // 0
    }

    let mut ca = Integer::from(0);
    let mut cb = Integer::from(0);
    _content(&mut ca, aa);
    _content(&mut cb, bb);

    let (cq, r) = ca.div_rem_ref(&cb).complete();
    if !r.is_zero() {
        return false; // 0
    }

    let mut a = ZZX::new();
    let mut b = ZZX::new();
    divide(&mut a, aa, bb);
    divide_with_integer(&mut b, bb, &cb);

    if !a.lead_coeff().is_divisible(&b.lead_coeff()) {
        return false; // 0
    }

    if !a.const_term().is_divisible(&b.const_term()) {
        return false; // 0
    }

    let coeff_bnd = a.max_bits() as i64
        + ((Integer::from(da + 1).significant_bits() + 1) / 2) as i64
        + (da - db);

    let bp = b.coeffs.clone();
    let lc = bp[db as usize].clone();

    let lc_is_one = lc == 1;

    let mut xp = a.coeffs.clone();

    let dq = da - db;
    let mut q = ZZX::new();
    q.set_length(dq as usize + 1);

    let mut t = Integer::new();
    for i in (0..=dq).rev() {
        if !lc_is_one {
            let (q, r) = &xp[i as usize + db as usize].div_rem_ref(&lc).complete();
            t = q.clone();
            if !r.is_zero() {
                return false;
            }
        } else {
            t = xp[i as usize + db as usize].clone();
        }

        if t.significant_bits() as i64 > coeff_bnd {
            return false;
        }

        q.coeffs[i as usize] = t.clone();
        for j in (0..db).rev() {
            let s = Integer::from(&t * &bp[j as usize]);
            xp[i as usize + j as usize] = Integer::from(&xp[i as usize + j as usize] - &s);
        }
    }

    for i in 0..db {
        if !xp[i as usize].is_zero() {
            return false;
        }
    }

    mul_with_integer(qq, &q, &cq);
    true
}

fn plain_divide_(a: &ZZX, b: &ZZX) -> bool {
    if b.deg() == 0 {
        divide_(a, b)
    } else {
        let mut q = ZZX::new();
        plain_divide(&mut q, a, b)
    }
}

fn divide_with_integer(q: &mut ZZX, a: &ZZX, b: &Integer) -> bool {
    if b == &0 {
        if a.is_zero() {
            q.clear();
            return true; // 1
        } else {
            return false; // 0
        }
    }

    if b == &1 {
        q.coeffs = a.coeffs.clone();
        return true; // 1
    }

    if b == &-1 {
        negate(q, a);
        return true; // 1
    }

    let n = a.coeffs.len();
    let mut res: Vec<Integer> = Vec::with_capacity(n);
    for i in 0..n {
        let (q, r) = a.coeffs[i].clone().div_rem(b.clone());
        if !r.is_zero() {
            return false;
        }
        res[i] = q;
    }

    q.coeffs = res;
    true
}

fn divide_with_integer_(a: &ZZX, b: &Integer) -> bool {
    if b == &0 {
        return a.is_zero();
    }
    if b == &1 || b == &-1 {
        return true;
    }

    let n = a.coeffs.len();
    for i in 0..n {
        if !a.coeffs[i].is_divisible(b) {
            return false;
        }
    }

    true
}

fn mul_with_integer(x: &mut ZZX, a: &ZZX, b: &Integer) {
    if b == &0 {
        x.clear();
        return;
    }

    let t = b.clone();
    let da = a.deg();
    x.set_length(da as usize + 1);

    for i in 0..da + 1 {
        x.coeffs[i as usize] = Integer::from(&a.coeffs[i as usize] * &t);
    }
}

fn hom_divide(q: &mut ZZX, a: &ZZX, b: &ZZX) -> bool {
    if b.is_zero() {
        if a.is_zero() {
            q.clear();
            return true; // 1
        } else {
            return false; // 0
        }
    }

    if a.is_zero() {
        q.clear();
        return true; // 1
    }

    if b.deg() == 0 {
        return divide_with_integer(q, a, &b.const_term());
    }

    if a.deg() < b.deg() {
        return false; // 0
    }

    let mut ca = Integer::new();
    let mut cb = Integer::new();
    _content(&mut ca, a);
    _content(&mut cb, b);

    let (_, r) = ca.div_rem_ref(&cb).complete();
    if !r.is_zero() {
        return false; // 0
    }

    let mut aa = ZZX::new();
    let mut bb = ZZX::new();
    divide_with_integer(&mut aa, a, &ca);
    divide_with_integer(&mut bb, b, &cb);

    if !aa.lead_coeff().is_divisible(&bb.lead_coeff()) {
        return false; // 0
    }

    if !aa.const_term().is_divisible(&bb.const_term()) {
        return false; // 0
    }

    todo!()
}

fn hom_divide_(a: &ZZX, b: &ZZX) -> bool {
    if b.deg() == 0 {
        divide_with_integer_(a, &b.const_term())
    } else {
        let mut q = ZZX::new();
        hom_divide(&mut q, a, b)
    }
}

/// c = content of f, sign(c) = sign(f.lead_coeff())
fn _content(c: &mut Integer, f: &ZZX) {
    let mut res = Integer::from(0);
    for i in 0..f.coeffs.len() {
        res = res.gcd(&f.coeffs[i]);
        if res == 1 {
            break;
        }
    }
    if f.lead_coeff().is_negative() {
        res = -res;
    }
    *c = res;
}

fn content(f: &ZZX) -> Integer {
    let mut res = Integer::from(0);
    _content(&mut res, f);
    res
}

/// x = a ^ 2
fn _sqr(c: &mut ZZX, a: &ZZX) {
    if a.is_zero() {
        c.set_length(0);
        return;
    }

    let maxa = a.max_size();
    let k = maxa;
    let s = a.deg() + 1;

    // if s == 1 || (k == 1 && s < 50) || (k == 2 && s < 25) || (k == 3 && s < 25) || (k == 4 && s < 10) {
    //     plain_sqr(c, a);
    // } else if s < 80 || (k < 30 && s < 150) {
    //     kar_sqr(c, a);
    // } else if choose_ss(a.deg(), a.max_bits(), a.deg(), a.max_bits()) {
    //     ss_sqr(c, a);
    // } else {
    //     hom_sqr(c, a);
    // }

    // TODO: uncomment the above.
    plain_sqr(c, a);
}

/// x = a ^ 2
pub fn sqr(a: &ZZX) -> ZZX {
    let mut c = ZZX::new();
    _sqr(&mut c, a);
    c
}

fn plain_mul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    if a == b {
        plain_sqr(c, a);
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

fn kar_mul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    if a.is_zero() || b.is_zero() {
        c.clear();
        return;
    }

    if a == b {
        kar_sqr(c, a);
        return;
    }

    let sa = a.coeffs.len();
    let sb = b.coeffs.len();

    c.set_length(sa + sb - 1);

    let maxa = a.max_bits();
    let maxb = b.max_bits();
    let xover = 2;

    if sa < xover || sb < xover {
        // _plain_mul_in_(&c.coeffs, &a.coeffs, sa, &b.coeffs, sb);
        todo!()
    } else {
        /* karatsuba */
        todo!()
    }

    c.normalize();
}

fn ss_mul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    todo!("impl `void ss_mul(ZZX& c, const ZZX& a, const ZZX& b)` func");
}

fn hom_mul(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    todo!("impl `void hom_mul(ZZX& c, const ZZX& a, const ZZX& b)` func");
}

fn choose_ss(a_deg: i64, a_max_bits: u32, b_deg: i64, b_max_bits: u32) -> bool {
    todo!("impl `bool choose_ss(int a_deg, int a_max_bits, int b_deg, int b_max_bits)` func");
}

fn plain_sqr(c: &mut ZZX, a: &ZZX) {
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

fn kar_sqr(c: &mut ZZX, a: &ZZX) {
    todo!("impl `void kar_sqr(ZZX& c, const ZZX& a)` func");
}

fn ss_sqr(c: &mut ZZX, a: &ZZX) {
    todo!("impl `void ss_sqr(ZZX& c, const ZZX& a)` func");
}

fn hom_sqr(c: &mut ZZX, a: &ZZX) {
    todo!("impl `void hom_sqr(ZZX& c, const ZZX& a)` func");
}

/// x = a % X^m
fn _trunc(x: &mut ZZX, a: &ZZX, m: i64) {
    if m < 0 {
        panic!("trunc: m < 0");
    }

    let m = m as usize;

    let n = m.min(a.coeffs.len());
    x.set_length(n);

    for i in 0..n {
        x.coeffs[i] = a.coeffs[i].clone();
    }
    x.normalize();
}

pub fn trunc(a: &ZZX, m: i64) -> ZZX {
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

    x.set_length(da - n + 1);

    for i in 0..=(da - n) {
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
}

pub fn left_shift(a: &ZZX, n: i64) -> ZZX {
    let mut x = ZZX::new();
    _left_shift(&mut x, a, n);
    x
}

impl Shl<i64> for ZZX {
    type Output = ZZX;

    fn shl(self, shift: i64) -> Self::Output {
        left_shift(&self, shift)
    }
}

impl ShlAssign<i64> for ZZX {
    fn shl_assign(&mut self, shift: i64) {
        *self = left_shift(self, shift);
    }
}

impl Shr<i64> for ZZX {
    type Output = ZZX;

    fn shr(self, shift: i64) -> Self::Output {
        right_shift(&self, shift)
    }
}

impl ShrAssign<i64> for ZZX {
    fn shr_assign(&mut self, shift: i64) {
        *self = right_shift(self, shift);
    }
}

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

pub fn add(x: &mut ZZX, a: &ZZX, b: &ZZX) {
    let da = a.deg();
    let db = b.deg();
    let maxab = da.max(db);

    x.set_length(maxab as usize + 1);

    for i in 0..maxab as usize + 1 {
        let a = if i <= da as usize {
            a.coeffs[i].clone()
        } else {
            Integer::from(0)
        };
        let b = if i <= db as usize {
            b.coeffs[i].clone()
        } else {
            Integer::from(0)
        };
        x.coeffs[i] = a + b;
    }

    x.normalize();
}

pub fn sub(x: &mut ZZX, a: &ZZX, b: &ZZX) {
    let da = a.deg();
    let db = b.deg();
    let maxab = da.max(db);

    x.set_length(maxab as usize + 1);

    for i in 0..maxab as usize + 1 {
        let a = if i <= da as usize {
            a.coeffs[i].clone()
        } else {
            Integer::from(0)
        };
        let b = if i <= db as usize {
            b.coeffs[i].clone()
        } else {
            Integer::from(0)
        };
        x.coeffs[i] = a - b;
    }

    x.normalize();
}

fn negate(q: &mut ZZX, a: &ZZX) {
    q.coeffs = a.coeffs.clone();
    for i in 0..q.coeffs.len() {
        q.coeffs[i].neg_assign();
    }
}

impl Add for ZZX {
    type Output = ZZX;

    fn add(self, rhs: Self) -> Self::Output {
        let mut output = ZZX::new();
        add(&mut output, &self, &rhs);
        output
    }
}

impl Add<&ZZX> for ZZX {
    type Output = ZZX;

    fn add(self, rhs: &ZZX) -> Self::Output {
        let mut output = ZZX::new();
        add(&mut output, &self, rhs);
        output
    }
}

impl Add<ZZX> for &ZZX {
    type Output = ZZX;

    fn add(self, rhs: ZZX) -> Self::Output {
        let mut output = ZZX::new();
        add(&mut output, self, &rhs);
        output
    }
}

impl Add<&ZZX> for &ZZX {
    type Output = ZZX;

    fn add(self, rhs: &ZZX) -> Self::Output {
        let mut output = ZZX::new();
        add(&mut output, self, rhs);
        output
    }
}

impl AddAssign for ZZX {
    fn add_assign(&mut self, rhs: Self) {
        let mut tmp = ZZX::new();
        add(&mut tmp, self, &rhs);
        *self = tmp;
    }
}

impl Sub for ZZX {
    type Output = ZZX;

    fn sub(self, rhs: Self) -> Self::Output {
        let mut output = ZZX::new();
        sub(&mut output, &self, &rhs);
        output
    }
}

impl Sub<&ZZX> for ZZX {
    type Output = ZZX;

    fn sub(self, rhs: &ZZX) -> Self::Output {
        let mut output = ZZX::new();
        sub(&mut output, &self, rhs);
        output
    }
}

impl Sub<ZZX> for &ZZX {
    type Output = ZZX;

    fn sub(self, rhs: ZZX) -> Self::Output {
        let mut output = ZZX::new();
        sub(&mut output, self, &rhs);
        output
    }
}

impl Sub<&ZZX> for &ZZX {
    type Output = ZZX;

    fn sub(self, rhs: &ZZX) -> Self::Output {
        let mut output = ZZX::new();
        sub(&mut output, self, rhs);
        output
    }
}

impl SubAssign for ZZX {
    fn sub_assign(&mut self, rhs: Self) {
        let mut tmp = ZZX::new();
        sub(&mut tmp, self, &rhs);
        *self = tmp;
    }
}

impl Neg for ZZX {
    type Output = ZZX;

    fn neg(self) -> Self::Output {
        let mut output = ZZX::new();
        negate(&mut output, &self);
        output
    }
}

impl NegAssign for ZZX {
    fn neg_assign(&mut self) {
        let mut output = ZZX::new();
        negate(&mut output, self);
        *self = output;
    }
}

/// x = a^{-1} % X^m
fn _inv_trunc(c: &mut ZZX, a: &ZZX, e: i64) {
    if e < 0 {
        panic!("inv_trunc: m < 0");
    }

    if e == 0 {
        c.clear();
        return;
    }

    // TODO: add ntl overflow check `NTL_OVERFLOW(e, 1, 0)`

    newton_inv_trunc(c, a, e);
}

pub fn inv_trunc(a: &ZZX, e: i64) -> ZZX {
    let mut x = ZZX::new();
    _inv_trunc(&mut x, a, e);
    x
}

fn newton_inv_trunc(c: &mut ZZX, a: &ZZX, e: i64) {
    let x = if a.const_term() == 1 {
        Integer::from(1)
    } else if a.const_term() == -1 {
        Integer::from(-1)
    } else {
        panic!("inv_trunc: non-invertible constant term");
    };

    if e == 1 {
        c.coeffs = ZZX::new_with_val(x).coeffs;
        return;
    }

    let mut _e = Vec::new();
    _e.push(e);
    let mut e = e;
    while e > 1 {
        e = (e + 1) / 2;
        _e.push(e);
    }

    let _l = _e.len();

    let mut g = ZZX::new();
    g.set_max_length(_e[0] as usize);

    let mut g0 = ZZX::new();
    g0.set_max_length(_e[0] as usize);

    let mut g1 = ZZX::new();
    g1.set_max_length(((3 * _e[0] + 1) / 2) as usize);

    let mut g2 = ZZX::new();
    g2.set_max_length(_e[0] as usize);

    let mut g = ZZX::new_with_val(x);

    for i in (1.._l).rev() {
        // lift from _e[i] to _e[i - 1]
        let k = _e[i];
        let l = _e[i - 1] - _e[i];

        g0 = trunc(a, k + l);
        g1 = g0 * g.clone();
        g1 = right_shift(&g1, k);
        g1 = trunc(&g1, l);

        g2 = g1 * g.clone();
        g2 = trunc(&g2, l);
        g2 = left_shift(&g2, k);

        g = g - g2;
    }

    c.coeffs = g.coeffs;
}

/// x = a * b % X^n
fn _mul_trunc(x: &mut ZZX, a: &ZZX, b: &ZZX, n: i64) {
    let mut t = ZZX::new();
    mul(&mut t, a, b);
    _trunc(x, &t, n);
}

pub fn mul_trunc(a: &ZZX, b: &ZZX, n: i64) -> ZZX {
    let mut x = ZZX::new();
    _mul_trunc(&mut x, a, b, n);
    x
}

/// x = a^2 % X^n
fn _sqr_trunc(x: &mut ZZX, a: &ZZX, n: i64) {
    let mut t = ZZX::new();
    _sqr(&mut t, a);
    _trunc(x, &t, n);
}

pub fn sqr_trunc(a: &ZZX, n: i64) -> ZZX {
    let mut x = ZZX::new();
    _sqr_trunc(&mut x, a, n);
    x
}

fn _reverse(x: &mut ZZX, a: &ZZX, hi: i64) {
    if hi < 0 {
        x.clear();
        return;
    }

    // TODO: check NTL_OVERFLOW(hi, 1, 0)

    if x == a {
        let mut tmp = ZZX::new();
        _copy_reverse(&mut tmp, a, hi as usize);
        x.coeffs = tmp.coeffs;
    } else {
        _copy_reverse(x, a, hi as usize);
    }
}

/// x[0..hi] = reverse(a[0..hi]), with zero fill
/// input may not alias output
fn _copy_reverse(x: &mut ZZX, a: &ZZX, hi: usize) {
    let n = hi + 1;
    let m = a.coeffs.len();

    x.set_length(n);

    for i in 0..n {
        let j = hi - i;
        if j >= m {
            x.coeffs[i] = Integer::from(0);
        } else {
            x.coeffs[i] = a.coeffs[j].clone();
        }
    }

    x.normalize();
}

pub fn reverse(a: &ZZX, hi: i64) -> ZZX {
    let mut x = ZZX::new();
    _reverse(&mut x, a, hi);
    x
}
