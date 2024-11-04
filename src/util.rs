use rug::Integer;

/// Custom clone of NTL::ZZX 
#[derive(Clone, Debug, Default)]
pub struct ZZX {
    coeffs: Vec<Integer>,
}

impl ZZX {
    pub fn new() -> Self {
        ZZX { coeffs: vec![] }
    }

    pub fn new_with_u64(n: u64) -> Self {
        ZZX { coeffs: vec![Integer::from(n)] }
    }

    pub fn new_with_integer(n: Integer) -> Self {
        ZZX { coeffs: vec![n] }
    }

    /*************************
     * Some utility functions
     *************************/
    pub fn set_coeff<T>(&mut self, i: usize, n: Option<T>) where T: Into<Integer> {
        self.coeffs[i] = n.map(Into::into).unwrap_or(Integer::from(1));
    }

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
    pub fn set_length(&mut self, len: usize) {
        self.coeffs.resize(len, Integer::from(0));
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
    todo!("impl `void Mul(ZZX& c, const ZZX& a, const ZZX& b)` func");
}

pub fn rem(c: &mut ZZX, a: &ZZX, b: &ZZX) {
    todo!("impl `void Rem(ZZX& r, const ZZX& a, const ZZX& b)` func");
}
