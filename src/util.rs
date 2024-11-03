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
}
