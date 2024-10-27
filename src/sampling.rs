use rug::{float::Constant, ops::Pow, Float};

pub const ROWS: usize = 128;
pub const COLS: usize = 40;

pub struct Sampling {
    p: Vec<Vec<i32>>,
    begin: Vec<i32>,
    precision: i32,
    tailcut: f32,
    sigma: Float,
    c: Float,
}

impl Sampling {
    pub fn new(precision: i32, tailcut: f32, sigma: Float, center: Float) -> Self {
        let sampling = Self {
            p: vec![],
            begin: vec![],
            precision,
            tailcut,
            sigma,
            c: center,
        };
        // RR::set_precision(precision as u64); 
        sampling.build_probability_matrix();
        sampling
    }

    pub fn knuth_yao(&self) -> i32 {
        todo!()
    }

    fn build_probability_matrix(&self) {
        todo!()
    }

    fn probability(&self, x: Float, sigma: Float, c: Float) -> Float {
        let pi = Float::with_val(32, Constant::Pi);
        let mut s: Float = sigma * (Float::with_val(32, 2) * pi).sqrt();
        let mut over_s: Float = 1 / s;

        if x == 0 {
            return over_s;
        }
        // over_s * exp(-(power((x - c) / sigma, 2)) / 2.0)
        todo!()
    }

    //  Method for computing the binary expansion of a given probability in [0, 1] 
    fn binary_expansion(
        &self,
        aux_p: &mut Vec<Vec<i32>>,
        mut probability: Float,
        precision: u64,
        index: usize,
    ) {
        let mut i = -1;
        let mut j: usize = 0;

        while probability > 0 && ((j as u64) < precision) {
            let pow = Float::with_val(32, 2).pow(i); // 2 ^ {i}
            if pow <= probability {
                aux_p[j][index] = 1;
                probability -= pow;
            } else {
                aux_p[j][index] = 0;
            }
            i -= 1;
            j += 1;
        }
    }

    // bit = 0 then return a
    fn select(a: i32, b: i32, bit: u32) -> i32 {
        let mask = -(bit as i32);
        (mask & (a ^ b)) ^ a
    }
}
