use rug::{
    float::{Constant, Round},
    ops::{DivAssignRound, Pow},
    Float,
};

use crate::util::randombits_i64;

#[derive(Debug, Clone)]
pub struct Sampling {
    p: Vec<Vec<i64>>,
    begin: Vec<i64>,
    precision: u32,
    tailcut: f32,
    sigma: Float,
    c: Float,
}

impl Sampling {
    pub fn new(precision: u32, tailcut: f32, sigma: Float, center: Float) -> Self {
        let mut sampling = Self {
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

    // Knuth-Yao algorithm to obtain a sample from the discrete gaussian
    pub fn knuth_yao(&self) -> i64 {
        let bound = (self.tailcut * self.sigma.clone().to_f32()).round() as usize;
        let center = self.c.to_f32().round() as i64;
        let p_num_rows = self.p.len(); // precision
        let p_num_cols = self.p[0].len();
        
        let mut random_bits: Vec<i64> = vec![0; p_num_rows];
        
        let length = 64; // sizeof(unsigned long)*8 // 64 bits
        
        let mut index = 0;
        for _ in 0..(p_num_rows / length + 1) {
            let mut r: u64 = rand::random::<u64>(); // RandomWord(); // It returns a word filled with pseudo-random bits
            let mut j = 0;
            while j < length && index < p_num_rows {
                random_bits[index] = (r & 1) as i64; // getting the least significant bit
                
                j += 1;
                index += 1;
                r = r >> 1;
            }
        }
        
        let mut d = 0; // distance
        let invalid_sample = bound + 1;
        let signal = 1 - 2 * randombits_i64(1) as i64; // Sample a random signal s
        let mut hit = false;

        let mut s = 0;
        for row in 0..p_num_rows {
            d = 2 * d + random_bits[row]; // Distance calculus
            for col in self.begin[row] as usize..p_num_cols {
                d = d - self.p[row][col];

                let enable = d == -1;

                // when enable & !hit becomes 1, "col" is added to "S";
                // e.g. enable = 1 and hit = 0
                s += select(invalid_sample as i64, col as i64, enable & !hit);
                hit |= enable & !hit;
            }
        }

        // Note: the "col" value is in [0, bound]. So, the invalid sample must be greater than bound.
        s = s % invalid_sample as i64;
        s = s - bound as i64 + center;
        s *= signal;

        s
    }

    fn build_probability_matrix(&mut self) {
        // RR::set_precision(to_long(self.precision));

        let mut aux_p: Vec<Vec<i64>> = vec![];
        let mut aux_begin: Vec<i64> = vec![];

        // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
        let mut prob_of_x: Vec<Float> = vec![];

        let bound = (self.tailcut * self.sigma.clone().to_f32()).round() as usize;
        prob_of_x.resize_with(bound + 1, || Float::with_val(self.precision, 0));
        aux_p.resize_with(self.precision as usize, || vec![]);

        for i in 0..aux_p.len() {
            aux_p[i].resize_with(bound + 1, || 0);
        }

        for x in (1..=bound).rev() {
            prob_of_x[bound - x] = self.probability(
                Float::with_val(self.precision, x) + self.c.clone(),
                self.sigma.clone(),
                self.c.clone(),
            );
        }

        prob_of_x[bound] = self.probability(
            Float::with_val(self.precision, 0) + self.c.clone(),
            self.sigma.clone(),
            self.c.clone(),
        );
        prob_of_x[bound].div_assign_round(Float::with_val(self.precision, 2), Round::Nearest);

        let mut i = -1;
        for j in 0..self.precision as usize {
            let pow: Float = Float::with_val(self.precision, 2).pow(i); // 2^{i}
            i -= 1;
            for x in (0..=bound).rev() {
                aux_p[j][bound - x] = 0;
                if prob_of_x[bound - x] >= pow.clone() {
                    aux_p[j][bound - x] = 1;
                    prob_of_x[bound - x] -= pow.clone();
                }
            }
        }

        self.p = aux_p;

        let p_num_cols = self.p[0].len();
        let p_num_rows = self.p.len();

        aux_begin.resize_with(p_num_rows, || 0);

        // computing in which position the non-zero values in P start and end
        for i in 0..p_num_rows {
            aux_begin[i] = p_num_cols as i64 - 1;

            for j in 0..p_num_cols {
                if self.p[i][j] == 1 {
                    aux_begin[i] = j as i64;
                    break;
                }
            }
        }

        self.begin = aux_begin;
    }

    fn probability(&self, x: Float, sigma: Float, c: Float) -> Float {
        let pi = Float::with_val(self.precision, Constant::Pi);
        let s: Float = sigma.clone() * (Float::with_val(self.precision, 2) * pi).sqrt();
        let over_s: Float = 1 / s;

        if x == 0 {
            return over_s;
        }

        // over_s * exp(-(power((x - c) / sigma, 2)) / 2.0)
        let tmp = (x - c) / sigma;
        let tmp2 = tmp.clone() * tmp;
        let tmp3 = -(tmp2 / Float::with_val(self.precision, 2.0));
        let tmp4 = tmp3.exp();
        over_s * tmp4
    }

    //  Method for computing the binary expansion of a given probability in [0, 1]
    fn binary_expansion(
        &self,
        aux_p: &mut Vec<Vec<i64>>,
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
}

// bit = 0 then return a
fn select(a: i64, b: i64, bit: bool) -> i64 {
    if bit {
        b
    } else {
        a
    }
}
