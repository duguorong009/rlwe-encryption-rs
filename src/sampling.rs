use rug::{float::{Constant, Round}, ops::{DivAssignRound, Pow}, Float};

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

    pub fn knuth_yao(&self) -> i32 {
        todo!()
    }

    fn build_probability_matrix(&mut self) {
        // RR::set_precision(to_long(self.precision));

        let mut aux_p: Vec<Vec<i32>> = vec![];
        let mut aux_begin: Vec<i32> = vec![];

        // The random variable consists of elements in [c-tailcut*sigma, c+tailcut*sigma]
        let mut prob_of_x: Vec<Float> = vec![];

        let bound = (self.tailcut * self.sigma.clone().to_f32()).round() as usize;
        prob_of_x.reserve_exact(bound + 1);
        aux_p.reserve_exact(self.precision as usize);

        for i in 0..aux_p.len() {
            aux_p[i].reserve_exact(bound + 1);
        }

        for x in (1..=bound).rev() {
            prob_of_x[bound - x] = self.probability(Float::with_val(32, x) + self.c.clone(), self.sigma.clone(), self.c.clone());
        }

        prob_of_x[bound] = self.probability(Float::with_val(32, 0) + self.c.clone(), self.sigma.clone(), self.c.clone());
        prob_of_x[bound].div_assign_round(Float::with_val(32, 2), Round::Nearest);

        let mut i = -1;
        for j in 0..self.precision as usize {
            let pow: Float = Float::with_val(32, 2).pow(i); // 2^{i}
            i -= 1;
            for x in (1..=bound).rev() {
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
        
        aux_begin.reserve_exact(p_num_rows);

        // computing in which position the non-zero values in P start and end
        for i in 0..p_num_rows {
            aux_begin[i] = p_num_cols as i32 - 1;

            for j in 0..p_num_cols {
                if self.p[i][j] == 1 {
                    aux_begin[i] = j as i32;
                    break;
                }
            }
        }

        self.begin = aux_begin;
    }

    fn probability(&self, x: Float, sigma: Float, c: Float) -> Float {
        let pi = Float::with_val(32, Constant::Pi);
        let mut s: Float = sigma.clone() * (Float::with_val(32, 2) * pi).sqrt();
        let mut over_s: Float = 1 / s;

        if x == 0 {
            return over_s;
        }

        // over_s * exp(-(power((x - c) / sigma, 2)) / 2.0)
        let tmp = (x - c) / sigma;
        let tmp2 = tmp.clone() * tmp;
        let tmp3 = -(tmp2 / Float::with_val(32, 2.0));
        let tmp4 = tmp3.exp();
        over_s * tmp4
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
