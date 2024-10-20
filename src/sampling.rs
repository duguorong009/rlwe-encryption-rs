pub const ROWS: usize = 128;
pub const COLS: usize = 40;

pub struct Sampling {
    p: Vec<Vec<i32>>,
    begin: Vec<i32>,
    precision: i32,
    tailcut: f32,
    sigma: RR,
    c: RR,
}

impl Sampling {
    pub fn new(precision: i32, tailcut: f32, sigma: RR, center: RR) -> Self {
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

    fn probability(&self, x: RR, sigma: RR, c: RR) -> RR {
        todo!()
    }

    fn binary_expansion(
        &self,
        aux_p: Vec<Vec<i32>>,
        probability: RR,
        precision: u64,
        index: usize,
    ) {
        todo!()
    }

    // bit = 0 then return a
    fn select(a: i32, b: i32, bit: u32) -> i32 {
        let mask = -(bit as i32);
        (mask & (a ^ b)) ^ a
    }
}
