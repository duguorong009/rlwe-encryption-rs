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

    //  Method for computing the binary expansion of a given probability in [0, 1] 
    fn binary_expansion(
        &self,
        aux_p: &mut Vec<Vec<i32>>,
        probability: RR,
        precision: u64,
        index: usize,
    ) {
        RR mut pow = todo!(); 
        let i = -1;
        let j = 0;

        while probability > 0 && j < precision {
            pow = power2_RR(i); // 2 ^ {i}
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
