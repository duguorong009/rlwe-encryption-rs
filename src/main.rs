use util::{randombits_i64, randombits_u64, zzx::*};

mod encryption_scheme;
mod sampling;
mod util;

// // RLWE
// const P: usize = 1024;
// const Q: usize = 11289;
// const SIGMA: f32 = 3.19;

// // ALTERNATE
// const P: usize = 14; // poly degree
// const Q: usize = 179424673; // 15485863, 8380417
// const SIGMA: f32 = 2.0;

// NTRU: NTRU Prime parameters
const P: usize = 761;
const Q: usize = 4591;
const SIGMA: f32 = 2.0;

fn main() {
    println!("Hello, world!");
}

fn random_poly() -> ZZX {
    let mut a = ZZX::new();
    a.set_length(P);
    for i in 0..P {
        a[i] = randombits_u64(((Q as f64).log2() as u32) / 4).into();
    }
    a
}

fn random_message() -> Vec<i32> {
    let mut a = Vec::new();
    for _ in 0..P {
        a.push(randombits_i64(1) as i32);
    }
    a
}
