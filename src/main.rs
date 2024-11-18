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
