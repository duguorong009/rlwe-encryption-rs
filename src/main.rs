use encryption_scheme::EncryptionScheme;
use rug::Float as RR;
use util::{randombits_i64, randombits_u64, zzx::*};

mod encryption_scheme;
mod sampling;
mod util;

// // RLWE
// const P: usize = 1024;
// const Q: usize = 11289;
// const SIGMA: f32 = 3.19;

// ALTERNATE
const P: usize = 14; // poly degree
const Q: usize = 179424673; // 15485863, 8380417
const SIGMA: f32 = 2.0;

// // NTRU: NTRU Prime parameters
// const P: usize = 761;
// const Q: usize = 4591;
// const SIGMA: f32 = 2.0;

const BENCH_LOOPS: usize = 10;

fn main() {
    println!(
        "----------------------------\nRing-LWE encryption scheme\n----------------------------\n"
    );

    let precision = 256;
    let center = RR::with_val(precision, 0);
    let tail_cut = 13.2;

    let mut total_errors = 0;

    let sigma = RR::with_val(precision, SIGMA);
    let es = EncryptionScheme::new(P as i32, Q as i32, precision, tail_cut, sigma, center);

    /* key generation */
    let a = random_poly();
    let mut r2 = ZZX::new();
    let mut p1 = ZZX::new();
    es.key_generation(&a, &mut r2, &mut p1);

    for _ in 0..BENCH_LOOPS {
        let m = random_message();

        // println!("Message being encrypted: {:?}\n", m);

        // encryption
        let mut mprime = ZZX::new();
        es.encode(&mut mprime, &m);

        let mut c1 = ZZX::new();
        let mut c2 = ZZX::new();
        es.encryption(&mut c1, &mut c2, &a, &p1, &mprime);

        // println!("Encrypted message: {:?}\n", c1);

        // decryption
        let mut moriginal = ZZX::new();
        es.decryption(&mut moriginal, &c1, &c2, &r2);

        let mut mdecoded = vec![0; P];
        es.decode(&mut mdecoded, &moriginal);

        // println!("Decrypted message: {:?}\n", mdecoded);

        let mut counter = 0;
        for i in 0..P {
            if mdecoded[i] != m[i] {
                // println!("{}: {} != {}", i, mdecoded[i], m[i]);
                counter += 1;
            }
        }
        // println!("Number of incorrect decodings: {}\n", counter);
        if counter > 0 {
            total_errors += 1;
        }

        if total_errors == 0 {
            println!("OK\n");
        } else {
            println!("{} executions have failed.", total_errors);
        }

        if total_errors == BENCH_LOOPS {
            println!("All executions have failed.");
        }
    }
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
