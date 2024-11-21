pub mod zzx;
pub mod zz_p;

use rand::Rng;

/// Equivalent to `NTL::RandomBits_ulong`
pub fn randombits_u64(bits: u32) -> u64 {
    assert!(bits > 0 && bits <= 64, "Bits must be between 1 and 64");

    // Generate a random u64
    let mut rng = rand::thread_rng();
    let random_value: u64 = rng.gen::<u64>();

    // Mask the result to ensure it fits within the specified number of bits
    let mask = if bits == 64 {
        u64::MAX
    } else {
        (1 << bits) - 1
    };
    random_value & mask
}

/// Equivalent to `NTL::RandomBits_long`
pub fn randombits_i64(bits: u32) -> i64 {
    assert!(bits > 0 && bits <= 64, "Bits must be between 1 and 64");

    // Generate a random u64
    let mut rng = rand::thread_rng();
    let random_value: i64 = rng.gen::<i64>();

    // Mask the result to ensure it fits within the specified number of bits
    let mask = if bits == 64 {
        i64::MAX
    } else {
        (1 << bits) - 1
    };
    random_value & mask
}
