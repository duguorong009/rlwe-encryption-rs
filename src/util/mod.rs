pub mod zzx;

use rand::Rng;

/// Equivalent to `NTL::RandomBits_ulong`
pub fn randombits_u64(bits: u8) -> u64 {
    assert!(bits > 0 && bits <= 64, "Bits must be between 1 and 64");

    // Generate a random u64
    let mut rng = rand::rng();
    let random_value: u64 = rng.random::<u64>();

    // Mask the result to ensure it fits within the specified number of bits
    let mask = if bits == 64 {
        u64::MAX
    } else {
        (1 << bits) - 1
    };
    random_value & mask
}

/// Equivalent to `NTL::RandomBits_long`
pub fn randombits_i64(bits: u8) -> i64 {
    assert!(bits > 0 && bits <= 64, "Bits must be between 1 and 64");

    // Generate a random u64
    let mut rng = rand::rng();
    let random_value: i64 = rng.random::<i64>();

    // Mask the result to ensure it fits within the specified number of bits
    let mask = if bits == 64 {
        i64::MAX
    } else {
        (1 << bits) - 1
    };
    random_value & mask
}
