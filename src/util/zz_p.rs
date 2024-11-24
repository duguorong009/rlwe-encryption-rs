use rug::Integer as ZZ;

struct ZZ_pFFTInfoT {
    num_primes: i64,
    max_root: i64,
    minus_mod_p: ZZ, // -M mod p, M = product of primes
    crt_struct: ZZ_CRTStructAdapter,
    rem_struct: ZZ_RemStructAdapter,

    // following arrays are indexed 0..num_primes-1
    // q[i] = FFTPrime[i]
    prime: Vec<i64>, // prime[i] = q[i]
    prime_recip: Vec<i128>, // prime_recip[i] = 1/double(q[i])
    u: Vec<i64>, // u[i] = (M/q[i])^{-1} mod q[i]
    uqinv: Vec<u64>, 

    reduct_struct: ZZ_ReductStructAdapter,
}

