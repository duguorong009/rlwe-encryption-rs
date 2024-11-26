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

struct ZZ_pInfoT {
    p: ZZ, // the modulus
    size: u64, // p.size()
    extended_modulus_size: u64,

    fft_info: ZZ_pFFTInfoT, // Lazy<ZZ_pFFTInfoT> FFTInfo; # in C++
}

struct ZZ_pTmpSpaceT {
    crt_tmp_vec: ZZ_TmpVecAdapter,
    rem_tmp_vec: ZZ_TmpVecAdapter,
}

struct ZZ_pContext {
    info: ZZ_pInfoT, // SmartPtr<ZZ_pInfoT> ptr; # in C++
}

impl ZZ_pContext {
    pub fn new(p: ZZ) {
        todo!();
    }

    pub fn save() {
        todo!();
    }

    pub fn restore() {
        todo!();
    }
}

struct ZZ_pBak {
    c: ZZ_pContext,
    must_restore: bool,
}

impl ZZ_pBak {
    pub fn new() -> Self {
        todo!();
    }

    pub fn save() {
        todo!();
    }

    pub fn restore() {
        todo!();
    }
}

struct ZZ_pPush {
    bak: ZZ_pBak,
}

impl ZZ_pPush {
    pub fn new() -> Self {
        todo!();
    }
}

thread_local! {
    // info for current modulus, initially null
    // plain pointer for faster TLS access
    pub static ZZ_pInfo: RefCell<Option<ZZ_pInfoT>> = RefCell::new(None);

    // space for temps associated with current modulus, 
    // plain pointer for faster TLS access
    pub static ZZ_pTmpSpace: RefCell<Option<ZZ_pTmpSpaceT>> = RefCell::new(None);

    // flag indicating if current modulus is fully installed
    pub static ZZ_pInstalled: Cell<bool> = Cell::new(false);
}
