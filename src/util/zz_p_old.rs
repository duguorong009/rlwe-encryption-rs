use std::{cell::{Cell, OnceCell, RefCell}, rc::Rc, sync::{Arc, LazyLock, Mutex, OnceLock}};
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
    size: usize, // p.size()
    extended_modulus_size: usize,

    fft_info: OnceLock<ZZ_pFFTInfoT>, // Lazy<ZZ_pFFTInfoT> FFTInfo; # in C++
}

impl ZZ_pInfoT {
    pub fn new(p: ZZ) -> Self {
        if p <= 1 {
            panic!("ZZ_pContext: p must be > 1");
        }

        let p = p;
        let size = p.significant_digits::<u64>();
        let extended_modulus_size = 2 * size + (64 + 64 - 1) / 64;

        Self {
            p,
            size,
            extended_modulus_size,
            fft_info: OnceLock::new(),
        }
    }
}

struct ZZ_pTmpSpaceT {
    crt_tmp_vec: ZZ_TmpVecAdapter,
    rem_tmp_vec: ZZ_TmpVecAdapter,
}

struct ZZ_pContext {
    info: Rc<ZZ_pInfoT>, // SmartPtr<ZZ_pInfoT> ptr; # in C++
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

pub static ZZ_P_INFO_STG: LazyLock<Mutex<Option<Arc<ZZ_pInfoT>>>> = LazyLock::new(|| Mutex::new(None));
pub static ZZ_P_TMP_SPACE_STG: LazyLock<Mutex<Option<Arc<ZZ_pTmpSpaceT>>>> = LazyLock::new(|| Mutex::new(None));

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
