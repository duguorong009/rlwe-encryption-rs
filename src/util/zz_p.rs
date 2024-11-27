use rug::Integer as ZZ;

#[derive(Debug, Clone)]
struct ZZ_pContext {
    p: ZZ,
}

impl ZZ_pContext {
    fn new(p: ZZ) -> Self {
        assert!(p > 1, "ZZ_pContext: p must be > 1");
        Self { p }
    }

    pub fn modulus(&self) -> &ZZ {
        &self.p
    }
}
