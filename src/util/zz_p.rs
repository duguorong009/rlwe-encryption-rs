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

#[derive(Debug, Clone)]
struct ZZ_p {
    value: ZZ,
    context: ZZ_pContext,
}

impl ZZ_p {
    fn new(value: ZZ, context: ZZ_pContext) -> Self {
        Self { value, context }
    }
}
