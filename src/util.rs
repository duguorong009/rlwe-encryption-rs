use rug::Integer;

/// Custom clone of NTL::ZZX 
#[derive(Clone, Debug, Default)]
struct ZZX {
    coeffs: Vec<Integer>,
}

