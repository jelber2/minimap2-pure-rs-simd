fn sketch(sequence: &[u8], is_hpc: bool) -> Vec<Mm128> {
    // Convert input sequence to PackedSeqVec
    let packed_seq = PackedSeqVec::from(sequence);

    // Initialize vector to hold results
    let mut results = Vec::new();

    // Call the appropriate simd_minimizers function based on is_hpc parameter
    let minimizers = if is_hpc {
        simd_minimizers::minimizer_positions(&packed_seq)
    } else {
        simd_minimizers::canonical_minimizer_positions(&packed_seq)
    };

    // Iterate over minimizers to get kmer values
    for (pos, (hash, rid, strand)) in minimizers.iter().enumerate() {
        let x = (hash << 8) | (pos as u64);
        let y = (rid << 32) | ((pos as u64) << 1) | (*strand as u64);
        results.push(Mm128 { x, y });
    }

    results
}

// Keep the test module intact
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sketch() {
        // Your test cases go here
    }
}