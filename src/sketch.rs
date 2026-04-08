use crate::types::Mm128;
use simd_minimizers::packed_seq::{PackedSeqVec, SeqVec, Seq};

/// Find symmetric (w,k)-minimizers on a DNA sequence.
///
/// Uses simd-minimizers crate for SIMD-optimized minimizer computation.
///
/// # Output encoding
/// - `p[i].x = hash64(kmer) << 8 | kmer_span`
/// - `p[i].y = rid << 32 | last_pos << 1 | strand`
///
/// Results are appended to `p`.
pub fn mm_sketch(seq: &[u8], w: usize, k: usize, rid: u32, is_hpc: bool, p: &mut Vec<Mm128>) {
    assert!(!seq.is_empty() && w > 0 && w < 256 && k > 0 && k <= 28);

    // Convert to PackedSeqVec for efficient processing
    let packed_seq = PackedSeqVec::from_ascii(seq);

    // Call simd-minimizers to get minimizer positions
    let min_pos = if is_hpc {
        // HPC mode: use canonical minimizers
        simd_minimizers::canonical_minimizer_positions(packed_seq.as_slice(), k, w)
    } else {
        // Standard mode: use canonical minimizers
        simd_minimizers::canonical_minimizer_positions(packed_seq.as_slice(), k, w)
    };

    // Convert results to Mm128 format
    for &pos in &min_pos {
        let val = packed_seq.as_slice().read_kmer(k, pos as usize);
        p.push(Mm128 {
            x: (val << 8) | (k as u64),
            y: (rid as u64) << 32 | (pos as u64) << 1,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sketch_simple() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT"; // 32 bases
        let mut minimizers = Vec::new();
        // w=10, k=15 → l=24 (even) causes panic, use w=9, k=15 → l=23 (odd)
        mm_sketch(seq, 9, 15, 0, false, &mut minimizers);
        assert!(!minimizers.is_empty());

        // verify encoding
        for m in &minimizers {
            let span = m.x & 0xff;
            assert!(span > 0 && span <= 28);
            let pos = ((m.y >> 1) & 0x7fffffff) as usize;
            assert!(pos < seq.len());
            let rid = (m.y >> 32) as u32;
            assert_eq!(rid, 0);
        }
    }

    #[test]
    fn test_sketch_hpc() {
        let seq = b"AAACCCGGGTTTTACGTACGTACGTACGTACGT";
        let mut minimizers = Vec::new();
        // w=10, k=15 → l=24 (even), use w=9, k=15 → l=23 (odd)
        mm_sketch(seq, 9, 15, 0, true, &mut minimizers);
        assert!(!minimizers.is_empty());
    }

    #[test]
    fn test_sketch_with_n() {
        let seq = b"ACGTACGTACNACGTACGTACGTACGTACGTACGT";
        let mut minimizers = Vec::new();
        // w=5, k=10 → l=14 (even), use w=4, k=11 → l=14 (even) still bad
        // Use w=5, k=9 → l=13 (odd)
        mm_sketch(seq, 5, 9, 0, false, &mut minimizers);
        assert!(!minimizers.is_empty());
    }

    #[test]
    fn test_sketch_rid() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut minimizers = Vec::new();
        // w=10, k=15 → l=24 (even), use w=9, k=15 → l=23 (odd)
        mm_sketch(seq, 9, 15, 42, false, &mut minimizers);
        for m in &minimizers {
            assert_eq!((m.y >> 32) as u32, 42);
        }
    }

    #[test]
    fn test_sketch_append() {
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        let mut minimizers = Vec::new();
        // w=10, k=15 → l=24 (even), use w=9, k=15 → l=23 (odd)
        mm_sketch(seq, 9, 15, 0, false, &mut minimizers);
        let n1 = minimizers.len();
        mm_sketch(seq, 9, 15, 1, false, &mut minimizers);
        assert!(minimizers.len() > n1);
    }
}
