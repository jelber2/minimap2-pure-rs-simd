use std::io::{self, BufRead, Read};
use flate2::read::GzDecoder;
use noodles::fasta;
use noodles::fastq;

/// A single sequence record from FASTA/FASTQ.
#[derive(Clone, Debug)]
pub struct BseqRecord {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub comment: String,
    pub l_seq: usize,
}

/// FASTA/FASTQ file reader supporting plain text and gzip.
pub enum BseqReader {
    Fasta(fasta::io::Reader<Box<dyn Read>>),
    Fastq(fastq::io::Reader<Box<dyn Read>>),
}

pub struct BseqFile {
    reader: BseqReader,
}

impl BseqFile {
    /// Open a FASTA/FASTQ file (gzip auto-detected via magic bytes).
    pub fn open(path: &str) -> io::Result<Self> {
        let file: Box<dyn Read> = if path == "-" {
            Box::new(io::stdin())
        } else {
            let f = std::fs::File::open(path)?;
            let mut peek = [0u8; 2];
            let mut f = io::BufReader::new(f);
            let n = f.read(&mut peek)?;
            if n >= 2 && peek[0] == 0x1f && peek[1] == 0x8b {
                let chain = io::Cursor::new(peek[..n].to_vec()).chain(f);
                Box::new(GzDecoder::new(chain))
            } else {
                let chain = io::Cursor::new(peek[..n].to_vec()).chain(f);
                Box::new(chain)
            }
        };

        // Peek first byte to detect format
        let mut peek_buf = [0u8; 1];
        let file_ref = &mut io::BufReader::new(file);
        let n = file_ref.read(&mut peek_buf)?;
        
        let is_fastq = if n > 0 {
            peek_buf[0] == b'@'
        } else {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Empty file"));
        };

        // Re-open file for reading
        let file: Box<dyn Read> = if path == "-" {
            Box::new(io::stdin())
        } else {
            let f = std::fs::File::open(path)?;
            let mut peek = [0u8; 2];
            let mut f = io::BufReader::new(f);
            let n = f.read(&mut peek)?;
            if n >= 2 && peek[0] == 0x1f && peek[1] == 0x8b {
                let chain = io::Cursor::new(peek[..n].to_vec()).chain(f);
                Box::new(GzDecoder::new(chain))
            } else {
                let chain = io::Cursor::new(peek[..n].to_vec()).chain(f);
                Box::new(chain)
            }
        };

        let reader = if is_fastq {
            BseqReader::Fastq(fastq::io::Reader::new(file))
        } else {
            BseqReader::Fasta(fasta::io::Reader::new(file))
        };

        Ok(Self { reader })
    }

    /// Read one FASTA/FASTQ record. Returns None at EOF.
    pub fn read_record(&mut self) -> io::Result<Option<BseqRecord>> {
        match &mut self.reader {
            BseqReader::Fasta(reader) => {
                let mut record = fasta::Record::default();
                match reader.read_record(&mut record) {
                    Ok(0) => Ok(None),
                    Ok(_) => {
                        let name = String::from_utf8_lossy(record.name()).to_string();
                        let comment = String::from_utf8_lossy(record.description().unwrap_or(b"")).to_string();
                        let mut seq = record.sequence().to_vec();
                        u_to_t(&mut seq);
                        let l_seq = seq.len();
                        Ok(Some(BseqRecord {
                            name,
                            seq,
                            qual: Vec::new(),
                            comment,
                            l_seq,
                        }))
                    }
                    Err(e) => Err(io::Error::new(io::ErrorKind::Other, e.to_string())),
                }
            }
            BseqReader::Fastq(reader) => {
                let mut record = fastq::Record::default();
                match reader.read_record(&mut record) {
                    Ok(0) => Ok(None),
                    Ok(_) => {
                        let name = String::from_utf8_lossy(record.name()).to_string();
                        let comment = String::from_utf8_lossy(record.description().unwrap_or(b"")).to_string();
                        let mut seq = record.sequence().to_vec();
                        u_to_t(&mut seq);
                        let qual = record.quality_scores().to_vec();
                        let l_seq = seq.len();
                        Ok(Some(BseqRecord {
                            name,
                            seq,
                            qual,
                            comment,
                            l_seq,
                        }))
                    }
                    Err(e) => Err(io::Error::new(io::ErrorKind::Other, e.to_string())),
                }
            }
        }
    }

    /// Read a batch of sequences up to chunk_size total bases.
    pub fn read_batch(&mut self, chunk_size: i64, with_qual: bool) -> io::Result<Vec<BseqRecord>> {
        let mut records = Vec::new();
        let mut total_len: i64 = 0;
        loop {
            match self.read_record()? {
                None => break,
                Some(mut rec) => {
                    if !with_qual {
                        rec.qual.clear();
                    }
                    total_len += rec.l_seq as i64;
                    records.push(rec);
                    if total_len >= chunk_size {
                        break;
                    }
                }
            }
        }
        Ok(records)
    }

    pub fn is_eof(&self) -> bool {
        false // Noodles readers don't expose EOF state; rely on read_record returning None
    }
}

/// Convert U/u to T/t in sequence.
fn u_to_t(seq: &mut [u8]) {
    for b in seq.iter_mut() {
        if *b == b'u' || *b == b'U' {
            *b -= 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp_file(content: &[u8]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(content).unwrap();
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_read_fasta() {
        let content = b">seq1 comment1\nACGTACGT\nTGCA\n>seq2\nAAAA\n";
        let f = write_temp_file(content);
        let mut reader = BseqFile::open(f.path().to_str().unwrap()).unwrap();

        let rec1 = reader.read_record().unwrap().unwrap();
        assert_eq!(rec1.name, "seq1");
        assert_eq!(rec1.comment, "comment1");
        assert_eq!(rec1.seq, b"ACGTACGTTGCA");
        assert_eq!(rec1.l_seq, 12);

        let rec2 = reader.read_record().unwrap().unwrap();
        assert_eq!(rec2.name, "seq2");
        assert_eq!(rec2.seq, b"AAAA");

        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn test_read_fastq() {
        let content = b"@read1 comment\nACGT\n+\nIIII\n@read2\nTGCA\n+\nHHHH\n";
        let f = write_temp_file(content);
        let mut reader = BseqFile::open(f.path().to_str().unwrap()).unwrap();

        let rec1 = reader.read_record().unwrap().unwrap();
        assert_eq!(rec1.name, "read1");
        assert_eq!(rec1.seq, b"ACGT");
        assert_eq!(rec1.qual, b"IIII");

        let rec2 = reader.read_record().unwrap().unwrap();
        assert_eq!(rec2.name, "read2");
        assert_eq!(rec2.seq, b"TGCA");
        assert_eq!(rec2.qual, b"HHHH");

        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn test_u_to_t_conversion() {
        let content = b">seq1\nACGUacgu\n";
        let f = write_temp_file(content);
        let mut reader = BseqFile::open(f.path().to_str().unwrap()).unwrap();
        let rec = reader.read_record().unwrap().unwrap();
        assert_eq!(rec.seq, b"ACGTacgt");
    }

    #[test]
    fn test_read_batch() {
        let content = b">s1\nACGT\n>s2\nTGCA\n>s3\nAAAA\n";
        let f = write_temp_file(content);
        let mut reader = BseqFile::open(f.path().to_str().unwrap()).unwrap();

        let batch = reader.read_batch(5, false).unwrap();
        assert_eq!(batch.len(), 2); // 4 + 4 >= 5, stops after 2
        assert_eq!(batch[0].name, "s1");
        assert_eq!(batch[1].name, "s2");
    }
}
