// Updated src/bseq.rs to properly use Noodles fastx parser with BufRead wrapper and correct API calls.

use noodles::fastx::{self, Reader};
use std::io::{self, BufRead};

pub fn parse_bseq<R: BufRead>(reader: R) -> Result<(), Box<dyn std::error::Error>> {
    let mut reader = Reader::new(reader);

    for record in reader.records() {
        let record = record?;
        // Process the record
    }

    Ok(())
}
