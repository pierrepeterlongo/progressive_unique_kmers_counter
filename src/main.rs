use clap::Parser; 
use futures::{SinkExt, StreamExt};
use std::path::PathBuf;
use std::sync::Arc;
use tokio::sync::{mpsc, Mutex};
use warp::ws::{Message, WebSocket};
use warp::Filter;
use fxread;
use std::error::Error;


/// Fast hash map
use rustc_hash::FxHashMap;

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Length of k-mers
    #[arg(short, long)]
    k: usize,

    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,

    /// Stop if acceleration is below this threshold
    #[arg(long, default_value_t = 10.0)]
    stop_acceleration: f32,
}


// enum RecordReader<R: Read> {
//     Fasta(fasta::Records<BufReader<R>>),
//     Fastq(fastq::Records<BufReader<R>>),
// }

// impl<R: Read> RecordReader<R> {
//     fn next_record(&mut self) -> Option<Result<Vec<u8>, Box<dyn std::error::Error>>> {
//         match self {
//             RecordReader::Fasta(reader) => reader.next().map(|r| {
//                 r.map(|rec| rec.seq().to_vec())
//                     .map_err(|e| e.into())
//             }),
//             RecordReader::Fastq(reader) => reader.next().map(|r| {
//                 r.map(|rec| rec.seq().to_vec())
//                     .map_err(|e| e.into())
//             }),
//         }
//     }
// }

fn encode_kmer_u64(kmer: &[u8]) -> Option<u64> {
    if kmer.len() > 32 {
        return None;
    }
    let mut value: u64 = 0;
    for &base in kmer {
        value <<= 2;
        match base {
            b'A' | b'a' => value |= 0,
            b'C' | b'c' => value |= 1,
            b'G' | b'g' => value |= 2,
            b'T' | b't' => value |= 3,
            _ => return None, // non-ACGT
        }
    }
    Some(value)
}

fn reverse_complement_u64(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 0b11;
        let comp = base ^ 0b11; // complement (A<->T, C<->G)
        rc |= comp << (2 * (k - i - 1));
    }
    rc
}

fn canonical_kmer_u64(kmer: u64, k: usize) -> u64 {
    let rc = reverse_complement_u64(kmer, k);
    std::cmp::min(kmer, rc)
}

/// Fast reverse complement for &[u8]
fn _reverse_complement(kmer: &[u8]) -> Vec<u8> {
    kmer.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => c,
        })
        .collect()
}

/// Return the canonical k-mer (lexicographically smallest between kmer and its reverse complement)
fn _canonical_kmer(kmer: &[u8]) -> Vec<u8> {
    let rc = _reverse_complement(kmer);
    if kmer <= rc.as_slice() {
        kmer.to_vec()
    } else {
        rc
    }
}

/// WebSocket handling
async fn handle_connection(ws: WebSocket, rx: Arc<Mutex<mpsc::Receiver<(u32, u32)>>>) {
    let (mut ws_tx, _) = ws.split();
    let mut rx = rx.lock().await;
    while let Some((reads, kmers)) = rx.recv().await {
        let message = format!("{} {}", reads, kmers);
        if ws_tx.send(Message::text(message)).await.is_err() {
            break;
        }
    }
}
//
// fn open_reader(path: &PathBuf) -> Result<RecordReaderWrapper<Box<dyn ReadSeek>>, Box<dyn Error>> {
//     let file = File::open(path)?;
//     let reader: Box<dyn ReadSeek> = if path.extension().map(|e| e == "gz").unwrap_or(false) {
//         Box::new(MultiGzDecoder::new(file))
//     } else {
//         Box::new(file)
//     };
//
//     RecordReaderWrapper::new(reader)
// }


// fn open_reader(path: &PathBuf) -> Result<RecordReader<impl Read>, Box<dyn std::error::Error>> {
//     let file = File::open(path)?;
//     let reader: Box<dyn Read> = if path.extension().map(|e| e == "gz").unwrap_or(false) {
//         Box::new(MultiGzDecoder::new(file))
//     } else {
//         Box::new(file)
//     };

//     let mut buffered = BufReader::new(reader);

//     // Peek at the first byte
//     let first_byte = {
//         let buf = buffered.fill_buf()?;
//         if buf.is_empty() {
//             return Err("Input file is empty".into());
//         }
//         buf[0]
//     };

//     // Decide format by first byte
//     if first_byte == b'>' {
//         Ok(RecordReader::Fasta(fasta::Reader::new(buffered).records()))
//     } else if first_byte == b'@' {
//         Ok(RecordReader::Fastq(fastq::Reader::new(buffered).records()))
//     } else {
//         Err(format!("Unknown file format: expected '>' or '@', got '{}'", first_byte as char).into())
//     }
// }

// function read_idx_records_and_get_position to read "idx" reads from a fastx file and return the position in the file
use fxread::initialize_reader;
fn read_idx_records_and_get_position(path: &PathBuf, idx: usize) -> Result<u64, Box<dyn Error>> {
    let mut reader = initialize_reader(path).unwrap();
    let mut position = 0;
    let mut current_idx = 0;
    
    while let Some(record) = reader.next() {
        // if fasta sum the id and seq lengths
        // if fastq sum the id, seq, and quality lengths
        if current_idx >= idx {
            break;
        }
        if record.is_fasta() {
            position += record.id().len() + record.seq().len() + 2; // +2 for newline characters
        } else if record.is_fastq() {
            position += record.id().len() + record.seq().len() + record.qual().unwrap().len() + 4; // +4 for newline characters and '+' line
        } else {
            return Err("Unknown file format".into());
        }
        current_idx += 1;
    }

    Ok(position as u64)
}



#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let k = args.k;
    let stop_acceleration = args.stop_acceleration;

    // let mut unique_kmers: FxHashMap<Vec<u8>, bool> = FxHashMap::default();
    let mut unique_kmers: FxHashMap<u64, bool> = FxHashMap::default();
    let mut unique_solid_kmers = 0;

    let (tx, rx) = mpsc::channel(100);
    let rx = Arc::new(Mutex::new(rx));

    // WebSocket route
    let ws_route = warp::path("ws")
        .and(warp::ws())
        .map(move |ws: warp::ws::Ws| {
            let rx = rx.clone();
            ws.on_upgrade(move |socket| handle_connection(socket, rx))
        });

    tokio::spawn(async move {
        warp::serve(ws_route)
            .run(([127, 0, 0, 1], 3030))
            .await;
    });

    let mut reader = fxread::initialize_reader(&args.input)?;
    let mut idx = 0;

    let mut prev_kmers = 0u32;
    let mut avg_growth = 0.0;
    let mut growth_history: Vec<i32> = Vec::new();
    let mut accel_history: Vec<i32> = Vec::new();

    while let Some(seq_result) = reader.next_record()? {
        let sequence = seq_result.as_str();
        // check the length of the sequence
        if sequence.len() < k {
            idx += 1;
            continue;
        }

        for i in 0..=(sequence.len() - k) {
            let kmer = &sequence[i..i + k];
            if let Some(encoded) = encode_kmer_u64(kmer.as_bytes()) {
                    let canonical = canonical_kmer_u64(encoded, k);
                    match unique_kmers.get_mut(&canonical) {
                        Some(seen) => {
                            if !*seen {
                                *seen = true;
                                unique_solid_kmers += 1;
                            }
                        }
                        None => {
                            unique_kmers.insert(canonical, false);
                        }
                    }
            }
        }

        if idx % 10000 == 0 {

            let reads = idx as u32;
            let kmers = unique_solid_kmers;
            let growth = kmers as i32 - prev_kmers as i32;

            growth_history.push(growth);
            if growth_history.len() > 10 {
                growth_history.remove(0);
            }

            // Compute acceleration only if we have at least 2 growth values
            if growth_history.len() >= 2 {
                let acceleration = growth_history[growth_history.len() - 1]
                    - growth_history[growth_history.len() - 2];
                accel_history.push(acceleration);
                if accel_history.len() > 10 {
                    accel_history.remove(0);
                }
            }

            // Compute averages
            avg_growth = growth_history.iter().sum::<i32>() as f32 / growth_history.len() as f32;
            let avg_accel = if !accel_history.is_empty() {
                accel_history.iter().sum::<i32>() as f32 / accel_history.len() as f32
            } else {
                0.0
            };

            println!(
                "Processed {} reads, unique k-mers: {}, Δ_avg: {:.1}, Δ²_avg: {:.1}",
                reads, kmers, avg_growth, avg_accel
            );

            // WebSocket message can include acceleration too if desired
            tx.send((reads, kmers)).await?;

            // Auto-stop condition
            if reads > 50000 && avg_accel.abs() < stop_acceleration {
                println!(
                    "Stopping early: acceleration average {:.1} < {} after {} reads.",
                    avg_accel, stop_acceleration, reads
                );
                break;
            }

            prev_kmers = kmers;
        }

        idx += 1;
    }
    println!(
        "Final unique k-mers: {}, after processing {} reads.",
        unique_solid_kmers, idx
    );


	// get the reached position in the file
    let reached_position = read_idx_records_and_get_position(&args.input, idx)?;
    // get the file size
    let file_size = args.input.metadata()?.len();
    // TODO: does not work for gzipped files
    // special case for gzipped files
    if args.input.extension().map(|e| e == "gz").unwrap_or(false) {
        println!("Input file is gzipped, file size may not be accurate.");
    }
    println!(
        "Reached position in file: {} (for {} reads), file size: {}",
        reached_position, idx, file_size
    );
    // estimate the total number of reads in the file
    let total_reads = (file_size as f64 / reached_position as f64 * idx as f64) as u64;
    // remaining reads
    let remaining_reads = total_reads - idx as u64;
    println!(
        "Estimated total reads: {}, remaining reads: {}",
        total_reads, remaining_reads
    );

    // we have the total number of remaining reads, given that we know the growth we can estimate the total number of remaining kmers
    // we use the last avg growth value to estimate the remaining kmers
    
    let estimated_remaining_kmers =  remaining_reads as f64 * avg_growth as f64 / 10000.0;
    println!(
        "Estimated remaining unique k-mers: {:.0} (based on last growth value of {:.1})",
        estimated_remaining_kmers, avg_growth
    );
    println!("Estimated total unique k-mers: {:.0}", unique_solid_kmers as f64 + estimated_remaining_kmers);


    Ok(())
}
