use clap::Parser; 
use futures::{SinkExt, StreamExt};
use core::num;
use std::path::PathBuf;
use std::sync::Arc;
use tokio::sync::{mpsc, Mutex};
use warp::ws::{Message, WebSocket};
use warp::Filter;
use needletail::{parse_fastx_file, Sequence};
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

    /// Give the number of reads
    #[arg(long, default_value_t = 0)]
    nb_reads: u64,

    /// Plot the intermediate results in plot.html
    #[arg(long, default_value_t = false)]
    plot: bool,
}


fn _encode_kmer_u64(kmer: &[u8]) -> Option<u64> {
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

fn _reverse_complement_u64(kmer: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for i in 0..k {
        let base = (kmer >> (2 * i)) & 0b11;
        let comp = base ^ 0b11; // complement (A<->T, C<->G)
        rc |= comp << (2 * (k - i - 1));
    }
    rc
}

fn _canonical_kmer_u64(kmer: u64, k: usize) -> u64 {
    let rc = _reverse_complement_u64(kmer, k);
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
// function read_idx_records_and_get_position to read "idx" reads from a fastx file and return the position in the file

fn read_idx_records_and_get_position(path: &PathBuf, idx: usize) -> Result<u64, Box<dyn Error>> {
    let mut reader = parse_fastx_file(&path).expect("valid path/file");

    // let mut reader = initialize_reader(path).unwrap();
    let mut position = 0;
    let mut current_idx = 0;
    // Iterate through the records until we reach the idx-th record
    // or until we reach the end of the file
    // We assume that the file is not empty and has at least idx records
    if idx == 0 {
        return Ok(0);
    }
    while let Some(record) = reader.next() {
        if current_idx >= idx {
            break;
        }
        // get the sequence record
        let seq_record = record.expect("valid record");
        // Since needletail does not provide is_fasta/is_fastq, we can check for quality
        if seq_record.qual().is_none() {
            // FASTA: no quality
            position += seq_record.id().len() + seq_record.seq().len() + 2; // +2 for newline characters
        } else {
            // FASTQ: has quality
            position += seq_record.id().len() + seq_record.seq().len() + seq_record.qual().unwrap().len() + 4; // +4 for newline characters and '+' line
        }
        current_idx += 1;
        // reader.next(); // Move to the next record (already done by while let Some(record) = reader.next())
    }

    Ok(position as u64)
}



#[tokio::main]
async fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let k = args.k;
    let stop_acceleration = args.stop_acceleration;

    let mut unique_kmers: FxHashMap<Vec<u8>, bool> = FxHashMap::default();
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

    
    let mut reader = parse_fastx_file(&args.input).expect("valid path/file");
    // let mut reader = fxread::initialize_reader(&args.input)?;
    let mut idx = 0;

    let mut prev_kmers = 0u32;
    let mut avg_growth = 0.0;
    let mut growth_history: Vec<i32> = Vec::new();
    let mut accel_history: Vec<i32> = Vec::new();
    let min_number_low_acceleration = 3; // Minimum number of low acceleration values to consider stopping
    // Initialize the growth history with a single value of 0
    let mut number_low_acceleration = 0;
    while let Some(record) = reader.next() {
        
        idx += 1;
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        let rc = norm_seq.reverse_complement();
        for (_, kmer, _) in norm_seq.canonical_kmers(k as u8, &rc) {
            // for debug purpose, print the k-mer in ascii, transforming the &[u8] to a string
            // println!("Processing k-mer: {:?}", String::from_utf8_lossy(kmer));

                match unique_kmers.get_mut(kmer) {
                        Some(seen) => {
                            if !*seen {
                                *seen = true;
                                unique_solid_kmers += 1;
                                // println!(
                                //     "Found new unique k-mer: {} (total: {})",
                                //     String::from_utf8_lossy(kmer),
                                //     unique_solid_kmers
                                // );
                            }
                        }
                        None => {
                            unique_kmers.insert(kmer.to_vec(), false);
                            // println!(
                            //     "Inserted new k-mer: {}",
                            //     String::from_utf8_lossy(kmer)
                            // );
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
                "Nb reads: {} Nb unique k-mers: {}, growth: {:.1}, acceleration: {:.1}",
                reads, kmers, avg_growth, avg_accel
            );

            
            // WebSocket message 
            if args.plot {
               tx.send((reads, kmers)).await?;
            }
            // Auto-stop condition
            if reads > 50000 && avg_accel.abs() < stop_acceleration {
                number_low_acceleration += 1;

                println!(
                    "Low acceleration average {:.1}  number {}/{}.",
                    avg_accel, number_low_acceleration, min_number_low_acceleration
                );
                if number_low_acceleration >= min_number_low_acceleration {
                    println!(
                        "Stopping early: acceleration average {:.1} < {} after {} reads.",
                        avg_accel, stop_acceleration, reads
                    );
                    break;
                }
            }
            prev_kmers = kmers;
        }

    }
    println!(
        "Final unique k-mers: {}, after processing {} reads.",
        unique_solid_kmers, idx
    );


    // special case for gzipped files
    if args.input.extension().map(|e| e == "gz").unwrap_or(false) && args.nb_reads == 0 {
        println!("Input file is gzipped, cannot estimated total number of kmers.");
        return Ok(());
    }


	// get the reached position in the file
    let reached_position = read_idx_records_and_get_position(&args.input, idx)?;
    // get the file size
    let file_size = args.input.metadata()?.len();

    println!(
        "Reached position in file: {} (for {} reads), file size: {}",
        reached_position, idx, file_size
    );
    let mut total_reads = args.nb_reads;

    // if the number of reads is given, we can estimate the total number of kmers
    if total_reads > 0 {
        println!("Total number of reads given: {}", total_reads);
        if idx > total_reads as usize {
            print!("Warning: idx ({}) is greater than total reads ({}). ", idx, total_reads);
            println!("I use {} reads as total", total_reads);
            idx = total_reads as usize;
        }
        // we can estimate the total number of kmers based on the number of reads and the reached position
    } else {
        total_reads = (file_size as f64 / reached_position as f64 * idx as f64) as u64;
        
        println!(
            "Estimated total reads: {}",
            total_reads
        );
    }

    // we have the total number of remaining reads, given that we know the growth we can estimate the total number of remaining kmers
    // we use the last avg growth value to estimate the remaining kmers
    
    // remaining reads
    let remaining_reads = total_reads - idx as u64;
    println!("Remaining reads: {}", remaining_reads);

    let estimated_remaining_kmers =  remaining_reads as f64 * avg_growth as f64 / 10000.0;
    println!(
        "Estimated remaining unique k-mers: {:.0} (based on last growth value of {:.1})",
        estimated_remaining_kmers, avg_growth
    );
    println!("Estimated total unique k-mers: {:.0}", unique_solid_kmers as f64 + estimated_remaining_kmers);


    Ok(())
}
