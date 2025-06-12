use clap::Parser; 
use std::path::PathBuf;
use needletail::{parse_fastx_file, Sequence};
use std::error::Error;
use std::usize;

// UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
// use tokio::sync::{mpsc, Mutex};
// use std::sync::Arc;
// use warp::ws::{Message, WebSocket};
// use warp::Filter;
// use futures::{SinkExt, StreamExt};


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

    /// Stop if acceleration is below or equal to this threshold
    #[arg(long, default_value_t = 10.0)]
    stop_acceleration: f64,

    /// Give the number of reads
    #[arg(long, default_value_t = 0)]
    nb_reads: u64,

    // // UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
    // /// Plot the intermediate results in plot.html
    // #[arg(long, default_value_t = false)]
    // plot: bool,

    /// Stop the evaluation when the number of low acceleration values is reached
    #[arg(long, default_value_t = 3)]
    min_number_low_acceleration: usize,
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

    // UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
// /// WebSocket handling
// async fn handle_connection(ws: WebSocket, rx: Arc<Mutex<mpsc::Receiver<(u32, u32)>>>) {
//     let (mut ws_tx, _) = ws.split();
//     let mut rx = rx.lock().await;
//     while let Some((reads, kmers)) = rx.recv().await {
//         let message = format!("{} {}", reads, kmers);
//         if ws_tx.send(Message::text(message)).await.is_err() {
//             break;
//         }
//     }
// }


// Given a set (f(x)=y) of (kmers seen, distinct solid kmers), fit a model f(x) = b*x + c*x^2 +d*x^3 of the last 50 values
// using a cubic polynomial regression
// and return the coefficients (b, c, d) as a tuple
fn fit_model(nb_kmers: &[u64], nb_dskmers: &[u64]) -> (f64, f64, f64) {
    let n = nb_kmers.len();
    // if less than 50 reads, we cannot fit a cubic model, return (0.0, 0.0, 0.0)
    if n < 50  {
        return (0.0, 0.0, 0.0); // Not enough data or mismatched lengths
    }

    // Sum the last 50 values for necessary calculations
    let nb_kmers = &nb_kmers[n - 50..];
    let nb_dskmers = &nb_dskmers[n - 50..];
    let n = nb_kmers.len() as f64;
    let sum_x = nb_kmers.iter().sum::<u64>() as f64;
    let sum_y = nb_dskmers.iter().sum::<u64>() as f64;
    let sum_x2 = nb_kmers.iter().map(|&x| (x * x) as f64).sum::<f64>();
    let sum_x3 = nb_kmers.iter().map(|&x| (x * x * x) as f64).sum::<f64>();
    let sum_x4 = nb_kmers.iter().map(|&x| (x * x * x * x) as f64).sum::<f64>();
    let sum_xy = nb_kmers.iter().zip(nb_dskmers).map(|(&x, &y)| (x * y) as f64).sum::<f64>();
    let sum_x2y = nb_kmers.iter().zip(nb_dskmers).map(|(&x, &y)| (x * x * y) as f64).sum::<f64>();

    // Solve the linear system using Cramer's rule or any other method
    let denom = n as f64 * (sum_x2 * sum_x4 - sum_x3 * sum_x3)
        - sum_x * (n as f64 * sum_x3 - sum_x2 * sum_x2);
    
    if denom == 0.0 {
        return (0.0, 0.0, 0.0); // Degenerate case
    }

    let b = (n as f64 * sum_xy - sum_x * sum_y) / denom;
    let c = (n as f64 * sum_x2y - sum_x2 * sum_y) / denom;
    let d = (sum_xy - b * sum_x - c * sum_x2) / n as f64;

    (b, c, d)
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
    
    let mut unique_kmers: FxHashMap<Vec<u8>, bool> = FxHashMap::default();
    let mut unique_solid_kmers = 0;

    // UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
    // let (tx, rx) = mpsc::channel(100);
    // let rx = Arc::new(Mutex::new(rx));

    // // WebSocket route
    // let ws_route = warp::path("ws")
    //     .and(warp::ws())
    //     .map(move |ws: warp::ws::Ws| {
    //         let rx = rx.clone();
    //         ws.on_upgrade(move |socket| handle_connection(socket, rx))
    //     });

    // tokio::spawn(async move {
    //     warp::serve(ws_route)
    //         .run(([127, 0, 0, 1], 3030))
    //         .await;
    // });

    
    let mut reader = parse_fastx_file(&args.input).expect("valid path/file");
    // let mut reader = fxread::initialize_reader(&args.input)?;
    let mut nb_kmers_seen: u32 = 0;
    let mut nb_reads_seen: usize = 0;


    // Initialize the growth history with a single value of 0
    let step = 500000; // step for printing the progress and computing acceleration

    // vector of nb_reads_seen and nb_kmers_seen 
    let mut nb_kmers_seen_history: Vec<u64> = Vec::new();
    let mut nb_dskmers_seen_history: Vec<u64> = Vec::new();
    // push 0 0 to these vectors
    nb_kmers_seen_history.push(0);
    nb_dskmers_seen_history.push(0);
    let (mut b, mut c, mut d) = (0.0, 0.0, 0.0);
    while let Some(record) = reader.next() {
        nb_reads_seen += 1;
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        let rc = norm_seq.reverse_complement();
        for (_, kmer, _) in norm_seq.canonical_kmers(k as u8, &rc) {
            // for debug purpose, print the k-mer in ascii, transforming the &[u8] to a string
            // println!("Processing k-mer: {:?}", String::from_utf8_lossy(kmer));
                nb_kmers_seen += 1;
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
                    if nb_kmers_seen % step == 0 {

                        // push the current values 
                        nb_kmers_seen_history.push(nb_kmers_seen as u64);
                        nb_dskmers_seen_history.push(unique_solid_kmers as u64);
                        print!(
                            "Processed {} reads, {} kmers seen, {} unique solid kmers. ",
                            nb_reads_seen, nb_kmers_seen, unique_solid_kmers
                        );
                        // fit the model to the last 50 values
                        (b, c, d) = fit_model(&nb_kmers_seen_history, &nb_dskmers_seen_history);

                        // if values are not (0, 0, 0)
                        if b != 0.0 || c != 0.0 || d != 0.0 {
                            println!(
                                "Fitted model: f(x) = {:.3} * x + {:.3} * x^2 + {:.3} * x^3",
                                b, c, d
                            );
                        } else {
                            println!("Not enough data to fit a model.");
                        }
                        if d.abs() < args.stop_acceleration {
                            println!("Stop the computation as acceleration is stable.");
                            break;
                        }
                    }
                }
        

    }
    println!(
        "Final unique k-mers: {}, after processing {} kmers.",
        unique_solid_kmers, nb_kmers_seen
    );


    println!(
        "Average distinct canonical kmers seen per kmer read: {:.2}",
        unique_solid_kmers as f64 / nb_kmers_seen as f64
    );
    let avg_kmers_seen_per_read = nb_kmers_seen as f64 / nb_reads_seen as f64;
    println!(
        "Average total kmers seen per read: {:.2}",
        avg_kmers_seen_per_read
    );

    // special case for gzipped files
    if args.input.extension().map(|e| e == "gz").unwrap_or(false) && args.nb_reads == 0 {
        println!("\x1b[93mInput file is gzipped, cannot estimated total number of kmers.\x1b[0m");
        return Ok(());
    }


	// get the reached position in the file
    let reached_position = read_idx_records_and_get_position(&args.input, nb_reads_seen)?;
    // get the file size
    let file_size = args.input.metadata()?.len();

    
    let mut total_reads = args.nb_reads;

    // if the number of reads is given, we can estimate the total number of distinct canonical kmers
    if total_reads > 0 {
        println!("Total number of reads given: {}", total_reads);
        if nb_reads_seen > total_reads as usize {
            print!("\x1b[93mWarning: idx ({}) is greater than total reads ({}). \x1b[0m", nb_reads_seen, total_reads);
            println!("\x1b[93mI use {} reads as total\x1b[0m", nb_reads_seen);
            nb_reads_seen = nb_reads_seen as usize;
        }
        // we can estimate the total number of kmers based on the number of reads and the reached position
    } else {
        println!(
        "Reached position in file: {} (for {} reads, {} kmers seen), file size: {}",
        reached_position, nb_reads_seen, nb_kmers_seen, file_size
        );
        total_reads = (file_size as f64 / reached_position as f64 * nb_reads_seen as f64) as u64;
        
        println!(
            "Estimated total reads: {}",
            total_reads
        );
    }

    // we have the total number of remaining reads, given that we know the growth we can estimate the total number of remaining kmers
    // we use the last avg growth value to estimate the remaining kmers
    
    // remaining reads
    let remaining_reads = total_reads - nb_reads_seen as u64;
    println!("Remaining reads: {}", remaining_reads);

    
    // Estimate the number of kmers we would have seen for the remaining reads
    let estimated_kmers_remaining = (remaining_reads as f64 * avg_kmers_seen_per_read) as u32;
    println!(
        "Estimated total kmers remaining: {} (based on average total kmers seen per read)",
        estimated_kmers_remaining
    );

    // As we know estimated_kmers_remaining, and we know the growth per step kmers, we can estimate the remaining unique kmers
    // We assume that the stable curve is in the form a +bx +cx^2
    // We can consider only bx + cx^2 for computing the number of expected remaining distinct solid kmers
    
    // b and c cannot be negative, so we take the max of 0 and the values
    b = b.max(0.0);
    c = c.max(0.0);
    // let estimated_remaining_solid_kmers =  estimated_kmers_remaining as f64 * avg_growth as f64 / step as f64;

    let estimated_remaining_solid_kmers = (b * estimated_kmers_remaining as f64
        + c * (estimated_kmers_remaining as f64).powi(2)) as u64;

    println!(
        "Estimated remaining unique k-mers: {:.0} (based on recent curve f(kmers) = {:.4}kmers + {:.4}kmers each kmer)",
        estimated_remaining_solid_kmers, b, c
    );
    let estimated_total_distinct_canonical_kmers = unique_solid_kmers + estimated_remaining_solid_kmers;
    println!("Estimated total distinct canonical k-mers: {:.0}", estimated_total_distinct_canonical_kmers);


    // show the span (|log2(unique_solid_kmers)|)
    println!(
        "Span of distinct canonical k-mers: {:.0}",
        estimated_total_distinct_canonical_kmers.ilog2()
    );

    Ok(())
}

