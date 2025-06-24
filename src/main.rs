pub mod velocity_acceleration;
use clap::Parser;
use std::path::PathBuf;
use needletail::{parse_fastx_file, Sequence};
use std::error::Error;
use std::usize;
use verbose_macros::{debug, verbose};
use bitnuc::as_2bit;

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
    // #[arg(long, default_value_t = 10.0)]
    // stop_acceleration: f64,
    
    /// Give the number of reads
    #[arg(long, default_value_t = 0)]
    nb_reads: usize,
    
    // // UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
    // /// Plot the intermediate results in plot.html
    // #[arg(long, default_value_t = false)]
    // plot: bool,
    
    /// Target number of stable curves to find before stopping
    #[arg(long, default_value_t = 20)]
    nb_stable_curves_found_target: usize,
    
    
    /// Show progress and details
    #[arg(short, long, default_value_t = false)]
    verbose: bool,
    
    /// Show debug information
    #[arg(short, long, default_value_t = false)]
    debug: bool,
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
// // }
// use plotters::prelude::*;
// /// plot the accumulated results y in function of x as well as the fitted (y = bx + cx^2) curve (from last value of xs to maxx) in a plot.png file
// fn _plot_curves(xs: &[f64], ys: &[f64], a: f64, b: f64, c: f64, maxx:u64) -> Result<(), Box<dyn Error>> {
//     let root = BitMapBackend::new("plot.png", (800, 600)).into_drawing_area();
//     root.fill(&WHITE)?;

//     let mut chart = ChartBuilder::on(&root)
//         .caption("K-mer Growth Curve", ("sans-serif", 30))
//         .margin(5)
//         .x_label_area_size(30)
//         .y_label_area_size(30)
//         // .build_cartesian_2d(0f64..*xs.last().unwrap(), 0f64..*ys.last().unwrap())?;
//         .build_cartesian_2d(0f64..*xs.last().unwrap(), 0f64..maxx as f64)?;

//     chart
//         .configure_mesh()
//         .x_desc("K-mers seen")
//         .y_desc("Distinct canonical k-mers")
//         .draw()?;

//     // Plot the data points
//     chart.draw_series(LineSeries::new(
//         xs.iter().zip(ys.iter()).map(|(&x, &y)| (x, y)),
//         &RED,
//     ))?;

//     // Plot the fitted curve
//     let last_xs = *xs.last().unwrap() as u64;
//     chart.draw_series(LineSeries::new(
//         (last_xs..=maxx).map(|x| x as f64).zip(
//             (last_xs..=maxx)
//                 .map(|x| a + b * x as f64 + c * (x as f64).powi(2))
//         ),
//         &BLUE,
//     ))?;


//     root.present()?;
//     Ok(())
// }



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

// UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
// #[tokio::main]
// async 
fn main() -> Result<(), Box<dyn std::error::Error>> {
    
    
    let args = Args::parse();
    let k = args.k;
    if k < 1 || k > 32 {
        return Err("k must be between 1 and 32".into());
    }
    
    let mut unique_kmers: FxHashMap<u64, bool> = FxHashMap::default();
    let mut unique_solid_kmers: usize = 0;
    
    // Set the debug and verbose flags
    verbose_macros::set_debug(args.debug);
    verbose_macros::set_verbose(args.verbose);
    
    
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
    let mut nb_kmers_seen: usize = 0;
    let mut nb_reads_seen: usize = 0;
    
    
    // Initialize the growth history with a single value of 0
    let step = 100000; // step for printing the progress and computing acceleration
    
    // vector of number of kmers seen at each step
    // and number of distinct solid kmers seen at each step
    let mut nb_kmers_seen_history: Vec<usize> = vec![0];
    let mut nb_dskmers_seen_history: Vec<usize> = vec![0];
    // vector of growth values
    let mut growth_history: Vec<f32> = vec![0.0];
    // let mut b: f32 = 0.0;
    // let mut c: f32 = 0.0;
    let mut nb_stable_curves_found: usize = 0;
    // let nb_stable_curves_found_target: usize = 20; // number of stable curves to find before stopping
    // TODO paralelize the processing of reads
    while let Some(record) = reader.next() {
        
        
        if nb_stable_curves_found == args.nb_stable_curves_found_target {
            break;
        }
        debug!(
            "Processed {} reads, {} kmers seen, {} unique solid kmers. ",
            nb_reads_seen, nb_kmers_seen, unique_solid_kmers
        );
        
        nb_reads_seen += 1;
        let seqrec = record.expect("invalid record");
        let norm_seq = seqrec.normalize(false);
        let rc = norm_seq.reverse_complement();
        
        for (_, kmer, _) in norm_seq.canonical_kmers(k as u8, &rc) {
            debug!("Processing k-mer: {:?}", String::from_utf8_lossy(kmer));
            // TODO transform the kmer to 2 bits per base (A=00, C=01, G=10, T=11)
            nb_kmers_seen += 1;
            match unique_kmers.get_mut(&as_2bit(kmer).expect("valid k-mer")) {
                Some(seen) => {
                    if !*seen {
                        *seen = true;
                        unique_solid_kmers += 1;
                        debug!(
                            "Found new unique k-mer: {} (total: {})",
                            String::from_utf8_lossy(kmer),
                            unique_solid_kmers
                        );
                    }
                }
                None => {
                    unique_kmers.insert(as_2bit(kmer).expect("valid k-mer"), false);
                    debug!(
                        "Inserted new k-mer: {}",
                        String::from_utf8_lossy(kmer)
                    );
                }
            }
            if nb_kmers_seen % step == 0 {
                
                // push the current values 
                nb_kmers_seen_history.push(nb_kmers_seen);
                nb_dskmers_seen_history.push(unique_solid_kmers);
                
                verbose!(
                    "Processed {} reads, {} kmers seen, {} unique solid kmers. ",
                    nb_reads_seen, nb_kmers_seen, unique_solid_kmers
                );
                if nb_kmers_seen_history.len() < 10 {
                    // verbose!("Not enough data to compute acceleration, waiting for more reads...");
                    continue; // Not enough data to compute acceleration
                }
                
                let (mut avg_growth, avg_acceleration) = velocity_acceleration::compute_accelerations(
                    &nb_kmers_seen_history,
                    &nb_dskmers_seen_history,
                    10
                );
                
                
                
                // b = avg_growth;
                // c = avg_acceleration;
                growth_history.push(avg_growth);
                
                verbose!("Poly2 curve found: f(kmers) = ? + {}*kmers + {}*kmers^2", avg_growth, avg_acceleration);
                if avg_acceleration > 1e-5 {
                    debug!("Non negligible positive quadratic coefficient detected ({}), waiting for more reads...", avg_acceleration);
                    continue; // Quadratic coefficient is negative, waiting for more reads
                }
                if avg_acceleration < 0.0 && avg_acceleration < -1e-4 {
                    debug!("Non negligible negative quadratic coefficient detected ({}), waiting for more reads...", avg_acceleration);
                    continue; // Quadratic coefficient is negative, waiting for more reads
                }
                
                if avg_growth < 0.0 {
                    if avg_growth.abs() > 1e-4 {
                        debug!("Non negligible negative linear coefficient detected, waiting for more reads...");
                        continue; // Linear coefficient is negative, waiting for more reads
                    }
                    else {
                        avg_growth = avg_growth.abs();
                    }
                }
                nb_stable_curves_found += 1;
                verbose!("Stable curve found {}/{}: f(kmers) = ? + {}*kmers + {}*kmers^2", 
                    nb_stable_curves_found, args.nb_stable_curves_found_target, avg_growth, avg_acceleration);
                
                break;
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
    let avg_kmers_seen_per_read: f64 = nb_kmers_seen as f64 / nb_reads_seen as f64;
    verbose!(
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
    
    
    let mut total_reads: usize = args.nb_reads;
    
    // if the number of reads is given, we can estimate the total number of distinct canonical kmers
    if total_reads > 0 {
        verbose!("Total number of reads given: {}", total_reads);
        if nb_reads_seen > total_reads as usize {
            print!("\x1b[93mWarning: idx ({}) is greater than total reads ({}). \x1b[0m", nb_reads_seen, total_reads);
            println!("\x1b[93mI use {} reads as total\x1b[0m", nb_reads_seen);
            nb_reads_seen = nb_reads_seen as usize;
        }
        // we can estimate the total number of kmers based on the number of reads and the reached position
    } else {
        verbose!(
            "Reached position in file: {} (for {} reads, {} kmers seen), file size: {}",
            reached_position, nb_reads_seen, nb_kmers_seen, file_size
        );
        total_reads = (file_size as f64 / reached_position as f64 * nb_reads_seen as f64) as usize;
        
        println!(
            "Estimated total reads: {}",
            total_reads
        );
    }
    
    // we have the total number of remaining reads, given that we know the growth we can estimate the total number of remaining kmers
    // we use the last avg growth value to estimate the remaining kmers
    
    // remaining reads
    let remaining_reads: usize = total_reads - nb_reads_seen;
    verbose!("Remaining reads: {}", remaining_reads);
    
    
    // Estimate the number of kmers we would have seen for the remaining reads
    let estimated_kmers_remaining: usize = (remaining_reads as f64 * avg_kmers_seen_per_read) as usize;
    verbose!(
        "Estimated total kmers remaining: {} (based on average total kmers seen per read)",
        estimated_kmers_remaining
    );

    println!("Estimated total kmers {}", nb_kmers_seen + estimated_kmers_remaining);
    
    
    // As we know estimated_kmers_remaining, and we know the growth per step kmers, we can estimate the remaining unique kmers
    // We assume that the stable curve is in the form a +bx +cx^2
    // We can consider only bx + cx^2 for computing the number of expected remaining distinct solid kmers
    
    // We consider the used growth as the last 20 values of the growth history
    // TRY LINEAR INTERPOLATION: use a better interpolation method to estimate the averag growth
    use interp::{interp, InterpMode};
    // x is the last 20 values of the nb_kmers_seen_history (or less if there are not enough values)
    let x: Vec<f32> = nb_kmers_seen_history.iter().rev().take(20).map(|&x| x as f32).collect();
    // y is the last 100 values of the unique_solid_kmers history (or less if there are not enough values)
    let y: Vec<f32> = nb_dskmers_seen_history.iter().rev().take(20).map(|&y| y as f32).collect();
    // We can use the interp function to compute the average growth
    let interpolated_number_of_kmers = interp(&x, &y, (nb_kmers_seen + estimated_kmers_remaining) as f32, &InterpMode::default());
    println!(
        "Estimated total distinct canonical using interp: {:.0}",
        interpolated_number_of_kmers
    );

    // // TRY POLYNOMIAL INTERPOLATION: use a polynomial interpolation to estimate the average growth
    // use polynominal_interpolation;
    // // create a vector xs composed of nb_kmers_seen_history as f64
    // let xs: Vec<f64> = nb_kmers_seen_history.iter().rev().take(20).map(|&x| x as f64).collect::<Vec<_>>().into_iter().rev().collect();
    // // create a vector ys composed of nb_dskmers_seen_history as f64
    // let ys: Vec<f64> = nb_dskmers_seen_history.iter().rev().take(20).map(|&x| x as f64).collect::<Vec<_>>().into_iter().rev().collect();
    // let f = polynominal_interpolation::newton_interpolation(xs, ys);
    // let estimated_distinct_kmers = f((nb_kmers_seen + estimated_kmers_remaining) as f64);
    // println!(
    //     "Estimated total distinct canonical using polynomial interpolation: {:.0}",
    //     estimated_distinct_kmers
    // );


    let used_avg_growth: f32 = growth_history.iter().rev().take(20).sum::<f32>() / 10.0;

    let estimated_remaining_solid_kmers: usize = (used_avg_growth * estimated_kmers_remaining as f32) as usize;
    
    verbose!(
        "Estimated remaining unique k-mers: {:.0} (based on recent curve f(kmers) = {}*kmers)",
        estimated_remaining_solid_kmers, used_avg_growth
    );
    let estimated_total_distinct_canonical_kmers = unique_solid_kmers + estimated_remaining_solid_kmers;
    println!(
        "Estimated total distinct canonical k-mers: {:.0}",
        estimated_total_distinct_canonical_kmers
    );
    
    // show the span (|log2(unique_solid_kmers)|)
    println!(
        "Span of distinct canonical k-mers: {:.0}",
        estimated_total_distinct_canonical_kmers.ilog2()
    );
    
    Ok(())
}
