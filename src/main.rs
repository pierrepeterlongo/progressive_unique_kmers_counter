use clap::Parser; 
use std::path::PathBuf;
use needletail::{parse_fastx_file, Sequence};
use std::error::Error;
use std::usize;
use nalgebra::{DMatrix, DVector, SVD};
use verbose_macros::{debug, verbose};


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
// }
use plotters::prelude::*;
/// plot the accumulated results y in function of x as well as the fitted (y = bx + cx^2) curves in a plot.png file
fn plot_curves(xs: &[f64], ys: &[f64], a: f64, b: f64, c: f64) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new("plot.png", (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("K-mer Growth Curve", ("sans-serif", 30))
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(30)
        .build_cartesian_2d(0f64..*xs.last().unwrap(), 0f64..*ys.last().unwrap())?;

    chart
        .configure_mesh()
        .x_desc("K-mers seen")
        .y_desc("Distinct canonical k-mers")
        .draw()?;

    // Plot the data points
    chart.draw_series(LineSeries::new(
        xs.iter().zip(ys.iter()).map(|(&x, &y)| (x, y)),
        &RED,
    ))?;

    // Plot the fitted curve
    chart.draw_series(LineSeries::new(
        xs.iter().map(|&x| (x, a + b * x + c * x.powi(2))),
        &BLUE,
    ))?;

    root.present()?;
    Ok(())
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

// Fit a quadratic a + bx + cx¨2 polynomial to the given x and y data points using least squares
// Returns the coefficients [a, b, c] if successful, or None if the input is invalid
fn fit_quadratic(xs: &[f64], ys: &[f64]) -> Option<Vec<f64>> {
    let n = xs.len();
    if n != ys.len() || n < 3 {
        return None; // Not enough data points or mismatched lengths
    }

    // Create the design matrix A
    let mut a = DMatrix::zeros(n, 3);
    for (i, &x) in xs.iter().enumerate() {
        a[(i, 0)] = x * x;
        a[(i, 1)] = x;
        a[(i, 2)] = 1.0;
    }

    // Create the observation vector y
    let y = DVector::from_column_slice(ys);

    // Solve the least squares problem using SVD
    let svd = SVD::new(a, true, true);
    let coefficients = svd.solve(&y, 1e-10).ok()?;

    let mut coeffs = coefficients.as_slice().to_vec();
    coeffs.reverse();
    Some(coeffs)
}


fn fit_cubic(xs: &[f64], ys: &[f64]) -> Option<Vec<f64>> {
    let n = xs.len();
    if n != ys.len() || n < 4 {
        return None; // Not enough data points or mismatched lengths
    }

    // Create the design matrix A
    let mut a = DMatrix::zeros(n, 4);
    for (i, &x) in xs.iter().enumerate() {
        a[(i, 0)] = x * x * x;
        a[(i, 1)] = x * x;
        a[(i, 2)] = x;
        a[(i, 3)] = 1.0;
    }

    // Create the observation vector y
    let y = DVector::from_column_slice(ys);

    // Solve the least squares problem using SVD
    let svd = SVD::new(a, true, true);
    let coefficients = svd.solve(&y, 1e-10).ok()?;

    
    let mut coeffs = coefficients.as_slice().to_vec();
    coeffs.reverse();
    Some(coeffs)
}


// UNCOMMENT THIS IF YOU WANT TO USE WEBSOCKET
// #[tokio::main]
// async 
fn main() -> Result<(), Box<dyn std::error::Error>> {


    let args = Args::parse();
    let k = args.k;
    
    let mut unique_kmers: FxHashMap<Vec<u8>, bool> = FxHashMap::default();
    let mut unique_solid_kmers: u32 = 0;

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
    let mut nb_kmers_seen: u32 = 0;
    let mut nb_reads_seen: usize = 0;


    // Initialize the growth history with a single value of 0
    let step = 50000; // step for printing the progress and computing acceleration

    // vector of nb_reads_seen and nb_kmers_seen as arr1
    let mut nb_kmers_seen_history = vec![0.0];
    let mut nb_dskmers_seen_history = vec![0.0];
    let mut a: f64 = 0.0;
    let mut b: f64 = 0.0;
    let mut c: f64 = 0.0;
    while let Some(record) = reader.next() {
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
                nb_kmers_seen += 1;
                match unique_kmers.get_mut(kmer) {
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
                            unique_kmers.insert(kmer.to_vec(), false);
                            debug!(
                                "Inserted new k-mer: {}",
                                String::from_utf8_lossy(kmer)
                            );
                        }
                    }
                }

            if nb_reads_seen % step == 0 {

                // push the current values 
                nb_kmers_seen_history.push(nb_kmers_seen as f64);
                nb_dskmers_seen_history.push(unique_solid_kmers as f64);
                
                    verbose!(
                    "Processed {} reads, {} kmers seen, {} unique solid kmers. ",
                    nb_reads_seen, nb_kmers_seen, unique_solid_kmers
                );
                if nb_kmers_seen_history.len() < 10 {
                    // verbose!("Not enough data to compute acceleration, waiting for more reads...");
                    continue; // Not enough data to compute acceleration
                }
                if unique_solid_kmers < 1000000 {
                    verbose!("Not enough unique solid kmers, waiting for more reads...");
                    continue; // Not enough unique solid kmers to compute acceleration
                }
                // // Create a f64 vector of the last 50 elements of nb_kmers_seen_history
                // // let xs: Vec<f64> = nb_kmers_seen_history[nb_kmers_seen_history.len()-10..]
                // let xs: Vec<f64> = nb_kmers_seen_history[..]
                //     .iter()
                //     .map(|&x| x as f64)
                //     .collect();
                
                
                // // Create a f64 vector of the last 50 elements of nb_dskmers_seen_history
                // // let ys: Vec<f64> = nb_dskmers_seen_history[nb_kmers_seen_history.len()-10..]
                // let ys: Vec<f64> = nb_dskmers_seen_history[..]
                //     .iter()
                //     .map(|&y| y as f64)
                //     .collect();

                let coefs_cube = fit_cubic(&nb_kmers_seen_history[nb_kmers_seen_history.len()-10..], &nb_dskmers_seen_history[nb_kmers_seen_history.len()-10..]);
                if coefs_cube.is_none() {
                    verbose!("Not enough data to compute acceleration, waiting for more reads...");
                    continue; // Not enough data to compute acceleration
                }
                let coefs_cube = coefs_cube.unwrap();
                
                    // print polynomial coefficients
                verbose!("y = {} + {}x + {}x^2 + {}x^3", coefs_cube[0], coefs_cube[1], coefs_cube[2], coefs_cube[3]);
                
                // if coefs_cube[1] < 0.0 {
                //     if coefs_cube[1].abs() > 1e-5 {
                //     verbose!("Non negligible negative linear coefficient detected, waiting for more reads...");
                //     continue; // Negative acceleration, waiting for more reads
                //     }
                // }

                if coefs_cube[2] < 0.0 {
                    if coefs_cube[2].abs() > 1e-5 {
                    verbose!("Non negligible negative acceleration detected, waiting for more reads...");
                    continue; // Negative acceleration, waiting for more reads
                    }
                }
                if coefs_cube[3] < 0.0 {
                    if coefs_cube[3].abs() > 1e-5 {
                        verbose!("Non negligible negative cubic term detected, waiting for more reads...");
                        continue; // Cubic term is negligible, waiting for more reads
                    }
                }

                if coefs_cube[3] < args.stop_acceleration {
                    verbose!("Try to stop the computation as acceleration is stable.");
                    let coefs_quad = fit_quadratic(&nb_kmers_seen_history, &nb_dskmers_seen_history);
                    assert!(coefs_quad.is_some(), "Failed to fit quadratic model");
                    let coefs_quad = coefs_quad.unwrap();
                    // b is the linear coefficient, c is the quadratic coefficient
                    // b and c cannot be negative, so we take the max of 0 and the values
                    a = coefs_quad[0];
                    b = coefs_quad[1];
                    c = coefs_quad[2];
                    verbose!("Poly2 curve found: f(kmers) = {} + {}*kmers + {}*kmers^2", a, b, c);
                    if b < 0.0 {
                            if b.abs() > 1e-5 {
                            verbose!("Non negligible negative linear coefficient detected, waiting for more reads...");
                            continue; // Linear coefficient is negative, waiting for more reads
                        }
                        else {
                            b = b.abs();
                        }
                    }

                    c = coefs_quad[2];

                    if c < 0.0 {
                        if c.abs() > 1e-10 {
                            verbose!("Non negligible negative quadratic coefficient detected, waiting for more reads...");
                            continue; // Quadratic coefficient is negative, waiting for more reads
                        }
                        else {
                            c = c.abs();
                        }
                    }
                    
                    println!("Stable curve found: f(kmers) = {:.4} + {:.4}*kmers + {:.4}*kmers^2", a, b, c);

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
        verbose!(
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
    verbose!("Remaining reads: {}", remaining_reads);

    
    // Estimate the number of kmers we would have seen for the remaining reads
    let estimated_kmers_remaining = (remaining_reads as f64 * avg_kmers_seen_per_read) as u32;
    println!(
        "Estimated total kmers remaining: {} (based on average total kmers seen per read)",
        estimated_kmers_remaining
    );

    // As we know estimated_kmers_remaining, and we know the growth per step kmers, we can estimate the remaining unique kmers
    // We assume that the stable curve is in the form a +bx +cx^2
    // We can consider only bx + cx^2 for computing the number of expected remaining distinct solid kmers
    
   // let estimated_remaining_solid_kmers =  estimated_kmers_remaining as f64 * avg_growth as f64 / step as f64;

    let estimated_remaining_solid_kmers: u32 = (b * (estimated_kmers_remaining as f64)
        + c * (estimated_kmers_remaining as f64).powi(2)) as u32;

    println!(
        "Estimated remaining unique k-mers: {:.0} (based on recent curve f(kmers) = {:.4}kmers + {:.4}kmers^2)",
        estimated_remaining_solid_kmers, b, c
    );
    let estimated_total_distinct_canonical_kmers = unique_solid_kmers + estimated_remaining_solid_kmers;
    println!("Estimated total distinct canonical k-mers: {:.0}", estimated_total_distinct_canonical_kmers);


    // show the span (|log2(unique_solid_kmers)|)
    println!(
        "Span of distinct canonical k-mers: {:.0}",
        estimated_total_distinct_canonical_kmers.ilog2()
    );

    // plot the curves
    if nb_kmers_seen_history.len() > 10 {
        if let Err(e) = plot_curves(&nb_kmers_seen_history, &nb_dskmers_seen_history, a, b, c) {
            eprintln!("Error plotting curves: {}", e);
        } else {
            println!("Plot saved to plot.png");
        }
    } else {
        println!("Not enough data to plot curves.");
    }

    Ok(())
}

