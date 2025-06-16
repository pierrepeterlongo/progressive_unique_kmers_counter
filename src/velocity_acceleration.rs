use std::f32::NAN;

pub fn compute_accelerations(x_values: &[usize], y_values: &[usize], window_size: usize) -> (f32, f32) {
    let n = x_values.len();
    let mut velocities = Vec::with_capacity(n - 1);
    let mut accelerations = Vec::with_capacity(n - 2);
    if window_size < 2 {
        panic!("Window size must be at least 2 to compute accelerations.");
    }
    if n < window_size || y_values.len() != n || x_values.len() != n {
        panic!("Input vectors must have the same length and at least {} points.", window_size);
    }
    // Calculate velocities
    for i in n-window_size+1..n {
        if x_values[i] - x_values[i - 1] == 0 {
            return (NAN, NAN); // Avoid division by zero
        }
        let velocity = (y_values[i] - y_values[i - 1]) as f32 / (x_values[i] - x_values[i - 1]) as f32;
        velocities.push(velocity);
    }
   

    let sum = velocities.iter().sum::<f32>();
    let average_velocity = sum / velocities.len() as f32;

    // Calculate accelerations
    for i in 1..velocities.len() {
        let acceleration = velocities[i] - velocities[i - 1];
        accelerations.push(acceleration);
    }

    let sum = accelerations.iter().sum::<f32>();
    let average_acceleration =  sum / accelerations.len() as f32;
    (average_velocity, average_acceleration)
}
