use nalgebra::{Dyn, Matrix3xX, OVector};

/// State of the N-body system
///
/// # Fields
/// * `positions` - 3xN matrix of particle positions
/// * `velocities` - 3xN matrix of particle velocities
/// * `masses` - N-dimensional vector of particle masses
pub struct State {
    // 3 rows (x, y, z), N columns (particles)
    pub positions: Matrix3xX<f32>,
    pub velocities: Matrix3xX<f32>,
    pub masses: OVector<f32, Dyn>,
}

/// Observables of the N-body system
///
/// # Fields
/// * `kinetic` - Total kinetic energy
/// * `potential` - Total potential energy
/// * `total` - Total energy
/// * `pos_std` - Standard deviation of particle positions
pub struct Observables {
    pub kinetic: f32,
    pub potential: f32,
    pub total: f32,
    pub pos_std: f32,
}

/// Configuration for the simulation run
///
/// # Fields
/// * `n_steps` - Number of simulation steps
/// * `output_traj` - Output file path for trajectory
/// * `output_obs` - Output file path for observables
/// * `stride` - Interval for writing output
/// * `burn_in` - Number of burn-in steps
pub struct SimulateConfig<'a> {
    pub n_steps: usize,
    pub output_traj: &'a std::path::PathBuf,
    pub output_obs: &'a std::path::PathBuf,
    pub stride: usize,
    pub burn_in: usize,
}
