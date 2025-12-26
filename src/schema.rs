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

/// Display implementation for Observables
impl std::fmt::Display for Observables {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "E_tot: {:8.3e}, E_kin: {:8.3e}, E_pot: {:8.3e}, Pos_std: {:8.3e}",
            self.total, self.kinetic, self.potential, self.pos_std
        )
    }
}

/// Initial configuration parameters
///
/// # Fields
/// * `seed` - Random seed for initialization
/// * `n_particles` - Number of particles
/// * `scale_pos` - Scale factor for positions
/// * `scale_vel` - Scale factor for velocities
pub struct InitConfig {
    pub seed: u64,
    pub n_particles: usize,
    pub scale_pos: f32,
    pub scale_vel: f32,
}

/// Configuration for the simulation run
///
/// # Fields
/// * `init_config` - Initial configuration parameters
/// * `n_steps` - Number of simulation steps
/// * `output_traj` - Output file path for trajectory
/// * `output_obs` - Output file path for observables
/// * `stride` - Interval for writing output
/// * `burn_in` - Number of burn-in steps
pub struct Config<'a> {
    pub init_config: &'a InitConfig,
    pub n_steps: usize,
    pub output_traj: &'a std::path::PathBuf,
    pub output_obs: &'a std::path::PathBuf,
    pub stride: usize,
    pub burn_in: usize,
}
