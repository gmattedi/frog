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
pub struct Observables {
    pub kinetic: f32,
    pub potential: f32,
    pub total: f32,
}
