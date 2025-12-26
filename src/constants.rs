/// Softening parameter to avoid singularities
pub const EPS: f32 = 1e-6;
/// Time step for the simulation
pub const DT: f32 = 1e-3;

/// Morse potential depth (D_e)
pub const MORSE_D: f32 = 0.1;
/// Morse potential width parameter (a)
pub const MORSE_A: f32 = 50.0;
/// Morse equilibrium distance (r_e)
pub const MORSE_R0: f32 = 0.1;
