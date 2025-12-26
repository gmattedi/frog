use crate::constants;
use crate::io;
use crate::schema;
use nalgebra::{Dyn, Matrix3xX, OVector, Vector3};
use rand;
use rand::{Rng, SeedableRng};
use rand_chacha;

/// Morse potential function and derivative helpers
pub(crate) fn morse_potential(r: f32) -> f32 {
    let e = (-constants::MORSE_A * (r - constants::MORSE_R0)).exp();
    constants::MORSE_D * ((1.0 - e) * (1.0 - e) - 1.0)
}

/// Returns the scalar radial force (signed) = -dV/dr
pub(crate) fn morse_force_mag(r: f32) -> f32 {
    let e = (-constants::MORSE_A * (r - constants::MORSE_R0)).exp();
    -2.0 * constants::MORSE_D * constants::MORSE_A * (1.0 - e) * e
}

/// Apply minimum image convention to a displacement vector for a cubic box
pub(crate) fn minimum_image(diff: &Vector3<f32>, box_size: f32) -> Vector3<f32> {
    Vector3::new(
        diff[0] - box_size * (diff[0] / box_size).round(),
        diff[1] - box_size * (diff[1] / box_size).round(),
        diff[2] - box_size * (diff[2] / box_size).round(),
    )
}

/// Wrap all particle positions in-place into the [0, box_size) cubic box
pub(crate) fn wrap_positions(state: &mut schema::State, box_size: f32) {
    for i in 0..state.positions.ncols() {
        for k in 0..3 {
            let x = state.positions[(k, i)];
            state.positions[(k, i)] = x - box_size * (x / box_size).floor();
        }
    }
}

/// Initialize the system with random positions, velocities, and masses
///
/// # Arguments
/// * `rng` - Random number generator
/// * `n_particles` - Number of particles
/// * `scale_pos` - Scale factor for positions
/// * `scale_vel` - Scale factor for velocities
///
/// # Returns
/// The initialized state of the system
pub fn init_system(config: &schema::InitConfig) -> schema::State {
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(config.seed);

    let positions = Matrix3xX::from_fn(config.n_particles, |_, _| {
        (rng.random::<f32>() - 0.5) * 2.0 * config.scale_pos
    });
    let velocities = Matrix3xX::from_fn(config.n_particles, |_, _| {
        (rng.random::<f32>() - 0.5) * 2.0 * config.scale_vel
    });
    let masses = OVector::<f32, Dyn>::from_fn(config.n_particles, |_, _| {
        rng.random::<f32>() * config.scale_mass
    });

    schema::State {
        positions,
        velocities,
        masses,
    }
}

/// Compute the total kinetic and potential energy of the system
///
/// # Arguments
/// * `state` - The current state of the system
///
/// # Returns
/// The observables containing kinetic, potential, and total energy
pub fn get_observables(
    state: &schema::State,
    periodic: bool,
    box_size: f32,
) -> schema::Observables {
    let mut kinetic_energy = 0.0;
    let mut potential_energy = 0.0;
    let n_particles = state.positions.ncols();

    let positions_mean = state.positions.column_mean();
    let mut positions_var = 0.0;

    for i in 0..n_particles {
        let mut diff;
        // Kinetic energy: 0.5 * m * v^2
        let v = state.velocities.column(i);
        kinetic_energy += 0.5 * state.masses[i] * v.norm_squared();

        // Potential energy: Morse potential between pairs
        let p_i = state.positions.column(i);
        for j in (i + 1)..n_particles {
            let p_j = state.positions.column(j);
            diff = p_j - p_i;
            if periodic {
                // minimum image convention
                diff = minimum_image(&diff, box_size);
            }
            let dist = diff.norm() + constants::EPS; // Avoid division by zero
            // Morse potential: V(r) = D * ( (1 - exp(-a*(r - r0)))^2 - 1 )
            potential_energy += morse_potential(dist);
        }

        let diff = state.positions.column(i) - &positions_mean;
        positions_var += diff.norm_squared();
    }

    let pos_std = (positions_var / n_particles as f32).sqrt();

    schema::Observables {
        kinetic: kinetic_energy,
        potential: potential_energy,
        total: kinetic_energy + potential_energy,
        pos_std,
    }
}

/// Advance the system by one time step using the Leapfrog integration method
///
/// # Arguments
/// * `state` - The current state of the system
pub fn step(state: &mut schema::State, periodic: bool, box_size: f32) {
    let mut new_positions = state.positions.clone();
    let mut new_velocities = state.velocities.clone();
    let n_particles = state.positions.ncols();

    // 1. Update velocities by half step (Leapfrog)
    for i in 0..n_particles {
        let mut acceleration = Vector3::zeros();
        let p_i = state.positions.column(i);

        for j in 0..n_particles {
            if i != j {
                let p_j = state.positions.column(j);
                let mut diff = p_j - p_i;
                if periodic {
                    diff = minimum_image(&diff, box_size);
                }
                let r = diff.norm() + constants::EPS;

                let force_mag = morse_force_mag(r);

                // Acceleration on particle i: a_i = F / m_i (direction along diff)
                acceleration += (diff / r) * (force_mag / state.masses[i]);
            }
        }

        // column_mut(i) is 3x1, acceleration is 3x1.
        // We use & to satisfy the AddAssign trait for views.
        new_velocities
            .column_mut(i)
            .axpy(constants::DT / 2.0, &acceleration, 1.0);
    }

    // 2. Update positions by full step
    for (mut pos_col, vel_col) in new_positions
        .column_iter_mut()
        .zip(new_velocities.column_iter())
    {
        pos_col.axpy(constants::DT, &vel_col, 1.0);
        if periodic {
            for k in 0..3 {
                let x = pos_col[k];
                pos_col[k] = x - box_size * (x / box_size).floor();
            }
        }
    }

    state.positions = new_positions;
    state.velocities = new_velocities;
}

/// Run the N-body simulation for a given number of steps
///
/// # Arguments
/// * `state` - The current state of the system
/// * `n_steps` - Number of simulation steps
/// * `output_traj` - Output file path for trajectory
/// * `output_obs` - Output file path for observables
/// * `stride` - Interval for writing output
/// * `burn_in` - Number of burn-in steps
///
/// # Panics
/// Panics if writing to the output files fails
pub fn simulate(state: &mut schema::State, config: &schema::Config) {
    let mut output_traj = schema::OutputFile {
        file: std::fs::File::create(&config.output_traj)
            .expect("Failed to create trajectory output file"),
        has_header: false,
    };

    let mut output_obs = schema::OutputFile {
        file: std::fs::File::create(&config.output_obs)
            .expect("Failed to create observables output file"),
        has_header: false,
    };

    // If periodic boundaries are enabled, wrap initial positions into the box
    if config.periodic {
        wrap_positions(state, config.box_size);
    }

    let observables = get_observables(state, config.periodic, config.box_size);
    io::write_state(state, 0, &mut output_traj);
    io::write_observables(&observables, 0, &mut output_obs);

    let pbar = indicatif::ProgressBar::new(config.n_steps as u64);
    pbar.set_style(
        indicatif::ProgressStyle::with_template(
            "{spinner:.green} {msg:.bold} [{pos}/{len}] [{wide_bar:.cyan/blue}] {eta_precise} [{per_sec}]",
        )
        .unwrap()
        .progress_chars("#>-"),
    );

    for i in 1..=config.n_steps {
        if (i >= config.burn_in) && (i % config.stride == 0) {
            io::write_state(state, i, &mut output_traj);
            let observables = get_observables(state, config.periodic, config.box_size);
            pbar.set_message(format!("{}", observables));
            io::write_observables(&observables, i, &mut output_obs);
        }
        step(state, config.periodic, config.box_size);
        pbar.inc(1);
    }
    pbar.finish_with_message("Simulation complete");
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{Dyn, Matrix3xX, OVector, Vector3};

    const EPS_F: f32 = 1e-6;

    #[test]
    fn morse_potential_minimum() {
        let r0 = constants::MORSE_R0;
        let v = morse_potential(r0);
        assert!((v + constants::MORSE_D).abs() < EPS_F, "V(r0) should be -D");
    }

    #[test]
    fn morse_force_zero_at_equilibrium() {
        let r0 = constants::MORSE_R0;
        let f = morse_force_mag(r0);
        assert!(f.abs() < EPS_F, "Force at r0 should be zero");
    }

    #[test]
    fn two_particle_potential_matches() {
        let r0 = constants::MORSE_R0;
        // Two particles along x separated by r0
        let p1 = Vector3::new(0.0, 0.0, 0.0);
        let p2 = Vector3::new(r0, 0.0, 0.0);
        let positions = Matrix3xX::from_columns(&[p1, p2]);
        let velocities = Matrix3xX::from_columns(&[Vector3::zeros(), Vector3::zeros()]);
        let masses = OVector::<f32, Dyn>::from_column_slice(&[1.0, 1.0]);

        let state = schema::State {
            positions,
            velocities,
            masses,
        };

        let obs = get_observables(&state, false, 10.0);
        let expected = morse_potential(r0);
        assert!((obs.kinetic).abs() < EPS_F, "Kinetic should be zero");
        assert!(
            (obs.potential - expected).abs() < 1e-5,
            "Potential should match pair Morse potential"
        );
        assert!(
            (obs.total - expected).abs() < 1e-5,
            "Total equals potential"
        );
    }

    #[test]
    fn morse_force_signs() {
        let r0 = constants::MORSE_R0;
        let r_less = r0 * 0.9;
        let r_greater = r0 * 1.1;
        let f_less = morse_force_mag(r_less);
        let f_greater = morse_force_mag(r_greater);
        // r < r0 -> repulsive -> positive force_mag (points along diff)
        assert!(f_less > 0.0, "Force should be repulsive for r < r0");
        // r > r0 -> attractive -> negative force_mag
        assert!(f_greater < 0.0, "Force should be attractive for r > r0");
    }

    #[test]
    fn minimum_image_behavior() {
        let box_size = 5.0;
        let diff = Vector3::new(4.9, 0.0, -4.9);
        let mi = minimum_image(&diff, box_size);
        // 4.9 -> -0.1 after minimum image
        assert!((mi[0] + 0.1).abs() < 1e-6);
        // -4.9 -> 0.1 after minimum image
        assert!((mi[2] - 0.1).abs() < 1e-6);
    }

    #[test]
    fn wrap_positions_initial_state() {
        let mut state = init_system(&schema::InitConfig {
            seed: 123,
            n_particles: 10,
            scale_pos: 1.0,
            scale_vel: 0.1,
            scale_mass: 1.0,
        });
        // Ensure some positions are negative initially
        let mut has_negative = false;
        for i in 0..state.positions.ncols() {
            let x = state.positions[(0, i)];
            if x < 0.0 {
                has_negative = true;
                break;
            }
        }
        assert!(has_negative, "Sanity: initial positions should contain negative values");

        // Wrap positions and verify all are in [0, box)
        let box_size = 10.0;
        wrap_positions(&mut state, box_size);
        for i in 0..state.positions.ncols() {
            for k in 0..3 {
                let v = state.positions[(k, i)];
                assert!(v >= 0.0 && v < box_size, "Position must be inside box after wrapping");
            }
        }
    }
}
