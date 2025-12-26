use crate::constants;
use crate::io;
use crate::schema;
use nalgebra::{Dyn, Matrix3xX, OVector, Vector3};
use rand;

/// Initialize the system with random positions, velocities, and masses
///
/// # Arguments
/// * `rng` - Random number generator
/// * `n_particles` - Number of particles
/// * `scale` - Scale factor for positions and velocities
///
/// # Returns
/// The initialized state of the system
pub fn init_system<R: rand::Rng + ?Sized>(
    rng: &mut R,
    n_particles: usize,
    scale: f32,
) -> schema::State {
    // nalgebra is column-major, so Matrix3xX is very efficient
    let positions = Matrix3xX::from_fn(n_particles, |_, _| rng.random::<f32>() * scale);
    let velocities = Matrix3xX::from_fn(n_particles, |_, _| {
        (rng.random::<f32>() - 0.5) * 2.0 * scale
    });
    let masses = OVector::<f32, Dyn>::from_fn(n_particles, |_, _| rng.random::<f32>());

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
pub fn get_observables(state: &schema::State) -> schema::Observables {
    let mut kinetic = 0.0;
    let mut potential = 0.0;
    let n_particles = state.positions.ncols();

    // Kinetic energy
    for i in 0..n_particles {
        let v = state.velocities.column(i);
        kinetic += 0.5 * state.masses[i] * v.norm_squared();
    }

    // Potential energy
    for i in 0..n_particles {
        let p_i = state.positions.column(i);
        for j in (i + 1)..n_particles {
            let p_j = state.positions.column(j);
            let dist = (p_j - p_i).norm() + constants::EPS;
            potential -= constants::G * state.masses[i] * state.masses[j] / dist;
        }
    }

    // Stadard deviation of positions
    let mean_pos = state.positions.column_mean();
    let mut pos_var = 0.0;
    for i in 0..n_particles {
        let diff = state.positions.column(i) - &mean_pos;
        pos_var += diff.norm_squared();
    }
    let pos_std = (pos_var / n_particles as f32).sqrt();

    schema::Observables {
        kinetic,
        potential,
        total: kinetic + potential,
        pos_std,
    }
}

/// Advance the system by one time step using the Leapfrog integration method
///
/// # Arguments
/// * `state` - The current state of the system
pub fn step(state: &mut schema::State) {
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
                let diff = p_j - p_i;
                let dist_sq = diff.norm_squared() + constants::EPS;

                // Acceleration = G * m_j * r_ij / |r_ij|^3
                // Note: m_i cancels out because a = F / m_i
                acceleration +=
                    (constants::G * state.masses[j] / (dist_sq * dist_sq.sqrt())) * diff;
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
pub fn simulation(
    state: &mut schema::State,
    n_steps: usize,
    output_traj: &std::path::PathBuf,
    output_obs: &std::path::PathBuf,
    stride: usize,
    burn_in: usize,
) {
    let mut traj_handle = std::fs::File::create(output_traj).unwrap();
    let mut output_handle = std::fs::File::create(output_obs).unwrap();

    io::write_state(state, 0, &mut traj_handle);
    io::write_observables(state, 0, &mut output_handle);

    let pbar = indicatif::ProgressBar::new(n_steps as u64);
    pbar.set_style(
        indicatif::ProgressStyle::with_template(
            "{spinner:.green} {msg:.bold} [{pos}/{len}] [{wide_bar:.cyan/blue}] {eta_precise} [{per_sec}]",
        )
        .unwrap()
        .progress_chars("#>-"),
    );

    for i in 1..=n_steps {
        if (i >= burn_in) && (i % stride == 0) {
            io::write_state(state, i, &mut traj_handle);
            io::write_observables(state, i, &mut output_handle);
        }
        step(state);
        pbar.inc(1);
    }
    pbar.finish_with_message("Simulation complete");
}
