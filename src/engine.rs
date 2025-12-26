use crate::constants;
use crate::io;
use crate::schema;
use nalgebra::{Dyn, Matrix3xX, OVector, Vector3};
use rand;
use rand::{Rng, SeedableRng};
use rand_chacha;

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
pub fn get_observables(state: &schema::State) -> schema::Observables {
    let mut kinetic_energy = 0.0;
    let mut potential_energy = 0.0;
    let n_particles = state.positions.ncols();

    let positions_mean = state.positions.column_mean();
    let mut positions_var = 0.0;

    for i in 0..n_particles {
        // Kinetic energy: 0.5 * m * v^2
        let v = state.velocities.column(i);
        kinetic_energy += 0.5 * state.masses[i] * v.norm_squared();

        // Potential energy: Morse potential between pairs
        let p_i = state.positions.column(i);
        for j in (i + 1)..n_particles {
            let p_j = state.positions.column(j);
            let diff = p_j - p_i;
            let dist = diff.norm() + constants::EPS; // Avoid division by zero
            // Morse potential: V(r) = D * ( (1 - exp(-a*(r - r0)))^2 - 1 )
            let e = (-constants::MORSE_A * (dist - constants::MORSE_R0)).exp();
            let v = constants::MORSE_D * ((1.0 - e) * (1.0 - e) - 1.0);
            potential_energy += v;
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
                let r = diff.norm() + constants::EPS;

                // Morse force magnitude: F = -dV/dr, where
                // dV/dr = 2 * D * a * (1 - exp(-a*(r - r0))) * exp(-a*(r - r0))
                let e = (-constants::MORSE_A * (r - constants::MORSE_R0)).exp();
                let d_vdr = 2.0 * constants::MORSE_D * constants::MORSE_A * (1.0 - e) * e;
                let force_mag = -d_vdr;

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

    let observables = get_observables(state);
    io::write_state(state, 0, &mut output_traj, config.center_trajectory);
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
            io::write_state(state, i, &mut output_traj, config.center_trajectory);
            let observables = get_observables(state);
            pbar.set_message(format!("{}", observables));
            io::write_observables(&observables, i, &mut output_obs);
        }
        step(state);
        pbar.inc(1);
    }
    pbar.finish_with_message("Simulation complete");
}
