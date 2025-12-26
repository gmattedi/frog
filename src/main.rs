use clap::Parser;
mod constants;
mod engine;
use rand::SeedableRng;
use rand_chacha;
mod cli;
mod io;
mod schema;

/// Entry point of the simulation program.
/// Initializes the simulation state and runs the simulation loop.
fn main() {
    let args = cli::Args::parse();
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(args.seed);

    let mut state = engine::init_system(&mut rng, args.n_particles, 1.0);

    engine::simulation(
        &mut state,
        args.n_steps,
        &args.trajectory,
        &args.observables,
        args.stride,
        args.burn_in,
    );
}
