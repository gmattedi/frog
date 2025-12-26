use clap::Parser;
mod cli;
mod constants;
mod engine;
mod io;
mod schema;

/// Entry point of the simulation program.
/// Initializes the simulation state and runs the simulation loop.
fn main() {
    let args = cli::Args::parse();

    let init_config = schema::InitConfig {
        seed: args.seed,
        n_particles: args.n_particles,
        scale_pos: args.scale_pos,
        scale_vel: args.scale_vel,
    };

    let config = schema::Config {
        init_config: &init_config,
        output_traj: &args.trajectory,
        output_obs: &args.observables,
        n_steps: args.n_steps,
        stride: args.stride,
        burn_in: args.burn_in,
    };

    let mut state = engine::init_system(&config.init_config);

    engine::simulate(&mut state, &config);
}
