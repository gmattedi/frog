use clap::Parser;
use std::path::PathBuf;

fn get_styles() -> clap::builder::Styles {
    clap::builder::Styles::styled()
        .usage(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Cyan))),
        )
        .header(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Cyan))),
        )
        .literal(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .invalid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .error(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Red))),
        )
        .valid(
            anstyle::Style::new()
                .bold()
                .fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::Green))),
        )
        .placeholder(
            anstyle::Style::new().fg_color(Some(anstyle::Color::Ansi(anstyle::AnsiColor::White))),
        )
}

/// Simple N-body gravitational simulator
///
/// # Arguments
/// * `n_particles` - Number of particles in the simulation
/// * `n_steps` - Number of simulation steps
/// * `trajectory` - Path to trajectory output file
/// * `observables` - Path to observables output file
/// * `seed` - Random seed
/// * `burn_in` - Number of burn-in steps
/// * `stride` - Interval for writing output
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(author, version, about, long_about = None, verbatim_doc_comment)]
#[command(name = "frog")]
#[command(styles=get_styles())]
pub struct Args {
    // Random seed
    #[arg(short = 'r', long, default_value_t = 42, help = "Random seed")]
    pub seed: u64,
    // Number of particles
    #[arg(short = 'n', long, help = "Number of particles in the simulation")]
    pub n_particles: usize,
    // Number of simulation steps
    #[arg(short = 's', long, help = "Number of simulation steps")]
    pub n_steps: usize,
    // Scale for positions
    #[arg(
        long,
        default_value_t = 1.0,
        help = "Scale factor for initial positions"
    )]
    pub scale_pos: f32,
    // Scale for velocities
    #[arg(
        long,
        default_value_t = 1e-4,
        help = "Scale factor for initial velocities"
    )]
    pub scale_vel: f32,
    // Scale for masses
    #[arg(long, default_value_t = 1.0, help = "Scale factor for initial masses")]
    pub scale_mass: f32,
    // Burn-in steps
    #[arg(
        short = 'b',
        long,
        default_value_t = 0,
        help = "Number of burn-in steps"
    )]
    pub burn_in: usize,
    // Output file path
    #[arg(
        short,
        long,
        default_value = "trajectory.txt",
        help = "Trajectory output file path"
    )]
    pub trajectory: PathBuf,
    // Energy file path
    #[arg(short, long, help = "Observables output file path")]
    pub observables: PathBuf,
    // Stride for output
    #[arg(short = 'd', long, default_value_t = 1, help = "Stride for output")]
    pub stride: usize,
    // Periodic boundary conditions
    #[arg(
        long,
        default_value_t = true,
        help = "Enable periodic boundary conditions"
    )]
    pub periodic: bool,
    // Cubic box size for PBC
    #[arg(
        long,
        default_value_t = 10.0,
        help = "Box size for periodic boundary conditions"
    )]
    pub box_size: f32,
}
