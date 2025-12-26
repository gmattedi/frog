use crate::schema;
use std::io::Write;

/// Write the current state to the given file
///
/// # Arguments
/// * `state` - The current state of the system
/// * `step` - The current simulation step
/// * `file` - The file to write to
///
/// # Panics
/// Panics if writing to the file fails
pub fn write_state(state: &schema::State, step: usize, file: &mut std::fs::File) {
    for i in 0..state.positions.ncols() {
        let pos = state.positions.column(i);
        let mass = state.masses[i];
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            step, i, pos[0], pos[1], pos[2], mass
        )
        .unwrap();
    }
}

/// Write the current observables to the given file
///
/// # Arguments
/// * `observables` - The current observables of the system
/// * `step` - The current simulation step
/// * `file` - The file to write to
///
/// # Panics
/// Panics if writing to the file fails
pub fn write_observables(observables: &schema::Observables, step: usize, file: &mut std::fs::File) {
    writeln!(
        file,
        "{}\t{}\t{}\t{}\t{}",
        step, observables.kinetic, observables.potential, observables.total, observables.pos_std
    )
    .unwrap();
}
