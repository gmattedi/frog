# Frog

Basic single-core dynamics simulation with leap-frog integration, in Rust.
The particle act under the Morse potential, and the simulation uses 
Periodic Boundary Conditions (PBC)

## Usage

```
cargo build --release
target/release/frog \
  -n 100 -s 1000000 -d 1000 \
  --box-size 10 --scale-pos 5 --scale-vel 0.1 \
  -t trajectory.txt -o observables.txt
```

## Analysis

You can use [analysis.ipynb](analysis.ipynb) to visualise the trajectory.
It needs to be run in a Python environment with the appropriate dependencies,
which are not defined here to keep things simple.