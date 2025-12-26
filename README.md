# Frog

Basic single-core dynamics simulation with leap-frog integration, in Rust.

## Usage

```
cargo build --release
target/release/frog \
  -n 100 -s 11000 -b 1000 -d 100 \
  -t trajectory.txt -o observables.txt
```

## Analysis

You can use [analysis.ipynb](analysis.ipynb) to visualise the trajectory.
It needs to be run in a Python environment with the appropriate dependencies,
which are not defined here to keep things simple.