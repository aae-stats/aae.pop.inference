# aae.pop.inference (development version)

## API changes

- Arguments to `simulate` are now passed as `sim_args` (a list) in `define`
- Dots argument in `define` passed to `stat` rather than `simulate`

# aae.pop.inference 0.1.0

## Features

- Basic implementation of SMC-ABC for population models created with aae.pop: `inference`
- Helpers to calculate summary statistics from replicated population trajectories: `stat_abundance_trend`, `stat_abundance_moment`
- Helper to create a simulation function for `inference` from an existing matrix population model defined with aae.pop: `define`

