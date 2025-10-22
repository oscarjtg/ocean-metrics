# ocean-metrics
Functions for computing metrics from Oceananigans.jl model netCDF output.

[Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) is a GPU-accelerated ocean model in the Julia programming language. 
It solves the Navier Stokes equations and tracer transport equations 
using a finite volume method on a staggered grid.
Oceananigans.jl can output data as netCDF files using its `NetCDFWriter`.
This Python module takes these output files and computes diagnostics 
and metrics from them.

# Installation instructions

## (1) Use module directly

Clone the repository with 
```
git clone https://github.com/oscarjtg/ocean-metrics.git
```
and either add the directory to your Python path 
or copy the `oceanmetrics.py` module into your project.

## (2) Install in a virtual environment

Clone the repository into your Python `site-packages` directory, 
activate a virtual environment, navigate to the ocean-metrics directory,
and type
```
pip install -e .
```
This will install the package in your virtual environment in editable mode.

# TODO

* Compute vorticity
* Compute speed
* Compute errors with respect to known fields
* Compute tendencies (time derivatives)