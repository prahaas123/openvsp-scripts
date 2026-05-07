# openvsp-scripts

Python scripts to automate aerodynamic analysis and concept design workflows using [OpenVSP](http://openvsp.org/). All scripts drive OpenVSP through the CLI using generated AngelScript files.

---

## Setup

- Install [OpenVSP 3.47.0+](http://openvsp.org/)
- Create a venv and compile the openvsp Python API

---

### `vsp_design_space_explore.py`

Design space exploration — Monte Carlo analysis of the various possible parameter configurations.

- Varies geometric parameters (span, chord, sweep, taper, twist, etc.) across a defined grid
- Runs VSPAERO for each configuration and records CL, CD, and L/D
- Writes all results to `sweep_results.csv`
- Visualizes results using pareto fronts, parallel coordinate plots, and splom plots

---

### `vsp_optimization.py`

Single-objective aerodynamic optimization using OpenVSP + VSPAERO as the analysis backend.

- Optimizes wing geometry parameters using the differetial evolution method
- Maximizes a performance metric (can just be L/D or a weighted sum of various objectives)
- Records the full optimization history to `optimization_results.csv`

---

### `tail_sizer.py`

Automated horizontal tail sizing and pitch trim analysis using a secant method iteration loop.

- Sizes the horizontal tail area from a target tail volume coefficient (`V_H`)
- Iteratively adjusts tail incidence angle until the pitching moment coefficient `CMy ≈ 0` at the design alpha
- Reports final trim incidence, tail chord, and tail span

---

### `vsp_single_design_delta.py`

Runs a full aero + stability analysis on a tailless flying/delta wing configuration.

- Takes in a `Wing.vsp3` custom flying wing geometry as input
- Runs a VSPAERO VLM sweep across an alpha range and writes polar results to `aero_full.csv`
- Runs a stability derivative analysis and writes results to `stability.csv`

---

### `vsp_single_design_conventional.py`

Runs a full aero + stability analysis on a conventional aircraft configuration (wing + horizontal tail + vertical tail).

- Builds a conventional aircraft with a wing, tails, and ailerons/elevators/rudder control surfaces
- Runs a VSPAERO VLM sweep and writes polar results to `aero_full.csv`
- Runs a stability derivative analysis and writes results to `stability.csv`

---

### `constraint_analysis.py`

Generates an aircraft constraint diagram by evaluating the required Thrust-to-Weight ratio (T/W) against Wing Loading (W/S). It calculates and visualizes the performance boundaries for takeoff, cruise, and landing requirements, ultimately shading the feasible design space for initial aircraft sizing.

---

## Configuration

At the top of each script, set:

```python
airfoil_file = r"path\to\your\airfoil.dat"               # path to airfoil .dat file
```

The `Airfoils/` directory contains `.dat` airfoil coordinate files in Selig format.

Wing and tail geometry parameters are defined as dictionaries near the top of each script and can be freely modified.

---
