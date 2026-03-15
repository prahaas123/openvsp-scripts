# openvsp-scripts

Python scripts to automate aerodynamic analysis and concept design workflows using [OpenVSP](http://openvsp.org/). All scripts drive OpenVSP through the CLI using generated AngelScript files.

---

## OpenVSP Setup

- Install [OpenVSP 3.47.0+](http://openvsp.org/)
- Set `vsp_exe` path at the top of each script

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
- Visualizes the trimmed configuration

---

### `vsp_single_design_delta.py`

Runs a full aero + stability analysis on a tailless delta wing configuration.

- Builds a delta wing geometry with a custom airfoil (DAE-21) and elevons
- Runs a VSPAERO VLM sweep across an alpha range and writes polar results to `cfd_sweep.csv`
- Runs a stability derivative analysis and writes results to `vsp_derivatives.csv`
- Visualizes the generated STL using PyVista

---

### `vsp_single_design_conventional.py`

Runs a full aero + stability analysis on a conventional aircraft configuration (wing + horizontal tail + vertical tail).

- Builds a three-surface model with a wing, tails, and ailerons/elevators/rudder control surfaces
- Runs a VSPAERO VLM sweep and writes polar results to `cfd_sweep.csv`
- Runs a stability derivative analysis and writes results to `vsp_derivatives.csv`
- Visualizes the generated STL using PyVista

---

### `scoring_sensititivity.py`

Sensitivity analysis script — evaluates how changes in individual design parameters affect a scoring metric.

---

## How It Works

Each script generates temporary `.vspscript` files (AngelScript) on the fly, executes them via:

```Powershell
vsp.exe -script <script_name>.vspscript
```

and then parses the output files (`.polar`, `.stab`) directly. Temporary scripts and intermediate VSP files are deleted automatically after each run.

---

## Configuration

At the top of each script, set:

```python
vsp_exe    = r"C:\Program Files\OpenVSP-3.47.0\vsp.exe"  # path to your vsp.exe
airfoil_file = r"path\to\your\airfoil.dat"               # path to airfoil .dat file
```

The `Airfoils/` directory contains `.dat` airfoil coordinate files in Selig format.

Wing and tail geometry parameters are defined as dictionaries near the top of each script and can be freely modified.

---
