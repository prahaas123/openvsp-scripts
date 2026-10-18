# openvsp-scripts

Python scripts to automate aerodynamic analysis and concept design workflows using [OpenVSP](http://openvsp.org/). All scripts drive OpenVSP through the CLI using generated AngelScript files.

---

## Setup

- Install [OpenVSP 3.47.0+](http://openvsp.org/)
- Create a venv and compile the openvsp Python API
- `pip install neuralfoil`

---

### `vsp_design_space_explore.py`

Design space exploration over the same variables and constraints as `vsp_optimization.py`.

- Sobol samples of the design variables, then one-at-a-time sweeps around the median feasible design
- Reports constraint activity, drag sensitivity per variable, and tightened bounds for the optimizer
- Writes results to `dse_samples.csv` and `dse_oat.csv` (reruns resume from these files)
- Runs VSPAERO cases in parallel, one per CPU core

---

### `vsp_optimization.py`

Single-objective aerodynamic optimization using OpenVSP + VSPAERO as the analysis backend.

- Optimizes root chord, taper, sweep, washout, span and alpha using differential evolution (dihedral is a fixed input)
- Minimizes cruise drag at a fixed airspeed, subject to lift ≥ weight (from the built-in weight model), trim at the target static margin, and AR / tip chord limits
- Profile drag comes from NeuralFoil cd(cl, Re) integrated over the VSPAERO spanwise strips (`PROFILE_DRAG_MODEL = "flat_plate"` switches to a flat-plate estimate)
- Runs VSPAERO cases in parallel, one per CPU core
- Records feasible designs to a timestamped `optimization_results_<date>_<time>.csv`

---

### `tail_sizer.py`

Automated horizontal tail sizing and pitch trim analysis using a secant method iteration loop.

- Sizes the horizontal tail area from a target tail volume coefficient (`V_H`), with `l_H` measured wing AC to tail AC
- Places the CG at the target static margin ahead of the VSPAERO neutral point
- Iteratively adjusts tail incidence angle until the pitching moment coefficient `CMy ≈ 0` at the design alpha
- Writes trim incidence, tail chord and span, tail LE position, NP and CG to `tail_sizing_options.csv`

---

### `vsp_single_design_delta.py`

Runs a full aero + stability analysis on a tailless flying/delta wing configuration.

- Takes in a `Wing.vsp3` custom flying wing geometry as input
- Runs a VSPAERO VLM sweep across an alpha range and writes polar results to `aero_full.csv`
- Runs a stability derivative analysis and writes results to `stability.csv`

---

### `vsp_single_design_conventional.py`

Runs a full aero + stability analysis on a conventional aircraft configuration (wing + horizontal tail + vertical tail).

- Builds a conventional aircraft with a wing, horizontal tail and vertical tail
- Runs a VSPAERO VLM sweep and writes polar results to `aero_full.csv`
- Runs a stability derivative analysis and writes results to `stability.csv`

---

### `constraint_analysis.py`

Generates an aircraft constraint diagram by evaluating the required Thrust-to-Weight ratio (T/W) against Wing Loading (W/S). It calculates and visualizes the performance boundaries for takeoff, cruise, landing and stall speed requirements, ultimately shading the feasible design space for initial aircraft sizing. Takeoff and landing use ground-roll physics with rolling friction and aerodynamic drag.

---

### `airfoil_coords_parser`

Takes 2D cross sectional images of a given wing/body and creates a .dat airfoil file in Selig format. This is useful when trying to recreate complex CAD shapes in OpenVSP.

---

## Configuration

At the top of each script, set:

```python
airfoil_file = os.path.join("Airfoils", "your_airfoil.dat")   # path to airfoil .dat file
```

The `Airfoils/` directory contains `.dat` airfoil coordinate files in Selig format.

Wing and tail geometry parameters are defined as dictionaries near the top of each script and can be freely modified.

---
