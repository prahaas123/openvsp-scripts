# openvsp-scripts

Python scripts to automate aerodynamic analysis and concept design workflows using [OpenVSP](http://openvsp.org/) and VSPAERO.

## Setup

- Install [OpenVSP 3.47.0+](http://openvsp.org/)
- Create a venv and install the OpenVSP Python API into it
- `pip install numpy scipy pandas matplotlib plotly neuralfoil aerosandbox dash dash-bootstrap-components opencv-python`

Airfoil `.dat` files (Selig format) go in `Airfoils/`. Each script picks one with `airfoil_file` or `AIRFOIL_FILE` near the top.

## Scripts

### `constraint_analysis.py`

Draws a T/W vs. W/S constraint diagram from cruise, takeoff ground roll and stall speed requirements, and shades the feasible region.

- **Edit:** `V_cruise`, `s_to`, `V_stall_max`, `C_Lmax`, `C_D0`, `AR`, `e`, rolling friction `mu`.

### `airfoil_analyzer.py`

A local Dash web app for comparing airfoils. It computes polars for a filtered set of airfoils and shows Cl, Cd, Cm and L/D vs. alpha, the drag polar, a Cl_max vs. L/D Pareto plot and the airfoil shapes.

- **Edit:** `AIRFOILS` (filter by thickness and camber with `search_airfoils_by_geometry`, or by name with `search_airfoils_by_name`), `REYNOLDS`, `TOP_N`.
- **Run:** run the script, then open the URL it prints.
- **Output:** results stay in the browser. Use it to narrow airfoil choices.

### `airfoil_coords_parser.py`

Traces a 2D cross-section image of a wing or body and writes a Selig-format `.dat` airfoil file. Useful for recreating CAD shapes in OpenVSP.

- **Edit:** `INPUT_IMAGE`, `OUTPUT_NAME`, `NUM_POINTS`, `THRESHOLD`.

### `vsp_design_space_explore.py`

Run this before `vsp_optimization.py`. It Sobol-samples the design space, runs VSPAERO on each sample in parallel and prints tightened bounds that cover the feasible designs.

- **Edit:**
  - `N_SAMPLES`, `KEEP_PCTL`
  - `VELOCITY` (cruise speed), `AIRFOIL_FILE`, `DIHEDRAL` (fixed)
  - `BOUNDS` for the 6 design variables: root chord, taper, LE sweep, washout, span, alpha
  - the weight model (`M_FIXED`, `M_SPAR_PER_M`, `M_SKIN_PER_M2`)
  - constraint targets (`STATIC_MARGIN`, `CM_TOL`, `AR_MIN`/`AR_MAX`, `TIP_CHORD_MIN`)
- **Output:** `dse_samples.csv` and a pair plot of the feasible region. Copy the tightened bounds into `BOUNDS` before optimizing.

### `vsp_optimization.py`

Minimizes cruise drag with differential evolution, subject to lift ≥ weight, trim at the target static margin, and the AR and tip chord limits.

- **Edit:** `POPSIZE`, `MAXITER`, `N_WORKERS`, and `PROFILE_DRAG_MODEL` (`"strip"` or `"flat_plate"`). Bounds and constraints are edited in `vsp_design_space_explore.py`.
- **Output:** every feasible design goes to `optimization_results_<date>_<time>.csv`. The best wing is printed and exported as an STL.
- **Stop:** `Ctrl+C` stops the run cleanly and still reports the best feasible design logged so far.

### `tail_sizer.py`

Sizes the horizontal and vertical tails for a conventional layout using tail volume coefficients (`V_H`, `V_V`). For every combination of wing config × tail moment arm × horizontal-tail aspect ratio, it places the CG at the target static margin (`SM`) using the VSPAERO neutral point. It then adjusts the tail incidence with a secant search in VSPAERO until the pitching moment is trimmed.

- **Edit:**
  - `wing_configs` (one dict per candidate wing)
  - `tail_moment_arms`, `htail_aspect_ratios`
  - `htail_params`, `vtail_params`
  - `velocity`, `alpha`, `SM`, `airfoil_file`
- **Output:** `tail_sizing_options.csv`

### `vsp_single_design_conventional.py`

Builds one conventional wing + H-tail + V-tail aircraft and analyzes it in full: an alpha sweep at each velocity and a stability derivative run.

- **Edit:**
  - `wing_params`, `htail_params`, `vtail_params`
  - `x_cg`, `airfoil_file`, `stability_velocity`
  - `velocities`, `alphas`, `stability_targets`
- **Output:** `aero_full.csv`, `stability.csv` and an interactive Plotly dashboard that opens in your browser.

### `vsp_single_design_delta.py`

The same full analysis as the conventional script, for a tailless swept/delta flying wing.

- **Edit:** `wing_params` (root chord, taper, sweep, washout, span, dihedral, alpha), plus `x_cg`, `airfoil_file`, `stability_velocity`, `velocities`, `alphas`.
- **Output:** `aero_full.csv`, `stability.csv` and an interactive Plotly dashboard that opens in your browser.

### Stability (no OpenVSP needed)

All three read the `aero_full.csv` and `stability.csv` from a single-design script, in the same folder.

#### `3dof_longitudinal_stability.py`

It then applies a pitch disturbance from the trim conditions and integrates the nonlinear 3-DOF longitudinal equations, using `Cm_q` from `stability.csv`.

- **Edit:** `mass`, inertia properties, and the initial disturbance.
- **Output:** time history plots.

#### `3dof_lateral_stability.py`

Simulates the lateral-directional response to a at `V_trim`. It interpolates the derivatives by velocity from `stability.csv`.

- **Edit:** `mass`, `I_x`, `I_z`, `I_xz`, `V_trim`, and the initial disturbance.
- **Output:** time history plots.

#### `dynamic_modes.py`

Finds the trim condition, then builds the linearized longitudinal and lateral state-space matrices and computes their eigenvalues. It plots them alongside the classic approximations: short period, phugoid, Dutch roll, roll and spiral.

- **Edit:** `mass`, `Ix`, `Iy`, `Iz`, `Ixz`, and the assumed `CLq`, `CDq`.
- **Output:** trim state and eigenvalues printed to the terminal and a root-locus style plot. Left half-plane means stable.
