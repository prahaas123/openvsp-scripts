import os
import warnings
import pandas as pd
import numpy as np
import json
import functools

import dash
from dash import dcc, html, Input, Output, ALL
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import neuralfoil as nf
import aerosandbox as asb

# Inputs

AIRFOILS = lambda: search_airfoils_by_geometry(5.0, 15.0, 0.0, 15.0)
REYNOLDS = 100000   # Re for the NeuralFoil polars
TOP_N = 30          # keep only the top N by Cl_max (None = keep all)

AIRFOIL_JSON = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Airfoils", "airfoils.json")

# Core computation functions

def fetch_af_polar(af_name, re_num):
    print(f"  -> COMPUTING: {af_name} at Re={re_num} using NeuralFoil...")
    alphas = np.linspace(-5, 15, 101)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # asb warns for names missing from its database
        airfoil = asb.Airfoil(name=af_name)
    if airfoil.coordinates is None:
        print(f"  -> SKIPPED: {af_name} is not in the AeroSandbox airfoil database.")
        return (None, None, None, None), (None, None)
    try:
        aero = nf.get_aero_from_airfoil(
            airfoil=airfoil,
            alpha=alphas,
            Re=re_num
        )
        cl = aero['CL']
        cd = aero['CD']
        cm = aero['CM']
    except Exception as e:
        print(f"  -> ERROR: Failed to compute {af_name} with NeuralFoil: {e}")
        return (None, None, None, None), (None, None)

    ld = cl / cd
    cl_max = np.max(cl)
    ld_max = np.max(ld)
    return (alphas, cl, cd, cm), [cl_max, ld_max]

def build_airfoil_json(json_filename=AIRFOIL_JSON):
    db_dir = os.path.join(os.path.dirname(asb.__file__), "geometry", "airfoil", "airfoil_database")
    names = sorted(os.path.splitext(f)[0].lower() for f in os.listdir(db_dir) if f.endswith(".dat"))
    print(f"Building '{json_filename}' from {len(names)} AeroSandbox airfoils (one-off)...")
    airfoil_data = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name in names:
            try:
                af = asb.Airfoil(name=name)
                if af.coordinates is None:
                    continue
                thickness = round(100 * float(af.max_thickness()), 1)
                camber = round(100 * float(af.max_camber()), 1)
            except Exception:
                continue
            if np.isfinite(thickness) and np.isfinite(camber):
                airfoil_data[name] = {"max_thickness_percent": thickness, "max_camber_percent": camber}
    os.makedirs(os.path.dirname(json_filename), exist_ok=True)
    with open(json_filename, 'w', encoding='utf-8') as f:
        json.dump(airfoil_data, f, indent=4)
    print(f"  -> Wrote {len(airfoil_data)} airfoils.")
    return airfoil_data

@functools.lru_cache(maxsize=None)
def load_airfoil_data(json_filename=AIRFOIL_JSON):
    try:
        with open(json_filename, 'r', encoding='utf-8') as f:
            return json.load(f)
    except FileNotFoundError:
        return build_airfoil_json(json_filename)
    except json.JSONDecodeError:
        print(f"[WARNING] '{json_filename}' is not valid JSON, rebuilding it.")
        return build_airfoil_json(json_filename)

def search_airfoils_by_geometry(min_thick=0.0, max_thick=66.4, min_camber=0.0, max_camber=16.4, json_filename=AIRFOIL_JSON):
    print(f"Searching for Thickness: {min_thick}-{max_thick}% | Camber: {min_camber}-{max_camber}%...")
    airfoil_names = []
    for name, specs in load_airfoil_data(json_filename).items():
        thickness = specs.get("max_thickness_percent")
        camber = specs.get("max_camber_percent")
        if thickness is not None and camber is not None:
            if (min_thick <= thickness <= max_thick) and (min_camber <= camber <= max_camber):
                airfoil_names.append(name)
    print(f"  -> Found {len(airfoil_names)} airfoils matching the criteria.")
    return airfoil_names

def search_airfoils_by_name(*search_terms, json_filename=AIRFOIL_JSON):
    if not search_terms:
        print("  -> No search terms provided.")
        return []
    print(f"Searching for terms: {', '.join(search_terms)}...")
    terms_lower = [term.lower() for term in search_terms]
    matched_airfoils = [name for name in load_airfoil_data(json_filename)
                        if any(term in name.lower() for term in terms_lower)]
    print(f"  -> Found {len(matched_airfoils)} airfoils matching the search terms.")
    return matched_airfoils

def filter_cm_at_alpha(airfoils, properties, target_alpha, min_cm, max_cm):
    passing_names = set()
    for name, (alphas, cl, cd, cm) in airfoils.items():
        cm_at_target = np.interp(target_alpha, alphas, cm)
        if min_cm <= cm_at_target <= max_cm:
            passing_names.add(name)
    filtered_airfoils = {n: v for n, v in airfoils.items() if n in passing_names}
    filtered_properties = [row for row in properties if row[0] in passing_names]
    print(f"\n  -> {len(filtered_airfoils)} / {len(airfoils)} airfoils passed the Cm filter "
          f"({min_cm} <= Cm <= {max_cm} at alpha={target_alpha}deg).")
    return filtered_airfoils, filtered_properties

def filter_top_clmax(airfoils, properties, top_n=100):
    sorted_properties = sorted(properties, key=lambda x: x[1], reverse=True)
    top_properties = sorted_properties[:top_n]
    top_names = {prop[0] for prop in top_properties}
    filtered_airfoils = {name: data for name, data in airfoils.items() if name in top_names}
    print(f"\n  -> Filtered down to the top {len(filtered_airfoils)} airfoils by Cl_max.")
    return filtered_airfoils, top_properties

# Data loading

def load_and_compute_data(reynolds_number=REYNOLDS):
    names = AIRFOILS() if callable(AIRFOILS) else AIRFOILS
    names = list(dict.fromkeys(n.lower() for n in names))  # dedupe, keep order
    airfoils = {}
    properties = []
    for airfoil in names:
        (alphas, cl, cd, cm), (clmax, ldmax) = fetch_af_polar(airfoil, reynolds_number)
        if alphas is not None:
            airfoils[airfoil] = (alphas, cl, cd, cm)
            properties.append([airfoil, clmax, ldmax])
    if TOP_N is None:
        return airfoils, properties
    return filter_top_clmax(airfoils, properties, TOP_N)

@functools.lru_cache(maxsize=None)
def fetch_and_parse_airfoil_coords(airfoil_name):
    try:
        coords = asb.Airfoil(name=airfoil_name).coordinates
    except Exception:
        return None, None
    if coords is None or len(coords) < 3:
        return None, None
    le_idx = int(np.argmin(coords[:, 0]))
    upper_surface = coords[:le_idx + 1][::-1]   # LE -> TE
    lower_surface = coords[le_idx:]             # LE -> TE
    if len(upper_surface) < 2 or len(lower_surface) < 2:
        return None, None
    return upper_surface, lower_surface

# App initialization

print("Computing airfoil aerodynamics via NeuralFoil...")
all_airfoils, all_properties = load_and_compute_data()

_colors_array = plt.cm.turbo(np.linspace(0, 1, max(len(all_airfoils), 1)))
color_map = {
    name: mcolors.to_hex(_colors_array[i])
    for i, name in enumerate(all_airfoils.keys())
}

# Pre-build every geometry figure once at startup.
# Callbacks will only toggle CSS display — no figure is ever rebuilt on interaction.
print("Pre-fetching airfoil geometry...")
_geo_figures = {}
for _name in all_airfoils.keys():
    _upper, _lower = fetch_and_parse_airfoil_coords(_name)
    if _upper is not None and _lower is not None:
        _color = color_map[_name]
        _r, _g, _b = mcolors.hex2color(_color)
        _fill = f"rgba({int(_r*255)},{int(_g*255)},{int(_b*255)},0.2)"
        _xf = np.concatenate([_upper[:, 0], _lower[::-1, 0]]).tolist()
        _yf = np.concatenate([_upper[:, 1], _lower[::-1, 1]]).tolist()
        _fig = go.Figure()
        _fig.add_trace(go.Scatter(
            x=_xf, y=_yf,
            fill="toself", fillcolor=_fill,
            line=dict(color=_color, width=1.5),
            mode="lines", showlegend=False,
        ))
        _fig.update_layout(
            template="plotly_white",
            paper_bgcolor="#ffffff",
            plot_bgcolor="#f8f9fa",
            title=dict(
                text=_name.upper(),
                font=dict(size=13, color=_color, family="Arial Black"),
                x=0.5, xanchor="center",
            ),
            margin=dict(l=10, r=10, t=36, b=10),
            height=150,
            xaxis=dict(
                showticklabels=False, showgrid=True,
                gridcolor="rgba(255,255,255,0.08)",
                scaleanchor="y", constrain="domain",
            ),
            yaxis=dict(
                showticklabels=False, showgrid=True,
                gridcolor="rgba(255,255,255,0.08)",
                constrain="domain",
            ),
        )
        _geo_figures[_name] = _fig
    else:
        _geo_figures[_name] = None

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.DARKLY])
server = app.server

# Layout helpers

PLOT_CONFIG = {"displayModeBar": False}
POLAR_HEIGHT = 280

PLOTLY_LAYOUT_BASE = dict(
    template="plotly_dark",
    paper_bgcolor="rgba(0,0,0,0)",
    plot_bgcolor="rgba(0,0,0,0)",
    font=dict(size=10),
    margin=dict(l=45, r=10, t=40, b=40),
    showlegend=False,
    xaxis=dict(showgrid=True, gridcolor="rgba(255,255,255,0.15)", gridwidth=0.5),
    yaxis=dict(showgrid=True, gridcolor="rgba(255,255,255,0.15)", gridwidth=0.5),
)

def build_sidebar():
    rows = [
        html.H5("Controls & Legend", className="mb-3 mt-1"),
        html.H6("Airfoil Selection", className="mb-2 text-muted"),
    ]
    for af_name in all_airfoils.keys():
        color = color_map[af_name]
        rows.append(
            html.Div([
                html.Div(style={
                    "width": "18px", "height": "18px",
                    "backgroundColor": color,
                    "borderRadius": "4px",
                    "border": "1px solid rgba(255,255,255,0.2)",
                    "marginRight": "3px", "marginTop": "3px", "flexShrink": "0",
                }),
                dcc.Checklist(
                    id={"type": "airfoil-check", "index": af_name},
                    options=[{"label": f" {af_name.upper()}", "value": af_name}],
                    value=[af_name],
                    inputStyle={"marginRight": "4px", "cursor": "pointer"},
                    labelStyle={"cursor": "pointer", "fontSize": "13px", "color": "white"},
                ),
            ], style={"display": "flex", "alignItems": "flex-start", "marginBottom": "5px"})
        )
    return html.Div(rows, style={
        "backgroundColor": "#16213e",
        "padding": "16px 12px",
        "minHeight": "100vh",
        "borderRight": "1px solid rgba(255,255,255,0.08)",
        "overflowY": "auto",
    })

def polar_graph(graph_id):
    return dcc.Graph(id=graph_id, config=PLOT_CONFIG, style={"height": f"{POLAR_HEIGHT}px"})

def build_geo_panel():
    children = []
    for af_name in all_airfoils.keys():
        fig = _geo_figures.get(af_name)
        if fig is not None:
            content = dcc.Graph(figure=fig, config=PLOT_CONFIG, style={"marginBottom": "8px"})
        else:
            content = dbc.Alert(
                f"Could not retrieve coordinate data for '{af_name}' from the UIUC database.",
                color="warning", className="py-2",
            )
        children.append(html.Div(
            content,
            id={"type": "geo-wrapper", "index": af_name},
            style={"display": "block"},
        ))
    return children

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col(build_sidebar(), width=1, style={"padding": "0"}),

        dbc.Col([
            html.H2("Airfoil Analyzer", className="my-3"),
            dbc.Row([
                # Left half: polars + pareto
                dbc.Col([
                    html.H4("Aerodynamic Polars", className="mb-2"),
                    dbc.Row([
                        dbc.Col(polar_graph("plot-cl-alpha"), width=6),
                        dbc.Col(polar_graph("plot-cd-alpha"), width=6),
                    ], className="g-2 mb-2"),
                    dbc.Row([
                        dbc.Col(polar_graph("plot-cm-alpha"), width=6),
                        dbc.Col(polar_graph("plot-ld-alpha"), width=6),
                    ], className="g-2 mb-2"),
                    dbc.Row([
                        dbc.Col(polar_graph("plot-drag-polar"), width=6),
                    ], className="g-2"),
                    html.Hr(className="my-3"),
                    html.H4("Pareto Front: Lift vs. Efficiency", className="mb-2"),
                    dcc.Graph(id="plot-pareto", config=PLOT_CONFIG, style={"height": "420px"}),
                ], width=6),

                # Right half: all geometry pre-rendered, shown/hidden via callback
                dbc.Col([
                    html.H4("Airfoil Geometry Profiles", className="mb-2"),
                    html.Div(build_geo_panel()),
                ], width=6),
            ], className="g-3"),
        ], width=11),
    ], className="g-0"),
], fluid=True, style={"padding": "0"})

# Callbacks
def _get_selected(values_list):
    selected = []
    for v in values_list:
        if v:
            selected.extend(v)
    return selected

def _make_polar_fig(selected, x_key, y_key, title, xlabel, ylabel):
    fig = go.Figure()
    for name in selected:
        if name not in all_airfoils:
            continue
        alpha, cl, cd, cm = all_airfoils[name]
        ld = cl / cd
        data = {"alpha": alpha, "cl": cl, "cd": cd, "cm": cm, "ld": ld}
        fig.add_trace(go.Scatter(
            x=data[x_key], y=data[y_key],
            mode="lines", name=name.upper(),
            line=dict(color=color_map[name], width=1.5),
        ))
    fig.update_layout(**PLOTLY_LAYOUT_BASE)
    fig.update_layout(
        title=dict(text=title, font=dict(size=11), x=0.5, xanchor="center"),
        xaxis_title=dict(text=xlabel, font=dict(size=9)),
        yaxis_title=dict(text=ylabel, font=dict(size=9)),
    )
    return fig

def _make_pareto_fig(selected):
    filtered_props = [p for p in all_properties if p[0] in selected]
    fig = go.Figure()
    fig.update_layout(**PLOTLY_LAYOUT_BASE)
    if not filtered_props:
        return fig

    df = pd.DataFrame(filtered_props, columns=["Name", "Cl_max", "LD_max"])
    df_sorted = df.sort_values(by=["LD_max", "Cl_max"], ascending=[False, False]).reset_index(drop=True)

    pareto_front = []
    max_cl_seen = -np.inf
    for _, row in df_sorted.iterrows():
        if row["Cl_max"] > max_cl_seen:
            pareto_front.append(row)
            max_cl_seen = row["Cl_max"]
    pareto_df = pd.DataFrame(pareto_front).sort_values(by="LD_max")

    fig.add_trace(go.Scatter(
        x=df["LD_max"], y=df["Cl_max"],
        mode="markers+text",
        marker=dict(color="steelblue", size=8, opacity=0.7),
        text=df["Name"], textposition="top right", textfont=dict(size=9),
        name="Airfoils",
    ))
    fig.add_trace(go.Scatter(
        x=pareto_df["LD_max"], y=pareto_df["Cl_max"],
        mode="lines+markers",
        line=dict(color="crimson", width=2),
        marker=dict(size=8, color="crimson"),
        name="Pareto Front",
    ))
    fig.update_layout(
        showlegend=True,
        legend=dict(font=dict(size=10)),
        xaxis_title="Maximum L/D Ratio (Efficiency)",
        yaxis_title="Maximum Lift Coefficient (Cl_max)",
        margin=dict(l=55, r=15, t=20, b=50),
    )
    return fig

@app.callback(
    Output("plot-cl-alpha",   "figure"),
    Output("plot-cd-alpha",   "figure"),
    Output("plot-cm-alpha",   "figure"),
    Output("plot-ld-alpha",   "figure"),
    Output("plot-drag-polar", "figure"),
    Output("plot-pareto",     "figure"),
    Input({"type": "airfoil-check", "index": ALL}, "value"),
)
def update_plots(values_list):
    selected = _get_selected(values_list)
    return (
        _make_polar_fig(selected, "alpha", "cl", "Cl vs Alpha",           "Alpha (deg)", "Cl"),
        _make_polar_fig(selected, "alpha", "cd", "Cd vs Alpha",           "Alpha (deg)", "Cd"),
        _make_polar_fig(selected, "alpha", "cm", "Cm vs Alpha",           "Alpha (deg)", "Cm"),
        _make_polar_fig(selected, "alpha", "ld", "L/D vs Alpha",          "Alpha (deg)", "L/D"),
        _make_polar_fig(selected, "cd",    "cl", "Drag Polar (Cl vs Cd)", "Cd",          "Cl"),
        _make_pareto_fig(selected),
    )

@app.callback(
    Output({"type": "geo-wrapper", "index": ALL}, "style"),
    Input({"type": "airfoil-check", "index": ALL}, "value"),
)
def update_geometry_visibility(values_list):
    selected = set(_get_selected(values_list))
    return [
        {"display": "block"} if name in selected else {"display": "none"}
        for name in all_airfoils.keys()
    ]

if __name__ == "__main__":
    app.run()