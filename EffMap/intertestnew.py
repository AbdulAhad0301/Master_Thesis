import pandas as pd
import numpy as np
from scipy.interpolate import interp1d, griddata, LinearNDInterpolator
import math

# --- Settings: Operating Envelope & Target Pressure Ratio ---
mass_flow_min = 2
mass_flow_max = 6
design_mass_flow = 4       # kg/s at design point
default_target_pr = 2.5    # Default fixed pressure ratio (or use parameter)

# --- Load Speed Curves & Build Interpolators ---
speed_files = [
    ("SpeedLine1.csv", 1500),
    ("SpeedLine2.csv", 2000),
    ("SpeedLine3.csv", 2500),
    ("SpeedLine4.csv", 3000),
    ("SpeedLine5.csv", 3500),
    ("SpeedLine6.csv", 4000),  # Reference = 100% speed
    ("SpeedLine7.csv", 4500)
]

speed_interp_funcs = {}

for filename, rpm in speed_files:
    df = pd.read_csv(filename)
    df.columns = ['mass_flow', 'pressure_ratio']
    # Use only data within the operating envelope
    df_operating = df[(df['mass_flow'] >= mass_flow_min) & (df['mass_flow'] <= mass_flow_max)]
    if len(df_operating) >= 2:
        # Use linear interpolation; consider spline if data is smooth
        f_interp = interp1d(df_operating['mass_flow'], df_operating['pressure_ratio'],
                            kind='linear', fill_value="extrapolate")
        speed_interp_funcs[rpm] = f_interp

# --- Load Efficiency Map Data ---
eff_files = [
    "Eff0_65.csv",
    "Eff0_70.csv",
    "Eff0_72.csv",
    "Eff0_74.csv",
    "Eff0_75.csv"
]

eff_points = []  # Each row: [mass_flow, pressure_ratio]
eff_values = []  # Efficiency values

for file in eff_files:
    # Extract nominal efficiency from filename (e.g., "Eff0_65.csv" -> 0.65)
    try:
        eff_val = float(file.replace("Eff", "").replace("_", ".").replace(".csv", ""))
    except:
        continue
    df = pd.read_csv(file)
    df.columns = ['mass_flow', 'pressure_ratio']
    for _, row in df.iterrows():
        if not pd.isna(row['mass_flow']) and not pd.isna(row['pressure_ratio']):
            eff_points.append([row['mass_flow'], row['pressure_ratio']])
            eff_values.append(eff_val)
eff_points = np.array(eff_points)
eff_values = np.array(eff_values)

# --- Compressor Model Function ---
def compressor_model_fixed_pr(mass_flow, pressure_ratio_query=default_target_pr):
    """
    Given a mass flow (kg/s) and target pressure ratio, this function:
      - Normalizes the mass flow relative to design_mass_flow.
      - Uses stored speed curve interpolators to gather pressure ratio data at various speeds.
      - Linearly interpolates in the rpm dimension to find the required speed that yields the target pressure ratio.
      - Uses griddata on the efficiency map points to estimate compressor efficiency at (mass_flow, pressure_ratio_query).
      
    If the efficiency is below 50%, it enforces a minimum efficiency of 50%.
    
    Returns:
      rpm_required: the required compressor speed in rpm.
      predicted_eff: predicted efficiency (in fraction, i.e., 0.65 means 65%).
    """
    if mass_flow < mass_flow_min or mass_flow > mass_flow_max:
        raise ValueError(f"Mass flow {mass_flow} is outside the operating envelope [{mass_flow_min}, {mass_flow_max}].")
    
    # Normalize mass flow in kg/s relative to design value
    # (Note: this normalization step is needed if the maps were built on absolute mass flows)
    # For this code we assume the map's x-axis is given in kg/s.
    mf_norm = mass_flow  # Here, mf_norm = mass_flow because the data is assumed in kg/s.
    
    rpms = []
    pr_list = []
    for rpm in sorted(speed_interp_funcs.keys()):
        pr_val = speed_interp_funcs[rpm](mf_norm)
        rpms.append(rpm)
        pr_list.append(pr_val)
    rpms = np.array(rpms)
    pr_list = np.array(pr_list)
    
    # Bracket the target pressure ratio (assume monotonic relationship)
    if pressure_ratio_query < pr_list.min() or pressure_ratio_query > pr_list.max():
        # Extrapolate if needed
        if pressure_ratio_query < pr_list.min():
            i_low, i_high = 0, 1
        else:
            i_low, i_high = -2, -1
    else:
        i_low = None
        for i in range(len(pr_list) - 1):
            if pr_list[i] <= pressure_ratio_query <= pr_list[i+1]:
                i_low = i
                break
        if i_low is None:
            # Try in reverse order if descending
            for i in range(len(pr_list) - 1):
                if pr_list[i] >= pressure_ratio_query >= pr_list[i+1]:
                    i_low = i
                    break
        if i_low is None:
            raise ValueError("Cannot bracket the target pressure ratio with the given data.")
        i_high = i_low + 1
    
    rpm_low = rpms[i_low]
    rpm_high = rpms[i_high]
    pr_low = pr_list[i_low]
    pr_high = pr_list[i_high]
    
    if pr_high == pr_low:
        rpm_required = rpm_low
    else:
        rpm_required = rpm_low + (pressure_ratio_query - pr_low) * (rpm_high - rpm_low) / (pr_high - pr_low)
    
    # Interpolate efficiency using griddata
    query_point = np.array([mass_flow, pressure_ratio_query])
    predicted_eff = griddata(eff_points, eff_values, query_point, method='linear')
    if predicted_eff is None or np.isnan(predicted_eff):
        predicted_eff = griddata(eff_points, eff_values, query_point, method='nearest')
    predicted_eff = float(predicted_eff)
    
    # Enforce a minimum efficiency of 50%
    if predicted_eff < 0.5:
        predicted_eff = 0.5
    
    return rpm_required, predicted_eff

# ----------------------------
# Example Usage
# ----------------------------
test_mass_flow = 4.5  # kg/s
try:
    predicted_rpm, predicted_eff = compressor_model_fixed_pr(test_mass_flow, default_target_pr)
    print(f"For mass flow {test_mass_flow} kg/s and target PR = {default_target_pr}:")
    print(f"  Predicted compressor speed: {predicted_rpm:.1f} rpm")
    print(f"  Predicted efficiency: {predicted_eff * 100:.1f}%")
except ValueError as e:
    print(e)
