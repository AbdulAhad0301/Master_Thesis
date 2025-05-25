import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, griddata

# ---------------------------
# Settings: Operating Envelope & Fixed Pressure Ratio
# ---------------------------
mass_flow_min = 2
mass_flow_max = 6
design_mass_flow = 4
fixed_pressure_ratio = 2.5

# ---------------------------
# Load Speed Curves and Build Interpolation Functions
# ---------------------------
# Each tuple is (filename, rpm)
speed_files = [
    ("SpeedLine1.csv", 1500),
    ("SpeedLine2.csv", 2000),
    ("SpeedLine3.csv", 2500),
    ("SpeedLine4.csv", 3000),
    ("SpeedLine5.csv", 3500),
    ("SpeedLine6.csv", 4000),
    ("SpeedLine7.csv", 4500)
]

# Dictionary mapping rpm to an interpolation function for pressure ratio vs mass flow.
speed_interp_funcs = {}

for filename, rpm in speed_files:
    df = pd.read_csv(filename)
    df.columns = ['x', 'y']  # x: mass flow, y: pressure ratio
    # Consider only points within the operating envelope.
    df_operating = df[(df['x'] >= mass_flow_min) & (df['x'] <= mass_flow_max)]
    if len(df_operating) >= 2:
        f_interp = interp1d(df_operating['x'], df_operating['y'], kind='linear', fill_value="extrapolate")
        speed_interp_funcs[rpm] = f_interp

# ---------------------------
# Load Efficiency Curves and Prepare Efficiency Data Points
# ---------------------------
# Efficiency files: the filename like "Eff0_65.csv" indicates an efficiency of 0.65.
eff_files = [
    "Eff0_65.csv",
    "Eff0_70.csv",
    "Eff0_72.csv",
    "Eff0_74.csv",
    "Eff0_75.csv"
]

# We'll gather all efficiency curve points into arrays.
eff_points = []  # Each point is [mass_flow, pressure_ratio]
eff_values = []  # Corresponding efficiency

for file in eff_files:
    # Extract efficiency value from filename: "Eff0_65.csv" -> 0.65
    eff_val = float(file.replace("Eff", "").replace("_", ".").replace(".csv", ""))
    df = pd.read_csv(file)
    df.columns = ['x', 'y']
    for _, row in df.iterrows():
        if not pd.isna(row['x']) and not pd.isna(row['y']):
            eff_points.append([row['x'], row['y']])
            eff_values.append(eff_val)
eff_points = np.array(eff_points)
eff_values = np.array(eff_values)

# ---------------------------
# Define the Compressor Model Function for Fixed Pressure Ratio
# ---------------------------
def compressor_model_fixed_pr(mass_flow, pr_fixed=fixed_pressure_ratio):
    """
    Given a mass flow (in m.sqrt(T)/p), this function finds the compressor speed (rpm)
    and efficiency such that the operating pressure ratio is fixed at pr_fixed.
    
    The function:
      - Checks the operating envelope (mass flow between mass_flow_min and mass_flow_max).
      - For each speed curve, it computes the predicted pressure ratio at the given mass flow.
      - It then linearly interpolates in the rpm dimension to find the speed that would yield pr_fixed.
      - Finally, it uses efficiency data to estimate the efficiency at (mass_flow, pr_fixed).
    """
    if mass_flow < mass_flow_min or mass_flow > mass_flow_max:
        raise ValueError(f"Mass flow {mass_flow} is outside the operating envelope [{mass_flow_min}, {mass_flow_max}].")
    
    # Gather the pressure ratio at the given mass flow for each available speed curve.
    rpms = []
    pr_values = []
    for rpm in sorted(speed_interp_funcs.keys()):
        pr_val = speed_interp_funcs[rpm](mass_flow)
        rpms.append(rpm)
        pr_values.append(pr_val)
    rpms = np.array(rpms)
    pr_values = np.array(pr_values)
    
    # Find the two speed curves that bracket the fixed pressure ratio.
    if pr_fixed < pr_values.min() or pr_fixed > pr_values.max():
        # If the fixed pressure ratio is outside the range, extrapolate using the two extreme points.
        if pr_fixed < pr_values.min():
            i_low, i_high = 0, 1
        else:
            i_low, i_high = -2, -1
    else:
        i_low = None
        for i in range(len(pr_values) - 1):
            # Assuming the pressure ratio increases with rpm:
            if pr_values[i] <= pr_fixed <= pr_values[i+1]:
                i_low = i
                break
        if i_low is None:
            # If not found, try the reverse (if the curves are descending).
            for i in range(len(pr_values) - 1):
                if pr_values[i] >= pr_fixed >= pr_values[i+1]:
                    i_low = i
                    break
        if i_low is None:
            raise ValueError("Cannot bracket the fixed pressure ratio with available data.")
        i_high = i_low + 1
    
    rpm_low = rpms[i_low]
    rpm_high = rpms[i_high]
    pr_low = pr_values[i_low]
    pr_high = pr_values[i_high]
    
    # Linear interpolation to compute the required rpm.
    if pr_high == pr_low:
        rpm_required = rpm_low
    else:
        rpm_required = rpm_low + (pr_fixed - pr_low) * (rpm_high - rpm_low) / (pr_high - pr_low)
    
    # Interpolate efficiency at (mass_flow, pr_fixed) using the efficiency data.
    query_point = np.array([mass_flow, pr_fixed])
    predicted_efficiency = griddata(eff_points, eff_values, query_point, method='linear')
    if np.isnan(predicted_efficiency):
        predicted_efficiency = griddata(eff_points, eff_values, query_point, method='nearest')
    predicted_efficiency = float(predicted_efficiency)
    
    return rpm_required, predicted_efficiency

# ---------------------------
# Example Usage
# ---------------------------
test_mass_flow = 4

try:
    predicted_rpm, predicted_eff = compressor_model_fixed_pr(test_mass_flow)
    print(f"For a mass flow of {test_mass_flow} and a fixed pressure ratio of {fixed_pressure_ratio}:")
    print(f"  Predicted Speed (rpm): {predicted_rpm:.1f}")
    print(f"  Predicted Efficiency: {predicted_eff:.3f}")
except ValueError as e:
    print(e)
