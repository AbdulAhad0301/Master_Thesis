import subprocess
import sys

def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

# Automatically install required libraries if not already installed
try:
    import pandas as pd
except ImportError:
    install("pandas")
    import pandas as pd

try:
    import numpy as np
except ImportError:
    install("numpy")
    import numpy as np

try:
    from scipy.interpolate import LinearNDInterpolator
except ImportError:
    install("scipy")
    from scipy.interpolate import LinearNDInterpolator

try:
    import matplotlib.pyplot as plt
except ImportError:
    install("matplotlib")
    import matplotlib.pyplot as plt


# -------- Load Speed Line Data --------
# Replace with your actual file paths if needed
speed_files = [
    ("SpeedLine1.csv", 1500),
    ("SpeedLine2.csv", 2000),
    ("SpeedLine3.csv", 2500),
    ("SpeedLine4.csv", 3000),
    ("SpeedLine5.csv", 3500),
    ("SpeedLine6.csv", 4000),  # reference line = 100%
    ("SpeedLine7.csv", 4500)
]

speed_df_list = []
for file, rpm in speed_files:
    df = pd.read_csv(file)
    df.columns = ['mass_flow', 'pressure_ratio']
    df['rpm'] = rpm
    speed_df_list.append(df)

speed_df = pd.concat(speed_df_list, ignore_index=True)

# Set 4000 rpm as 100%
reference_rpm = 4000
design_mass_flow = 4.0  # design point for mass flow

# Normalize
speed_df['mass_flow_%'] = (speed_df['mass_flow'] / design_mass_flow) * 100
speed_df['speed_%'] = (speed_df['rpm'] / reference_rpm) * 100

# -------- Load Efficiency Data --------
eff_files = [
    ("Eff0_65.csv", 0.65),
    ("Eff0_70.csv", 0.70),
    ("Eff0_72.csv", 0.72),
    ("Eff0_74.csv", 0.74),
    ("Eff0_75.csv", 0.75)
]

eff_df_list = []
for file, eff in eff_files:
    df = pd.read_csv(file)
    df.columns = ['mass_flow', 'pressure_ratio']
    df['efficiency'] = eff
    eff_df_list.append(df)

eff_df = pd.concat(eff_df_list, ignore_index=True)

# Normalize efficiency data
eff_df['mass_flow_%'] = (eff_df['mass_flow'] / design_mass_flow) * 100
eff_df['efficiency_%'] = eff_df['efficiency'] * 100

# -------- Build Interpolators --------
# Pressure ratio interpolator: (mass_flow_%, speed_%) -> pressure_ratio
pr_interp = LinearNDInterpolator(
    speed_df[['mass_flow_%', 'speed_%']].values,
    speed_df['pressure_ratio'].values
)

# Efficiency interpolator: (mass_flow_%, pressure_ratio) -> efficiency
eff_interp = LinearNDInterpolator(
    eff_df[['mass_flow_%', 'pressure_ratio']].values,
    eff_df['efficiency_%'].values
)

# -------- Compressor Model Function --------
def compressor_model(mass_flow, pressure_ratio_query):
    mf_pct = (mass_flow / design_mass_flow) * 100
    speeds_pct = np.linspace(80, 120, 1000)  # search range

    # Predict pressure ratios for all test speeds
    predicted_prs = [pr_interp(mf_pct, s) for s in speeds_pct]
    predicted_prs = np.array(predicted_prs)
    diff = np.abs(predicted_prs - pressure_ratio_query)

    if np.all(np.isnan(diff)):
        return None, None

    best_idx = np.nanargmin(diff)
    best_speed_pct = speeds_pct[best_idx]

    # Estimate efficiency at the selected operating point
    eff = eff_interp(mf_pct, pressure_ratio_query)
    if eff is None or np.isnan(eff):
        eff = 0

    return best_speed_pct, eff

# -------- Example Usage --------
if __name__ == "__main__":
    test_mass_flow = 4.5
    test_pressure_ratio = 2.5

    speed_percent, efficiency = compressor_model(test_mass_flow, test_pressure_ratio)

    if speed_percent is not None:
        print(f"At mass flow = {test_mass_flow} and PR = {test_pressure_ratio}:")
        print(f"  Required Speed: {speed_percent:.2f}% of 4000 rpm ({speed_percent/100*4000:.0f} rpm)")
        print(f"  Estimated Efficiency: {efficiency:.2f}%")
    else:
        print("Could not compute values for the given input.")
