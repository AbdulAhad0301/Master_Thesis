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
    from scipy.interpolate import LinearNDInterpolator, interp1d, griddata
except ImportError:
    install("scipy")
    from scipy.interpolate import LinearNDInterpolator, interp1d, griddata

try:
    import matplotlib.pyplot as plt
except ImportError:
    install("matplotlib")
    import matplotlib.pyplot as plt

import CoolProp.CoolProp as CP
import math

# ----------------------------
# Constants and Design Parameters
# ----------------------------
P_suc = 20e5         # Suction pressure (Pa) = 20 bar
P_disc = 70e5       # Discharge pressure (Pa) = 100 bar
# Target pressure ratio (P_disc/P_suc) = 5
target_pressure_ratio = P_disc / P_suc  
T_suc = 293.15       # Suction temperature (K)
R_const = 8.314      # Gas constant (J/mol·K)
k = 1.4              # Adiabatic index
eta_compressor = 0.75  # Design compressor efficiency (fraction)
fluid = 'Hydrogen'     # Fluid type

# Design parameters:
# In this revised approach, we consider the design mass flow directly in kg/s.
design_mass_flow = 4.0   # Design mass flow (kg/s) used in the efficiency map

# For the fixed-speed case, we now also use mass flow in kg/s directly.
# (In the previous version, Q_design was volumetric and needed density conversion.)

# Design speed for fixed-speed baseline (for reference only)
N_design = 4000      # Design speed (RPM) for fixed speed

# For normalization in the maps:
reference_rpm = 4000    # Set as 100% speed (from your speed map data)

# ----------------------------
# Build Efficiency Map Interpolator (Real Data)
# ----------------------------
# Load efficiency map data from CSV files.
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

# Normalize efficiency map data:
# Here, the 'mass_flow' in the CSVs is assumed in kg/s.
eff_df['mass_flow_%'] = (eff_df['mass_flow'] / design_mass_flow) * 100
eff_df['efficiency_%'] = eff_df['efficiency'] * 100

# Build the efficiency interpolator:
eff_interp = LinearNDInterpolator(
    eff_df[['mass_flow_%', 'pressure_ratio']].values,
    eff_df['efficiency_%'].values
)

# ----------------------------
# Load Speed Line Data (for Pressure Ratio Map)
# ----------------------------
speed_files = [
    ("SpeedLine1.csv", 1500),
    ("SpeedLine2.csv", 2000),
    ("SpeedLine3.csv", 2500),
    ("SpeedLine4.csv", 3000),
    ("SpeedLine5.csv", 3500),
    ("SpeedLine6.csv", 4000),  # Reference = 100%
    ("SpeedLine7.csv", 4500)
]

speed_df_list = []
for file, rpm in speed_files:
    df = pd.read_csv(file)
    df.columns = ['mass_flow', 'pressure_ratio']
    df['rpm'] = rpm
    speed_df_list.append(df)

speed_df = pd.concat(speed_df_list, ignore_index=True)

# Normalize speed line data: Assume that mass flow in these CSVs is also in kg/s.
speed_df['mass_flow_%'] = (speed_df['mass_flow'] / design_mass_flow) * 100
speed_df['speed_%'] = (speed_df['rpm'] / reference_rpm) * 100

# Build an interpolator for the pressure ratio from the speed map:
pr_interp = LinearNDInterpolator(
    speed_df[['mass_flow_%', 'speed_%']].values,
    speed_df['pressure_ratio'].values
)

def calculate_z_and_rho(P, T):
    """Calculate compressibility factor and density for the given pressure and temperature."""
    Z = CP.PropsSI('Z', 'T', T, 'P', P, 'Hydrogen')
    rho = CP.PropsSI('D', 'T', T, 'P', P, 'Hydrogen')
    return Z, rho

def compressor_work(P_suc, P_disc, T_suc, m_dot, R, z, eta_compressor, k):
    return (k / (k - 1)) * (z / eta_compressor) * T_suc * m_dot * R * ((P_disc / P_suc)**((k - 1) / k) - 1)
# ----------------------------
# Compressor Model Function Using Real Efficiency Map & Speed Map
# ----------------------------
def compressor_model(mass_flow, pressure_ratio_query):
    """
    Given a mass_flow (in kg/s) and a target pressure ratio (e.g., target_pressure_ratio),
    this function:
      - Normalizes the mass flow relative to design_mass_flow.
      - Searches over a range of speed percentages (80%-120% of reference) to find the speed
        that produces a predicted pressure ratio closest to pressure_ratio_query.
      - Retrieves efficiency from the efficiency map interpolator at the normalized mass flow and the target pressure ratio.
      - Enforces a minimum efficiency of 50%.
    
    Returns:
      best_speed_pct: required speed as a percentage relative to reference_rpm.
      eff: predicted efficiency (in fraction, e.g., 0.65 means 65%).
    """
    # Normalize mass flow
    mf_pct = (mass_flow / design_mass_flow) * 100
    speeds_pct = np.linspace(80, 120, 1000)  # search range
    predicted_prs = [pr_interp(mf_pct, s) for s in speeds_pct]
    predicted_prs = np.array(predicted_prs)
    diff = np.abs(predicted_prs - pressure_ratio_query)
    
    if np.all(np.isnan(diff)):
        return None, None
    best_idx = np.nanargmin(diff)
    best_speed_pct = speeds_pct[best_idx]
    
    # Query the efficiency map at the given (normalized mass flow, target pressure ratio)
    eff = eff_interp(mf_pct, pressure_ratio_query)
    if eff is None or np.isnan(eff):
        eff = eta_compressor * 100  # fallback in %, using design efficiency
    if eff < 50:
        eff = 50  # enforce minimum efficiency of 50%
    
    return best_speed_pct, eff / 100.0  # return efficiency as fraction

# ----------------------------
# Energy Saving Comparison Simulation
# ----------------------------
# New load profile: use fractions of design mass flow (4 kg/s) directly
# This will ensure that actual mass flow remains within +/- 20% of design (3.2 to 4.8 kg/s)
time = np.arange(0, 24, 1)  # Hours 0 to 23
load_profile = np.array([
    1.0, 0.95, 1.05, 1.1, 0.9, 1.0, 0.85, 1.15,
    1.0, 0.9, 1.05, 1.0, 0.95, 1.05, 1.0, 0.9,
    1.1, 1.0, 0.95, 1.0, 1.05, 0.9, 1.0, 1.0
])
# Now Q_demands represents mass flow in kg/s directly:
Q_demands = load_profile * design_mass_flow  # kg/s

# Lists to store simulation outputs (for fixed and variable speed)
fixed_speed_power = []
variable_speed_power = []
speeds = []          # Compressor speeds (RPM) for each hour (variable-speed case)
efficiencies = []    # Compressor efficiencies (fraction) for each hour

# ----------------------------
# Fixed-Speed Compressor Calculation (Fixed at design speed, N_design)
# ----------------------------
for m_dot in Q_demands:
    # For fixed-speed, we assume the compressor always operates at its design mass flow:
    # Thus, we use the design mass flow, not the current mass flow, to compute power consumption
    Z_suc, rho_suc = calculate_z_and_rho(P_suc, T_suc)
    # m_dot_fixed remains constant (design condition)
    m_dot_fixed = 4.8 
    W_fixed = compressor_work(P_suc, P_disc, T_suc, m_dot_fixed, R_const, Z_suc, eta_compressor, k)
    fixed_speed_power.append(W_fixed)
total_fixed_speed_energy = sum(fixed_speed_power)

# ----------------------------
# Variable-Speed Compressor Calculation using Real Efficiency Map
# ----------------------------
for m_dot in Q_demands:
    # Calculate the required speed using the compressor model
    speed_pct, eta_pct = compressor_model(m_dot, target_pressure_ratio)
    if speed_pct is None:
        # fallback to design values if interpolation fails
        speed_pct = 100  # 100% relative to reference_rpm
        eta_pct = eta_compressor
    # Convert the speed percentage back to actual RPM (using reference_rpm)
    N_variable = (speed_pct / 100) * reference_rpm
    # But also, ensure variable speed is capped at 4500 rpm
    N_variable = min(N_variable, 4500)
    speeds.append(N_variable)
    
    # In variable-speed mode, use the actual (current) mass flow (m_dot, kg/s)
    Z_suc, rho_suc = calculate_z_and_rho(P_suc, T_suc)
    m_dot_current = m_dot   # Already in kg/s; no need to multiply by density in this model since Q_demands are mass flows
    
    # Use the efficiency (eta_pct) from the map, already as a fraction
    eff_var = eta_pct  
    efficiencies.append(eff_var)
    
    # Calculate compressor work for the variable-speed operating point
    W_var = compressor_work(P_suc, P_disc, T_suc, m_dot_current, R_const, Z_suc, eff_var, k)
    variable_speed_power.append(W_var)
total_variable_speed_energy = sum(variable_speed_power)

# ----------------------------
# Print Energy Consumption Results
# ----------------------------
print(f"Total energy consumption (Fixed Speed): {total_fixed_speed_energy / 1e6:.2f} MW")
print(f"Total energy consumption (Variable Speed): {total_variable_speed_energy / 1e6:.2f} MW")
energy_savings = ((total_fixed_speed_energy - total_variable_speed_energy) / total_fixed_speed_energy) * 100
print(f"Energy Savings: {energy_savings:.2f}%")

# ----------------------------
# Plot 1: Power and Demand vs. Time (Dual Y-Axis)
# ----------------------------
fig, ax1 = plt.subplots(figsize=(12, 8))
color_power = 'tab:blue'
ax1.plot(time, fixed_speed_power, label="Fixed Speed Power (W)", color=color_power, marker='o', linestyle='-')
ax1.plot(time, variable_speed_power, label="Variable Speed Power (W)", color='tab:orange', marker='o', linestyle='-')
ax1.set_xlabel("Time (h)")
ax1.set_ylabel("Power (W)", color=color_power)
ax1.tick_params(axis='y', labelcolor=color_power)
ax1.set_xticks(time)
ax1.grid(True)

ax2 = ax1.twinx()
color_demand = 'tab:green'
ax2.bar(time, Q_demands, width=0.8, alpha=0.3, color=color_demand, label='Demand (kg/s)')
ax2.set_ylabel("Mass Flow (kg/s)", color=color_demand)
ax2.tick_params(axis='y', labelcolor=color_demand)
ax2.set_ylim(0, max(Q_demands)*1.2)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
plt.title("Hourly Demand and Compressor Power Consumption")
plt.show()

# ----------------------------
# Plot 2: Compressor Speed and Efficiency vs. Time (Dual Y-Axis)
# ----------------------------
# ----------------------------
# Plot 2: Compressor Speed and Efficiency vs. Time (Dual Y-Axis) with Annotations
# ----------------------------
fig, ax1 = plt.subplots(figsize=(12, 8))

# Plot compressor speed (left y-axis)
color_speed = 'tab:purple'
ax1.plot(time, speeds, label="Speed (RPM)", color=color_speed, marker='o', linestyle='-')
ax1.set_xlabel("Time (h)")
ax1.set_ylabel("Speed (RPM)", color=color_speed)
ax1.tick_params(axis='y', labelcolor=color_speed)
ax1.set_xticks(time)
ax1.grid(True)

# Plot efficiency (right y-axis)
color_eff = 'tab:red'
efficiencies_percent = np.array(efficiencies) * 100
ax2 = ax1.twinx()
ax2.plot(time, efficiencies_percent, label="Efficiency (%)", color=color_eff, marker='o', linestyle='-')
ax2.set_ylabel("Efficiency (%)", color=color_eff)
ax2.tick_params(axis='y', labelcolor=color_eff)
ax2.set_ylim(0, 110)

# Annotate each efficiency point with its value
for t, eff_val in zip(time, efficiencies_percent):
    ax2.text(t, eff_val + 2, f"{eff_val:.1f}%", ha='center', va='bottom', color=color_eff, fontsize=9)

# Combine legends from both axes
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.title("Hourly Compressor Speed and Efficiency")
plt.show()
