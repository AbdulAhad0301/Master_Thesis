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

import CoolProp.CoolProp as CP
import math

# ----------------------------
# Constants and Design Parameters
# ----------------------------
P_suc = 20e5         # Suction pressure (Pa)
P_disc = 70e5       # Discharge pressure (Pa)
T_suc = 293.15       # Suction temperature (K)
R_const = 8.314      # Gas constant (J/mol·K)
k = 1.4              # Adiabatic index
eta_compressor = 0.75  # Design compressor efficiency (fraction)
fluid = 'Hydrogen'     # Fluid type

Q_design = 500       # Design volumetric flow rate (m^3/s)
N_design = 4000      # Design speed (RPM) for the fixed-speed case

# For normalization in the maps:
reference_rpm = 4000       # Set as 100% speed reference (from the map)
# Instead of a fixed value for design mass flow, we compute it using the inlet density.
# Calculate design density at inlet conditions using T_suc, P_suc:
def calculate_z_and_rho(P, T):
    Z = CP.PropsSI('Z', 'T', T, 'P', P, fluid)
    rho = CP.PropsSI('D', 'T', T, 'P', P, fluid)
    return Z, rho

Z_design, rho_design = calculate_z_and_rho(P_suc, T_suc)
# design_mass_flow now in kg/s:
design_mass_flow = Q_design * rho_design

print(f"Design density (rho_design): {rho_design:.3f} kg/m³")
print(f"Design mass flow (design_mass_flow): {design_mass_flow:.3f} kg/s")

# ----------------------------
# Build Efficiency Map Interpolator (Real Data)
# ----------------------------
# Load Efficiency Data from CSV files. The efficiency values in the file names indicate performance.
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

# Normalize efficiency map data using design_mass_flow:
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
    ("SpeedLine6.csv", 4000),  # Reference: 100%
    ("SpeedLine7.csv", 4500)
]

speed_df_list = []
for file, rpm in speed_files:
    df = pd.read_csv(file)
    df.columns = ['mass_flow', 'pressure_ratio']
    df['rpm'] = rpm
    speed_df_list.append(df)

speed_df = pd.concat(speed_df_list, ignore_index=True)
# Normalize speed data using design_mass_flow and reference_rpm:
speed_df['mass_flow_%'] = (speed_df['mass_flow'] / design_mass_flow) * 100
speed_df['speed_%'] = (speed_df['rpm'] / reference_rpm) * 100

# Build interpolator for pressure ratio:
pr_interp = LinearNDInterpolator(
    speed_df[['mass_flow_%', 'speed_%']].values,
    speed_df['pressure_ratio'].values
)

# ----------------------------
# Compressor Model Function Using Real Efficiency Map and Speed Lines
# ----------------------------
def compressor_model(mass_flow, pressure_ratio_query):
    """
    mass_flow: actual mass flow (kg/s)
    pressure_ratio_query: target pressure ratio (from system requirements)
    
    The function:
      - Normalizes the mass flow (kg/s) to a percentage based on design_mass_flow.
      - Searches over a range of speed percentages (80%-120% of design speed) to find the speed 
        that produces a compressor map pressure ratio close to pressure_ratio_query.
      - Retrieves efficiency from the real efficiency map at the normalized mass flow and a reference pressure ratio (target_pressure_ratio).
      - Enforces a minimum efficiency of 50%.
      
    Returns:
      best_speed_pct: required speed in % relative to reference_rpm (4000 rpm).
      eff: interpolated efficiency in %.
    """
    mf_pct = (mass_flow / design_mass_flow) * 100
    speeds_pct = np.linspace(80, 120, 1000)  # search range (80%-120%)
    predicted_prs = [pr_interp(mf_pct, s) for s in speeds_pct]
    predicted_prs = np.array(predicted_prs)
    diff = np.abs(predicted_prs - pressure_ratio_query)
    
    if np.all(np.isnan(diff)):
        return None, None
    
    best_idx = np.nanargmin(diff)
    best_speed_pct = speeds_pct[best_idx]
    
    # For efficiency, use a reference pressure ratio, here set to 2.1 as design reference.
    target_pressure_ratio = 2.1  
    eff = eff_interp(mf_pct, target_pressure_ratio)
    if (eff is None) or np.isnan(eff):
        eff = eta_compressor * 100  # fallback to design efficiency (in %)
    # Enforce a minimum efficiency of 50%
    if eff < 50:
        eff = 50
    return best_speed_pct, eff

# ----------------------------
# Simulation: Hourly Load Profile and Variable-Speed Compressor Modeling
# ----------------------------
time = np.arange(0, 24, 1)  # Hours 0 to 23

# Example load profile (fraction of design volumetric flow, Q_design, in m^3/s)
load_profile = np.array([0.9, 0.8, 0.7, 0.6, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                          1.0, 0.8, 0.7, 0.6, 0.7, 0.8, 0.9, 0.8, 0.7, 0.6,
                          0.95, 0.88, 0.76, 0.65])
# Q_demands in m^3/s from volumetric load profile:
Q_demands = load_profile * Q_design

# Precompute design mass flow (from design volumetric flow Q_design) for normalization in variable-speed loop:
# Note: using design conditions P_suc and T_suc
_, rho_design = calculate_z_and_rho(P_suc, T_suc)
m_dot_design = Q_design * rho_design  # in kg/s; this replaces our earlier design_mass_flow value
print(f"Updated design mass flow for efficiency map normalization: {m_dot_design:.3f} kg/s")

# Lists to store outputs (fixed and variable speed cases)
fixed_speed_power = []
variable_speed_power = []
speeds = []          # Compressor speeds (RPM) for each hour
efficiencies = []    # Efficiency (fraction) for each hour

# ----------------------------
# Fixed-Speed Compressor Calculation (Assumed fixed at design speed, N_design)
# ----------------------------
for Q in Q_demands:
    Z_suc, rho_suc = calculate_z_and_rho(P_suc, T_suc)
    m_dot_des_fixed = Q_design * rho_suc  # mass flow (kg/s) at design volumetric flow under current density
    W_design = compressor_work(P_suc, P_disc, T_suc, m_dot_des_fixed, R_const, Z_suc, eta_compressor, k)
    fixed_speed_power.append(W_design)
total_fixed_speed_energy = sum(fixed_speed_power)

# ----------------------------
# Variable-Speed Compressor Calculation with Real Efficiency Map and Relaxed Speed Clipping
# ----------------------------
for Q in Q_demands:
    # Calculate ideal speed from flow demand (based on volumetric flow scaling)
    N = (Q / Q_design) * N_design  
    # Relax speed clipping: allow speeds up to 4500 rpm
    N = min(N, 4500)
    speeds.append(N)
    
    Z_suc, rho_suc = calculate_z_and_rho(P_suc, T_suc)
    m_dot_current = Q * rho_suc  # current mass flow (kg/s)
    
    # Use actual mass flow (kg/s) for efficiency map normalization:
    # Compute relative flow percentage using mass flow values:
    relative_flow_percent = (m_dot_current / m_dot_design) * 100.0
    relative_speed_percent = (N / N_design) * 100.0
    
    # Get required speed (as % of reference) and efficiency from the real map.
    speed_pct, eta_pct = compressor_model(m_dot_current, 2.1)  # Use target PR = 2.1
    if speed_pct is None:
        eta_pct = eta_compressor * 100  # fallback if map fails
    eff_var = eta_pct / 100.0  # Convert efficiency from % to fraction
    efficiencies.append(eff_var)
    
    # Calculate compressor work at variable-speed conditions
    W_actual = compressor_work(P_suc, P_disc, T_suc, m_dot_current, R_const, Z_suc, eff_var, k)
    variable_speed_power.append(W_actual)
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
ax2.bar(time, Q_demands, width=0.8, alpha=0.3, color=color_demand, label='Demand (m³/s)')
ax2.set_ylabel("Demand (m³/s)", color=color_demand)
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
fig, ax1 = plt.subplots(figsize=(12, 8))
color_speed = 'tab:purple'
ax1.plot(time, speeds, label="Speed (RPM)", color=color_speed, marker='o', linestyle='-')
ax1.set_xlabel("Time (h)")
ax1.set_ylabel("Speed (RPM)", color=color_speed)
ax1.tick_params(axis='y', labelcolor=color_speed)
ax1.set_xticks(time)
ax1.grid(True)

ax2 = ax1.twinx()
color_eff = 'tab:red'
efficiencies_percent = np.array(efficiencies) * 100
ax2.plot(time, efficiencies_percent, label="Efficiency (%)", color=color_eff, marker='o', linestyle='-')
ax2.set_ylabel("Efficiency (%)", color=color_eff)
ax2.tick_params(axis='y', labelcolor=color_eff)
ax2.set_ylim(0, 110)

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
plt.title("Hourly Compressor Speed and Efficiency")
plt.show()
