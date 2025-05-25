import pandas as pd
import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d, griddata
import math
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP
from config import NUMBER_OF_STAGES

# ----------------------------
# Helper Functions
# ----------------------------
def calculate_z_and_rho(P, T):
    """Calculate compressibility factor and density for hydrogen at given pressure (Pa) and temperature (K)."""
    Z = CP.PropsSI('Z', 'T', T, 'P', P, 'Hydrogen')
    rho = CP.PropsSI('D', 'T', T, 'P', P, 'Hydrogen')
    return Z, rho

def compressor_work(P_suc, P_disc, T_suc, m_dot, R, z, eta, k, N):
    """Compute compressor work (W) based on thermodynamic relations."""
    return (k / (k - 1)) * N *  (z / eta) * T_suc * m_dot * R * ((P_disc / P_suc)**((k - 1) / (N*k)) - 1)

# ----------------------------
# Constants and Design Parameters (Absolute Units)
# ----------------------------
P_suc = 20e5         # Suction pressure (Pa) = 20 bar
P_disc = 70e5        # Discharge pressure (Pa) = 70 bar
# The target pressure ratio is defined as:
target_pressure_ratio = 2.3  # for 70 bar/20 bar = 3.5
T_suc = 293.15       # Suction temperature (K)
R_const = 4124      # Gas constant (J/kg·K)
k = 1.4              # Adiabatic index
eta_compressor = 0.75  # Design compressor efficiency (fraction)
fluid = 'Hydrogen'     # Fluid type

# Here, design_mass_flow is provided in kg/s (absolute units used by our maps).
design_mass_flow = 4.0   # kg/s (this is used only for reference; no normalization is applied)

# Design speed for fixed-speed baseline:
N_design = 4000     # RPM (for baseline fixed-speed case)
#reference_rpm = 4000  # Reference speed for the speed map (all speeds in map are in RPM)

# ----------------------------
# Build Efficiency Map Interpolator (Real Data from Subfolder)
# ----------------------------
# Assume efficiency data CSV files are stored in a subfolder named "efficiency_data"
eff_folder = "EffMap/"
eff_files = [
    "Eff0_65.csv",
    "Eff0_70.csv",
    "Eff0_72.csv",
    "Eff0_74.csv",
    "Eff0_75.csv"
]

eff_points = []  # Each row: [mass_flow (kg/s), pressure_ratio]
eff_values = []  # Efficiency values (as a fraction; e.g., 0.65)
for file in eff_files:
    path = eff_folder + file
    df = pd.read_csv(path)
    df.columns = ['mass_flow', 'pressure_ratio']
    try:
        # Extract the nominal efficiency from the filename (e.g., "Eff0_65.csv" -> 0.65)
        eff_val = float(file.replace("Eff", "").replace("_", ".").replace(".csv", ""))
    except:
        continue
    for _, row in df.iterrows():
        if not pd.isna(row['mass_flow']) and not pd.isna(row['pressure_ratio']):
            eff_points.append([row['mass_flow'], row['pressure_ratio']])
            eff_values.append(eff_val)
eff_points = np.array(eff_points)
eff_values = np.array(eff_values)

# Build the efficiency interpolator using absolute values:
eff_interp = LinearNDInterpolator(eff_points, eff_values)

# ----------------------------
# Load Speed Line Data (for Pressure Ratio Map) Using Absolute Units
# ----------------------------
speed_files = [
    ("SpeedLine1.csv", 1500),
    ("SpeedLine2.csv", 2000),
    ("SpeedLine3.csv", 2500),
    ("SpeedLine4.csv", 3000),
    ("SpeedLine5.csv", 3500),
    ("SpeedLine6.csv", 4000),  # Reference speed line
    ("SpeedLine7.csv", 4500)
]

speed_folder = "EffMap/"
speed_interp_funcs = {}
mass_flow_min = 2  # Lower bound for mass flow (kg/s)
mass_flow_max = 6  # Upper bound for mass flow (kg/s)

for filename, rpm in speed_files:
    path = speed_folder + filename
    df = pd.read_csv(path)
    df.columns = ['mass_flow', 'pressure_ratio']
    df_operating = df[(df['mass_flow'] >= mass_flow_min) & (df['mass_flow'] <= mass_flow_max)]
    if len(df_operating) >= 2:
        f_interp = interp1d(df_operating['mass_flow'], df_operating['pressure_ratio'],
                            kind='linear', fill_value="extrapolate")
        speed_interp_funcs[rpm] = f_interp

# ----------------------------
# Compressor Model Function (Without Normalization)
# ----------------------------
def compressor_model(mass_flow, pressure_ratio_query):
    """
    Given an absolute mass_flow (kg/s) and target pressure ratio,
    this function:
      - Uses the raw speed map data to compute the compressor pressure ratio at various speeds.
      - Determines, by linear interpolation, the compressor speed (in RPM) required to hit pressure_ratio_query.
      - Uses the efficiency map interpolator (eff_interp) to interpolate the efficiency at the given mass_flow and pressure_ratio_query.
      - Enforces a minimum efficiency of 50%.
    
    Returns:
      rpm_required: Required compressor speed in RPM.
      predicted_eff: Predicted efficiency (as a fraction, e.g., 0.65 for 65%).
    """
    rpms = []
    pr_list = []
    for rpm in sorted(speed_interp_funcs.keys()):
        pr_val = speed_interp_funcs[rpm](mass_flow)
        rpms.append(rpm)
        pr_list.append(pr_val)
    rpms = np.array(rpms)
    pr_list = np.array(pr_list)
    
    # Find two bracketing speeds
    if pressure_ratio_query < pr_list.min() or pressure_ratio_query > pr_list.max():
        if pressure_ratio_query < pr_list.min():
            i_low, i_high = 0, 1
        else:
            i_low, i_high = -2, -1
    else:
        i_low = None
        for i in range(len(pr_list)-1):
            if pr_list[i] <= pressure_ratio_query <= pr_list[i+1]:
                i_low = i
                break
        if i_low is None:
            for i in range(len(pr_list)-1):
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
    
    # Interpolate efficiency using griddata without any normalization:
    query_point = np.array([mass_flow, pressure_ratio_query])
    predicted_eff = griddata(eff_points, eff_values, query_point, method='linear')
    if predicted_eff is None or np.isnan(predicted_eff):
        predicted_eff = griddata(eff_points, eff_values, query_point, method='nearest')
    predicted_eff = float(predicted_eff.item())
    
    if predicted_eff < 0.5:
        predicted_eff = 0.5  # Enforce a minimum efficiency of 50%
        
    return rpm_required, predicted_eff

# ----------------------------
# Example Test of Compressor Model Function
# ----------------------------
#test_mass_flow = 3.5  # kg/s
#try:
#    predicted_rpm, predicted_eff = compressor_model(test_mass_flow, target_pressure_ratio)
#    print(f"At mass flow {test_mass_flow} kg/s and target PR = {target_pressure_ratio}:")
#    print(f"Predicted compressor speed: {predicted_rpm:.1f} rpm")
#    print(f"Predicted efficiency: {predicted_eff*100:.1f}%")
#except ValueError as e:
#    print(e)

## ----------------------------
## Energy Saving Comparison Simulation (Using Mass Flow Directly)
## ----------------------------
## Create a new load profile (fractions of design mass flow)
#time = np.arange(0, 24, 1)  # 24 hours
#load_profile = np.array([
#    0.80, 0.95, 0.97, 1, 0.9, 1.0, 0.85, 1.0,
#    0.85, 0.8, 1.0, 1.0, 0.95, 1.0, 1.0, 0.9,
#    1, 0.80, 0.95, 1.0, 1.05, 0.9, 1.0, 1.0
#])
## Since we now use mass flow directly, Q_demands is in kg/s:
#Q_demands = load_profile * design_mass_flow  # kg/s
#
#fixed_speed_power = []
#variable_speed_power = []
#speeds = []          # Variable-speed operation (RPM)
#efficiencies = []    # Efficiency (fraction)
#
## Fixed-Speed Compressor Calculation (assumed fixed at design speed, N_design)
#for m_dot in Q_demands:
#    Z_suc, rho_suc = calculate_z_and_rho(P_suc, T_suc)
#    m_dot_fixed = design_mass_flow  # Fixed operating mass flow, i.e., constant at 4.0 kg/s
#    W_fixed = compressor_work(P_suc, P_disc, T_suc, m_dot_fixed, R_const, Z_suc, eta_compressor, k, N = NUMBER_OF_STAGES)/1000
#    print(W_fixed)
#    fixed_speed_power.append(W_fixed)
#total_fixed_speed_energy = sum(fixed_speed_power)
#
# Variable-Speed Compressor Calculation using real efficiency map

#for m_dot in Q_demands:
#    # Get required speed and efficiency from compressor model
#    speed_rpm, eta_val = compressor_model(m_dot, target_pressure_ratio)
#    if speed_rpm is None:
#        speed_rpm = N_design  # fallback to design speed
#        eta_val = eta_compressor
#    # Optionally cap variable speed at 4500 rpm
#    speed_rpm = min(speed_rpm, 4500)
#    speeds.append(speed_rpm)
#
#    # Use the actual mass flow (m_dot, kg/s) to compute compressor work
#    Z_suc, rho_suc = calculate_z_and_rho(P_suc, T_suc)
#    m_dot_current = m_dot  # Already in kg/s
#    
#    # Calculate the base work using the compressor_work function
#    W_base = compressor_work(P_suc, P_disc, T_suc, m_dot_current, R_const, Z_suc, eta_val, k, N = NUMBER_OF_STAGES)
#    
#    # Apply the affinity law correction factor:
#    # If N_design is the design speed, then the adjusted power is:
#    affinity_factor = (speed_rpm / N_design) ** 3
#    #print("Affinity factor", affinity_factor)
#    W_var = W_base * affinity_factor/1000
#    
#    variable_speed_power.append(W_var)
#    efficiencies.append(eta_val)
#
#total_variable_speed_energy = sum(variable_speed_power)
#
# ----------------------------
# Results and Plots
# ----------------------------
#print(f"Total energy consumption (Fixed Speed): {total_fixed_speed_energy / 1e6:.2f} MWh")
#print(f"Total energy consumption (Variable Speed): {total_variable_speed_energy / 1e6:.2f} MWh")
#energy_savings = ((total_fixed_speed_energy - total_variable_speed_energy) / total_fixed_speed_energy) * 100
#print(f"Energy Savings: {energy_savings:.2f}%")
#
## Plot Power and Demand vs. Time
#fig, ax1 = plt.subplots(figsize=(12, 8))
#ax1.plot(time, fixed_speed_power, label="Fixed Speed consumption (Wh)", color='tab:blue', marker='o', linestyle='-')
#ax1.plot(time, variable_speed_power, label="Variable Speed consumption (Wh)", color='tab:orange', marker='o', linestyle='-')
#ax1.set_xlabel("Time (h)")
#ax1.set_ylabel("Energy consumed (kWh)")
#ax1.set_xticks(time)
#ax1.grid(True)
#ax2 = ax1.twinx()
#ax2.bar(time, Q_demands, width=0.8, alpha=0.3, color='tab:green', label='Demand (kg/s)')
#ax2.set_ylabel("Mass Flow (kg/s)")
#ax2.set_ylim(0, max(Q_demands)*1.2)
#lines1, labels1 = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
#plt.title("Hourly Demand and Compressor Energy Consumption")
#plt.show()
#
## Plot Compressor Speed and Efficiency vs. Time, with efficiency annotations.
#fig, ax1 = plt.subplots(figsize=(12, 8))
#ax1.plot(time, speeds, label="Speed (RPM)", color='tab:purple', marker='o', linestyle='-')
#ax1.set_xlabel("Time (h)")
#ax1.set_ylabel("Speed (RPM)", color='tab:purple')
#ax1.tick_params(axis='y', labelcolor='tab:purple')
#ax1.set_xticks(time)
#ax1.grid(True)
#ax2 = ax1.twinx()
#efficiencies_percent = np.array(efficiencies) * 100
#ax2.plot(time, efficiencies_percent, label="Efficiency (%)", color='tab:red', marker='o', linestyle='-')
#ax2.set_ylabel("Efficiency (%)", color='tab:red')
#ax2.tick_params(axis='y', labelcolor='tab:red')
#ax2.set_ylim(70, 80)
#for t, eff_val in zip(time, efficiencies_percent):
#    ax2.text(t, eff_val + 1 , f"{eff_val:.1f}%", ha='center', va='bottom', color='tab:red', fontsize=9)
#lines1, labels1 = ax1.get_legend_handles_labels()
#lines2, labels2 = ax2.get_legend_handles_labels()
#ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')
#plt.title("Hourly Compressor Speed and Efficiency")
#plt.show()
