import math
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

# Fluid definition
fluid = 'Hydrogen'

# Constants for base conditions
R = 8.314               # Ideal gas constant (J/(mol*K))
Pb_kPa = 101.325        # Base pressure in kPa for density calc
Tb_K = 288.15           # Base temperature in K
M_H2 = 0.002016         # kg/mol for H2

# Pipeline defaults
DEFAULT_TOTAL_LENGTH_M = 500e3
DEFAULT_SECTION_LENGTH_M = 1e3
D_mm = 100
D_m = D_mm / 1000.0
epsilon = 0.0178e-3
T_K = 288.15

# -------------------------------
# Helper Functions
# -------------------------------
def calculate_property(P_kPa, T_K, prop):
    try:
        return CP.PropsSI(prop, 'T', T_K, 'P', P_kPa * 1000, fluid)
    except ValueError:
        if prop == 'Z':
            return 1.0
        if prop == 'V':
            return 8.76e-6
        if prop == 'D':
            return (Pb_kPa * 1000 * M_H2) / (R * Tb_K)
        raise

def calculate_velocity(P_kPa, Qb_Sm3_per_day, Z):
    return 14.737 * (Pb_kPa / Tb_K) * (Z * T_K / P_kPa) * (Qb_Sm3_per_day / (D_mm**2))

def calculate_pressure_drop_section(P_start_kPa, Qb_Sm3_per_day, sec_len_m):
    Z = calculate_property(P_start_kPa, T_K, 'Z')
    mu = calculate_property(P_start_kPa, T_K, 'V')
    rho = (P_start_kPa * 1000 * M_H2) / (R * T_K)
    v = calculate_velocity(P_start_kPa, Qb_Sm3_per_day, Z)
    Re = (rho * v * D_m) / mu

    f = 0.02
    for _ in range(30):
        new_f = 1.0 / (-2.0 * math.log10(epsilon/(3.7*D_m) + 2.51/(Re*math.sqrt(f))))**2
        if abs(new_f - f) < 1e-6:
            f = new_f
            break
        f = new_f

    return f * (sec_len_m / D_m) * (rho * v**2 / 2)

def calculate_total_pressure_drop(P1_bar, mass_flow_kg_per_day,
                                  total_length_m=DEFAULT_TOTAL_LENGTH_M,
                                  section_length_m=DEFAULT_SECTION_LENGTH_M):
    P1_kPa = P1_bar * 100.0
    rho_b = calculate_property(Pb_kPa, Tb_K, 'D')
    Qb_Sm3_day = mass_flow_kg_per_day / rho_b

    num_sections = int(total_length_m / section_length_m)
    P_current_kPa = P1_kPa
    total_delta_P = 0.0

    for _ in range(num_sections):
        dp = calculate_pressure_drop_section(P_current_kPa, Qb_Sm3_day, section_length_m)
        total_delta_P += dp
        P_current_kPa -= dp / 1000.0
        if P_current_kPa <= 0:
            raise RuntimeError("Pressure dropped below zero in section iteration")

    return total_delta_P, P_current_kPa

def solve_P2_for_flow(P1_bar, mass_flow_kg_per_day,
                      total_length_m=DEFAULT_TOTAL_LENGTH_M,
                      section_length_m=DEFAULT_SECTION_LENGTH_M):
    _, P_final_kPa = calculate_total_pressure_drop(
        P1_bar, mass_flow_kg_per_day,
        total_length_m, section_length_m
    )
    rho_b = calculate_property(Pb_kPa, Tb_K, 'D')
    Qb_Sm3_day = mass_flow_kg_per_day / rho_b
    return P_final_kPa/100.0, Qb_Sm3_day

# -------------------------------
# Main With Plotting
# -------------------------------
def main():
    pressures = [70, 80, 90, 100]  # bar
    base_mass_flow = 10000     # kg/day
    flow_factors = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]

    plt.figure(figsize=(10, 6))

    for P1 in pressures:
        flow_rates = []
        outlet_pressures = []
        for factor in flow_factors:
            test_flow = base_mass_flow * factor
            try:
                P2_bar, _ = solve_P2_for_flow(P1, test_flow,
                                              total_length_m=DEFAULT_TOTAL_LENGTH_M,
                                              section_length_m=DEFAULT_SECTION_LENGTH_M)
                flow_rates.append(test_flow)
                outlet_pressures.append(P2_bar)
            except:
                continue

        plt.plot(flow_rates, outlet_pressures, marker='o', label=f"Inlet P1 = {P1} bar")

    plt.xlabel("Mass Flow Rate (kg/day)")
    plt.ylabel("Outlet Pressure P2 (bar)")
    plt.title("Outlet Pressure vs Mass Flow Rate")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
