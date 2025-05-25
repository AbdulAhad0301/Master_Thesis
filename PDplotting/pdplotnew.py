import math
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt

# Fluid definition
fluid = 'Hydrogen'

# Constants for base conditions
R = 8.314
Pb_kPa = 101.325
Tb_K = 288.15
M_H2 = 0.002016

# Pipeline parameters
DEFAULT_TOTAL_LENGTH_M = 500e3
DEFAULT_SECTION_LENGTH_M = 25e3
D_mm = 100
D_m = D_mm / 1000.0
epsilon = 0.0178e-3
T_K = 288.15

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

def simulate_pressure_along_distance(P1_bar, mass_flow_kg_per_day,
                                     total_length_m=DEFAULT_TOTAL_LENGTH_M,
                                     section_length_m=DEFAULT_SECTION_LENGTH_M):
    P1_kPa = P1_bar * 100.0
    rho_b = calculate_property(Pb_kPa, Tb_K, 'D')
    Qb_Sm3_day = mass_flow_kg_per_day / rho_b

    num_sections = int(total_length_m / section_length_m)
    P_current_kPa = P1_kPa
    pressures_bar = [P_current_kPa / 100.0]
    distances_km = [0]

    for i in range(1, num_sections + 1):
        dp = calculate_pressure_drop_section(P_current_kPa, Qb_Sm3_day, section_length_m)
        P_current_kPa -= dp / 1000.0
        if P_current_kPa <= 0:
            pressures_bar.append(0.0)
            distances_km.append(i * section_length_m / 1000.0)
            break
        pressures_bar.append(P_current_kPa / 100.0)
        distances_km.append(i * section_length_m / 1000.0)

    return distances_km, pressures_bar

def main():
    inlet_pressure_bar = 70
    flow_rates = list(range(10_000, 35_001, 5_000))  # [10k, 15k, ..., 35k]

    plt.figure(figsize=(10, 6))

    for flow in flow_rates:
        distances_km, pressures_bar = simulate_pressure_along_distance(inlet_pressure_bar, flow)
        plt.plot(distances_km, pressures_bar, marker='o', label=f"{flow:,} kg/day")

    plt.xlabel("Distance Along Pipeline (km)")
    plt.ylabel("Outlet Pressure (bar)")
    plt.title("Pressure Drop Along 500 km Pipeline")
    plt.grid(True)
    plt.legend(title="Mass Flow Rate")
    plt.xlim(0, 500)
    plt.ylim(bottom=0)
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    main()
