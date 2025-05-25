import math
import CoolProp.CoolProp as CP

# Fluid definition
fluid = 'Hydrogen'

# Constants for base conditions
R = 8.314               # Ideal gas constant (J/(mol*K))
Pb_kPa = 101.325        # Base pressure in kPa for density calc
Tb_K = 288.15           # Base temperature in K
M_H2 = 0.002016         # kg/mol for H2

# Pipeline defaults (can be overridden per-call)
DEFAULT_TOTAL_LENGTH_M = 200e3       # Total segment length in m (default 500 km)
DEFAULT_SECTION_LENGTH_M = 20e3      # Section length in m (default 20 km)
D_mm = 895.3                         # Pipe diameter in mm
D_m = D_mm / 1000.0                 # Pipe diameter in m
epsilon = 0.0178e-3                # Pipe roughness in m
T_K = 288.15                       # Flow temperature in K

# -------------------------------
# Helper Functions
# -------------------------------
def calculate_property(P_kPa, T_K, prop):
    """Get property 'Z','V','D' from CoolProp or fallback."""
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
    """Eq. (4): velocity (m/s)."""
    return 14.737 * (Pb_kPa / Tb_K) * (Z * T_K / P_kPa) * (Qb_Sm3_per_day / (D_mm**2))


def calculate_pressure_drop_section(P_start_kPa, Qb_Sm3_per_day, sec_len_m):
    """Pressure drop for one section of length sec_len_m."""
    Z = calculate_property(P_start_kPa, T_K, 'Z')
    mu = calculate_property(P_start_kPa, T_K, 'V')
    rho = (P_start_kPa * 1000 * M_H2) / (R * T_K)
    v = calculate_velocity(P_start_kPa, Qb_Sm3_per_day, Z)
    Re = (rho * v * D_m) / mu

    # Colebrook-White iteration
    f = 0.02
    for _ in range(30):
        new_f = 1.0 / (-2.0 * math.log10(epsilon/(3.7*D_m) + 2.51/(Re*math.sqrt(f))))**2
        if abs(new_f - f) < 1e-6:
            f = new_f
            break
        f = new_f

    # ΔP = f*(L/D)*(ρv²/2)
    return f * (sec_len_m / D_m) * (rho * v**2 / 2)


def calculate_total_pressure_drop(P1_bar, mass_flow_kg_per_day,
                                  total_length_m=DEFAULT_TOTAL_LENGTH_M,
                                  section_length_m=DEFAULT_SECTION_LENGTH_M):

    P1_kPa = P1_bar * 100.0
    # Standard volume flow (m³/day)
    rho_b = calculate_property(Pb_kPa, Tb_K, 'D')
    Qb_Sm3_day = mass_flow_kg_per_day / rho_b

    num_sections = int(total_length_m / section_length_m)
    P_current_kPa = P1_kPa
    total_delta_P = 0.0

    for _ in range(num_sections):
        dp = calculate_pressure_drop_section(P_current_kPa, Qb_Sm3_day, section_length_m)
        total_delta_P += dp
        P_current_kPa -= dp / 1000.0  # Pa→kPa
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


def main():
    pressures = [70, 80, 90, 100]  # bar
    base_mass_flow = 4_278_788     # kg/day (original)
    flow_factors = [1.0]  # lower, original, higher

    # Header
    print(f"{'P1 (bar)':>10} | {'Flow (%)':>8} | {'Flow (kg/day)':>14} | {'P2 (bar)':>9}")
    print('-' * 55)

    for P1 in pressures:
        for factor in flow_factors:
            test_flow = base_mass_flow * factor
            try:
                # Directly use solve_P2_for_flow to get P2_bar
                P2_bar, _ = solve_P2_for_flow(P1, test_flow,
                                             total_length_m=DEFAULT_TOTAL_LENGTH_M,
                                             section_length_m=DEFAULT_SECTION_LENGTH_M)
                print(f"{P1:10.1f} | {factor*100:8.0f}% | {test_flow:14,.0f} | {P2_bar:9.2f}")
            except Exception as e:
                print(f"{P1:10.1f} | {factor*100:8.0f}% | {test_flow:14,.0f} | {'Error':>9}")

if __name__ == '__main__':
    main()

