import math
import pandas as pd
from rich.console import Console
from rich.table import Table
from rich.progress import track
from config import (
    PSUC_BAR, PDISC_BAR,
    TFLOW_K,
    R_SPEC, K_GAS, REF_RPM,
    PER_STAGE_PR,
    ETA,
    NUMBER_OF_STAGES, MOTOR_Eff,
    DESIGN_MASS_FLOW_4kg,
    LOAD_CSV
)
from pd import calculate_total_pressure_drop
from Effmap import cal_Z_Rho, compModel, comp_power




# Simulation of Energy
def simulate_energy(m_dot_design, units_init, n_enroute_stations, section_length_km, booster_to_storage_km):
    df = pd.read_csv(LOAD_CSV)
    fracs = df["NormalizedLoad"].to_numpy()
    sum_fix_i = sum_var_i = sum_fix_e = sum_var_e = storage_kWh = sum_mass = 0.0
    max_storage_units = 0

    for frac in track(fracs, total=len(fracs), description="Simulating hours"):
    
        m_dot = frac * m_dot_design
        sum_mass += m_dot * 3600
        # pressure drop for section
        _, P2 = calculate_total_pressure_drop(PDISC_BAR, m_dot * 86400, total_length_m=section_length_km*1e3)
        Z_en, _ = cal_Z_Rho(P2*1e5, TFLOW_K)
        Z_in, _ = cal_Z_Rho(PSUC_BAR*1e5, TFLOW_K)


        # initial station energy
        needed_init = min(math.ceil(m_dot / DESIGN_MASS_FLOW_4kg), units_init)
        flow_init = needed_init * DESIGN_MASS_FLOW_4kg
        excess = max(flow_init - m_dot, 0)
        if excess > 0:
            _ , Pst = calculate_total_pressure_drop(PDISC_BAR, flow_init*86400, total_length_m=booster_to_storage_km*1e3)
            Z, _ = cal_Z_Rho(Pst*1e5, TFLOW_K)
            W_s = comp_power(Pst*1e5, 200e5, TFLOW_K, excess, R_SPEC, Z, eta=ETA, k=K_GAS, N=NUMBER_OF_STAGES)/1000/MOTOR_Eff
            storage_kWh += W_s
            max_storage_units = max(max_storage_units, math.ceil(excess / DESIGN_MASS_FLOW_4kg))

        # fixed-speed initial
        Wf_i = comp_power(PSUC_BAR*1e5, PDISC_BAR*1e5, TFLOW_K, flow_init, R_SPEC, Z_in, eta=ETA, k=K_GAS, N=NUMBER_OF_STAGES)/1000/MOTOR_Eff
        sum_fix_i += Wf_i
        # variable-speed initial
        Wv_i = 0.0
        if needed_init > 0:
            each = m_dot / needed_init
            rpm, eff = compModel(each, PER_STAGE_PR)
            Wb = comp_power(PSUC_BAR*1e5, PDISC_BAR*1e5, TFLOW_K, each, R_SPEC, Z_in, eff, k=K_GAS, N=NUMBER_OF_STAGES)/1000/MOTOR_Eff
            Wv_i = Wb * (rpm / REF_RPM)**3 * needed_init
        sum_var_i += Wv_i

        # fixed-speed enroute
        Wf_e = comp_power(P2*1e5, PDISC_BAR*1e5, TFLOW_K, DESIGN_MASS_FLOW_4kg, R_SPEC, Z_en, eta=ETA, k=K_GAS, N=NUMBER_OF_STAGES)/1000/MOTOR_Eff
        sum_fix_e += Wf_e * needed_init * n_enroute_stations
        # variable-speed enroute
        Wv_e = 0.0
        if needed_init > 0:
            each_e = m_dot / needed_init
            rpm_e, eff_e = compModel(each_e, PER_STAGE_PR)
            Wb_e = comp_power(P2*1e5, PDISC_BAR*1e5, TFLOW_K, each_e, R_SPEC,Z_en, eff_e, k=K_GAS, N=NUMBER_OF_STAGES)/1000/MOTOR_Eff
            Wv_e = Wb_e * (rpm_e / REF_RPM)**3 * needed_init * n_enroute_stations
        sum_var_e += Wv_e
    return sum_fix_i, sum_fix_e, sum_var_i, sum_var_e, storage_kWh, max_storage_units, sum_mass, Z_in,Z_en