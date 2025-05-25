import math
import pandas as pd
from rich.console import Console
from rich.table import Table
from config import (
    RATED_POWER_KW,
    CEPCI_2007,
    CEPCI_2019,
    CEPCI_2024,
    EXCHANGE_RATE,
    DISCOUNT_RATE,
    COMP_LIFE_YRS,
    AVAILABILITY,
    DAYS_PER_YEAR,
    LABOR_RATE,
    BASE_ANNUAL_HOURS,
    OM_REPAIRS_RATE,
    INSURANCE_TAX_RATE,
    COST_FROM_HDSAM_USD,
    SCALING_EXPONENT,
    INSTALLATION_FACTOR,
    INDIRECT_FRACTION,
    VFD_COST_PER_KW_USD,
    INDIRECT_LABOR_FACTOR,
    ELEC_TARIFF,
    Q_STD,
    CONV_FACTOR
)

# --- Economic Helpers ---
def convert_USD2007_to_EUR2024(usd2007: float) -> float:
    return usd2007 * (CEPCI_2024 / CEPCI_2007) / EXCHANGE_RATE


def convert_USD_to_EUR(usd: float) -> float:
    return usd / EXCHANGE_RATE


def capital_recovery_factor(r: float, n: int) -> float:
    return r * (1 + r) ** n / ((1 + r) ** n - 1)

# --- Geological Storage Helper ---
def storage_capex_eur(Vcav_m3: float) -> float:
    """HDSAM V5 cavern + misc cost → 2024 EUR."""
    # eqn (22)
    C_cavern = (2_193.5 * Vcav_m3 / 9_000_000 + 7_000_000) * Vcav_m3 / 9_000_000
    # eqn (23)
    C_misc = (-0.0191 * Vcav_m3 / 9_000_000 + 282_603) * Vcav_m3 / 9_000_000
    C_2019USD = C_cavern + C_misc  # eqn (24)
    # CEPCI+FX → 2024 EUR
    return C_2019USD * (CEPCI_2024 / CEPCI_2019) / EXCHANGE_RATE

# --- Pipeline Economics Calculator using variable load CSV ---
def calculate_pipeline_economics(
    compressors_inlet: int,
    compressors_enroute: int,
    enroute_stations: int,
    m_dot_design: float,
    energy_inlet_fixed_kWh: float,
    energy_enroute_fixed_kWh: float,
    energy_inlet_var_kWh: float,
    energy_enroute_var_kWh: float,
    electricity_price_per_kWh: float
) -> dict:
    # Load hourly load fractions
    df = pd.read_csv("industrial_load.csv")
    fractions = df["NormalizedLoad"].to_numpy()
    hours_per_year = DAYS_PER_YEAR * 24
    factor = hours_per_year / len(fractions)
    annual_mass_kg = fractions.sum() * m_dot_design * 3600.0 * factor

    # Annual electricity costs
    annual_elec_fixed_inlet = energy_inlet_fixed_kWh * electricity_price_per_kWh
    annual_elec_fixed_enroute = energy_enroute_fixed_kWh * electricity_price_per_kWh
    annual_elec_var_inlet = energy_inlet_var_kWh * electricity_price_per_kWh
    annual_elec_var_enroute = energy_enroute_var_kWh * electricity_price_per_kWh

    # CapEx factors
    CRF = capital_recovery_factor(DISCOUNT_RATE, COMP_LIFE_YRS)
    UC_comp_unit = convert_USD2007_to_EUR2024(COST_FROM_HDSAM_USD) * (RATED_POWER_KW ** SCALING_EXPONENT)

    # Labor per station
    ann_hours = BASE_ANNUAL_HOURS * ((annual_mass_kg / 1e5) ** INDIRECT_LABOR_FACTOR)
    direct_lab = ann_hours * LABOR_RATE
    indirect_lab = direct_lab * INDIRECT_LABOR_FACTOR
    labor_per_station = direct_lab + indirect_lab

    def per_kg(cap, oandm, labor, elec):
        return {
            'capex_per_kg': cap / annual_mass_kg,
            'opex_per_kg': (labor + oandm) / annual_mass_kg,
            'energy_per_kg': elec / annual_mass_kg,
            'lcoh_per_kg': cap / annual_mass_kg + (labor + oandm) / annual_mass_kg + elec / annual_mass_kg
        }

    # --- Fixed-Speed Calculations ---
    cap_inlet = compressors_inlet * UC_comp_unit
    cap_enroute = compressors_enroute * enroute_stations * UC_comp_unit
    TIC_inlet = cap_inlet * INSTALLATION_FACTOR
    TIC_enroute = cap_enroute * INSTALLATION_FACTOR
    TCI_inlet = TIC_inlet * (1 + INDIRECT_FRACTION)
    TCI_enroute = TIC_enroute * (1 + INDIRECT_FRACTION)
    AnnCap_inlet = TCI_inlet * CRF
    AnnCap_enroute = TCI_enroute * CRF
    OM_inlet = OM_REPAIRS_RATE * TIC_inlet + INSURANCE_TAX_RATE * TCI_inlet
    OM_enroute = OM_REPAIRS_RATE * TIC_enroute + INSURANCE_TAX_RATE * TCI_enroute

    inlet_fixed = {
        'AnnualCapEx_EUR': AnnCap_inlet,
        'Electricity_EUR': annual_elec_fixed_inlet,
        'Labor_EUR': labor_per_station,
        'Fixed_O&M_EUR': OM_inlet,
        **per_kg(AnnCap_inlet, OM_inlet, labor_per_station, annual_elec_fixed_inlet)
    }
    enroute_fixed = {
        'AnnualCapEx_EUR': AnnCap_enroute,
        'Electricity_EUR': annual_elec_fixed_enroute,
        'Labor_EUR': labor_per_station * enroute_stations,
        'Fixed_O&M_EUR': OM_enroute,
        **per_kg(AnnCap_enroute, OM_enroute, labor_per_station * enroute_stations, annual_elec_fixed_enroute)
    }
    total_fixed = {k: inlet_fixed[k] + enroute_fixed[k] for k in inlet_fixed}
    total_fixed.update(per_kg(
        total_fixed['AnnualCapEx_EUR'],
        total_fixed['Fixed_O&M_EUR'],
        total_fixed['Labor_EUR'],
        total_fixed['Electricity_EUR']
    ))

    # --- Inject Geological Storage for Fixed-Speed Only ---
    Vcav = 140e6  # cavern volume in m³
    cap_storage = storage_capex_eur(Vcav)
    TIC_storage = cap_storage * INSTALLATION_FACTOR
    TCI_storage = TIC_storage * (1 + INDIRECT_FRACTION)
    AnnCap_storage = TCI_storage * CRF
    OM_storage = OM_REPAIRS_RATE * TIC_storage + INSURANCE_TAX_RATE * TCI_storage

    total_fixed['AnnualCapEx_EUR'] += AnnCap_storage
    total_fixed['Fixed_O&M_EUR'] += OM_storage
    total_fixed.update(per_kg(
        total_fixed['AnnualCapEx_EUR'],
        total_fixed['Fixed_O&M_EUR'],
        total_fixed['Labor_EUR'],
        total_fixed['Electricity_EUR']
    ))

    # --- Variable-Speed Calculations ---
    vfd_cost_unit = convert_USD_to_EUR(RATED_POWER_KW * VFD_COST_PER_KW_USD)
    cap_vfd_inlet = compressors_inlet * vfd_cost_unit
    cap_vfd_enroute = compressors_enroute * enroute_stations * vfd_cost_unit
    cap_var_inlet = cap_inlet + cap_vfd_inlet
    cap_var_enroute = cap_enroute + cap_vfd_enroute
    TIC_vi = cap_var_inlet * INSTALLATION_FACTOR
    TIC_ve = cap_var_enroute * INSTALLATION_FACTOR
    TCI_vi = TIC_vi * (1 + INDIRECT_FRACTION)
    TCI_ve = TIC_ve * (1 + INDIRECT_FRACTION)
    AnnCap_vi = TCI_vi * CRF
    AnnCap_ve = TCI_ve * CRF
    OM_vi = OM_REPAIRS_RATE * TIC_vi + INSURANCE_TAX_RATE * TCI_vi
    OM_ve = OM_REPAIRS_RATE * TIC_ve + INSURANCE_TAX_RATE * TCI_ve

    inlet_var = {
        'AnnualCapEx_EUR': AnnCap_vi,
        'Electricity_EUR': annual_elec_var_inlet,
        'Labor_EUR': labor_per_station,
        'Fixed_O&M_EUR': OM_vi,
        **per_kg(AnnCap_vi, OM_vi, labor_per_station, annual_elec_var_inlet)
    }
    enroute_var = {
        'AnnualCapEx_EUR': AnnCap_ve,
        'Electricity_EUR': annual_elec_var_enroute,
        'Labor_EUR': labor_per_station * enroute_stations,
        'Fixed_O&M_EUR': OM_ve,
        **per_kg(AnnCap_ve, OM_ve, labor_per_station * enroute_stations, annual_elec_var_enroute)
    }
    total_var = {k: inlet_var[k] + enroute_var[k] for k in inlet_var}
    total_var.update(per_kg(
        total_var['AnnualCapEx_EUR'],
        total_var['Fixed_O&M_EUR'],
        total_var['Labor_EUR'],
        total_var['Electricity_EUR']
    ))

    return {
        'fixed_speed': {'inlet': inlet_fixed, 'enroute': enroute_fixed, 'total': total_fixed},
        'variable_speed': {'inlet': inlet_var, 'enroute': enroute_var, 'total': total_var}
    }

# --- Example Usage ---
if __name__ == '__main__':
    compressors_inlet = 13
    compressors_enroute = 13
    enroute_stations = 2
    m_dot_design = Q_STD * CONV_FACTOR / 86400.0

    results = calculate_pipeline_economics(
        compressors_inlet,
        compressors_enroute,
        enroute_stations,
        m_dot_design,
        energy_inlet_fixed_kWh=1_042_899_479.51 + 22_947_423.19,
        energy_enroute_fixed_kWh=169_895_399.81,
        energy_inlet_var_kWh=624_344_557.14,
        energy_enroute_var_kWh=102_263_034.11,
        electricity_price_per_kWh=ELEC_TARIFF
    )

    console = Console()
    table = Table(title='Pipeline Economics Overview')
    table.add_column('Metric', style='cyan')
    table.add_column('Fix Inlet', justify='right')
    table.add_column('Fix Enroute', justify='right')
    table.add_column('Fix Total', justify='right')
    table.add_column('Var Inlet', justify='right')
    table.add_column('Var Enroute', justify='right')
    table.add_column('Var Total', justify='right')

    metrics = ['AnnualCapEx_EUR', 'Electricity_EUR', 'Labor_EUR', 'Fixed_O&M_EUR',
               'capex_per_kg', 'opex_per_kg', 'energy_per_kg', 'lcoh_per_kg']

    fix = results['fixed_speed']
    var = results['variable_speed']

    rows = []
    for m in metrics:
        table.add_row(
            m,
            f"{fix['inlet'][m]:,.2f}",
            f"{fix['enroute'][m]:,.2f}",
            f"{fix['total'][m]:,.2f}",
            f"{var['inlet'][m]:,.2f}",
            f"{var['enroute'][m]:,.2f}",
            f"{var['total'][m]:,.2f}"
        )
        rows.append({
            'Metric': m,
            'Fix Inlet': round(fix['inlet'][m], 2),
            'Fix Enroute': round(fix['enroute'][m], 2),
            'Fix Total': round(fix['total'][m], 2),
            'Var Inlet': round(var['inlet'][m], 2),
            'Var Enroute': round(var['enroute'][m], 2),
            'Var Total': round(var['total'][m], 2)
        })

    console.print(table)

    # Export to CSV
    df_out = pd.DataFrame(rows)
    csv_path = 'pipeline_economics_100.csv'
    df_out.to_csv(csv_path, index=False)
    print(f"→ Results exported to {csv_path}")
