import math
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
    INDIRECT_LABOR_FACTOR

)

# --- Economic Helpers ---
def convert_USD2007_to_EUR2024(usd2007: float) -> float:
    """Convert 2007 USD to 2024 EUR."""
    return usd2007 * (CEPCI_2024 / CEPCI_2007) / EXCHANGE_RATE

def convert_USD2019_to_EUR2024(usd2019: float) -> float:
    """Convert 2019 USD to 2024 EUR."""
    return usd2019 * (CEPCI_2024 / CEPCI_2019) / EXCHANGE_RATE

def capital_recovery_factor(r: float, n: int) -> float:
    """Capital Recovery Factor for rate r over n years."""
    return r * (1 + r)**n / ((1 + r)**n - 1)

def storage_cost_usd2019(volume_m3: float) -> float:
    """Salt‐cavern storage cost in 2019 USD for given working volume."""
    norm = volume_m3 / 9e6
    cavern = (2193.5 * norm + 7_000_000) * norm
    misc   = (-0.0191 * norm + 282_603) * norm
    return cavern + misc

def convert_USD_to_EUR(usd: float) -> float:
    return usd / EXCHANGE_RATE

# --- Pipeline Economics Calculator ---
def calculate_pipeline_economics(
    compressors_per_station: int,
    enroute_stations: int,
    mass_flow_kg_per_day: float,
    annual_elec_cost_fixed: float,
    annual_elec_cost_var: float,
    storage_volume_m3: float,
) -> dict:
    """
    Returns 'fixed_speed' and 'variable_speed' dicts with:
      - UC, TIC, TCI, Annualized CapEx
      - Electricity (user-supplied), labor, fixed O&M
      - capex/kg, opex/kg, energy/kg, LCOH (EUR/kg)
    """
    CRF = capital_recovery_factor(DISCOUNT_RATE, COMP_LIFE_YRS)
    UC_comp_unit = convert_USD2007_to_EUR2024(COST_FROM_HDSAM_USD) * (RATED_POWER_KW ** SCALING_EXPONENT)

    # Compressor counts
    main_cnt = compressors_per_station
    enr_cnt  = compressors_per_station * enroute_stations

    # Annual throughput (kg H2/yr)
    throughput = AVAILABILITY * mass_flow_kg_per_day * DAYS_PER_YEAR

    # Labor (EUR/yr)
    ann_hours = BASE_ANNUAL_HOURS * ((mass_flow_kg_per_day / 1e5) ** 0.25)
    direct_labor_per_station = ann_hours * LABOR_RATE
    indirect_labor_per_station =  direct_labor_per_station * INDIRECT_LABOR_FACTOR
    labor_per_station = direct_labor_per_station + indirect_labor_per_station
    labor_total = labor_per_station * (1 + enroute_stations)

    # --- Fixed-Speed (with storage) ---
    # Bare-module CapEx
    cap_comp_main = main_cnt * UC_comp_unit
    cap_comp_enr  = enr_cnt  * UC_comp_unit
    stor_usd2019  = storage_cost_usd2019(storage_volume_m3)
    cap_stor_eur  = convert_USD2019_to_EUR2024(stor_usd2019)
    UC_fixed      = cap_comp_main + cap_comp_enr + cap_stor_eur

    # Installed & annualized CapEx
    TIC_fixed    = UC_fixed * INSTALLATION_FACTOR
    TCI_fixed    = TIC_fixed * (1 + INDIRECT_FRACTION)
    AnnCap_fixed = TCI_fixed * CRF

    # Fixed O&M
    OM_cap_fixed = OM_REPAIRS_RATE * TIC_fixed
    OM_tax_fixed= INSURANCE_TAX_RATE * TCI_fixed
    OandM_fixed  = OM_cap_fixed + OM_tax_fixed

    # LCOH metrics
    capex_kg_fixed  = AnnCap_fixed / throughput
    opex_kg_fixed   = (labor_total + OandM_fixed) / throughput
    energy_kg_fixed = annual_elec_cost_fixed / throughput
    LCOH_fixed      = capex_kg_fixed + opex_kg_fixed + energy_kg_fixed

    # --- Variable-Speed (with VFD) ---

    # Bare-module CapEx
    vfd_usd          = RATED_POWER_KW * VFD_COST_PER_KW_USD   # e.g. 16 000 kW × 79.93 $/kW
    vfd_eur_per_unit = convert_USD_to_EUR(vfd_usd)
    cap_vfd_main     = vfd_eur_per_unit * main_cnt
    cap_vfd_enr      = vfd_eur_per_unit * enr_cnt
    UC_var           = cap_comp_main + cap_comp_enr + cap_vfd_main + cap_vfd_enr

    # Installed & annualized CapEx
    TIC_var    = UC_var * INSTALLATION_FACTOR
    TCI_var    = TIC_var * (1 + INDIRECT_FRACTION)
    AnnCap_var = TCI_var * CRF

    # Fixed O&M
    OM_cap_var = OM_REPAIRS_RATE    * TIC_var
    OM_tax_var = INSURANCE_TAX_RATE * TCI_var
    OandM_var  = OM_cap_var + OM_tax_var

    # LCOH metrics
    capex_kg_var  = AnnCap_var / throughput
    opex_kg_var   = (labor_total + OandM_var) / throughput
    energy_kg_var = annual_elec_cost_var / throughput
    LCOH_var      = capex_kg_var + opex_kg_var + energy_kg_var

    return {
        'fixed_speed': {
            'UC_EUR':          UC_fixed,
            'TIC_EUR':         TIC_fixed,
            'TCI_EUR':         TCI_fixed,
            'AnnualCapEx_EUR': AnnCap_fixed,
            'Electricity_EUR': annual_elec_cost_fixed,
            'Labor_EUR':       labor_total,
            'Fixed_O&M_EUR':   OandM_fixed,
            'CapEx_per_kg':    capex_kg_fixed,
            'Opex_per_kg':     opex_kg_fixed,
            'Energy_per_kg':   energy_kg_fixed,
            'LCOH_per_kg':     LCOH_fixed,
        },
        'variable_speed': {
            'UC_EUR':          UC_var,
            'TIC_EUR':         TIC_var,
            'TCI_EUR':         TCI_var,
            'AnnualCapEx_EUR': AnnCap_var,
            'Electricity_EUR': annual_elec_cost_var,
            'Labor_EUR':       labor_total,
            'Fixed_O&M_EUR':   OandM_var,
            'CapEx_per_kg':    capex_kg_var,
            'Opex_per_kg':     opex_kg_var,
            'Energy_per_kg':   energy_kg_var,
            'LCOH_per_kg':     LCOH_var,
        }
    }
if __name__ == "__main__":
    # Example parameters
    compressors_per_station = 13
    enroute_stations        = 2
    # Assume profile summed to 1,000,000 kg H₂ per year
    avg_daily_flow_kg       = 1_000_000 / DAYS_PER_YEAR
    # Example annual electricity costs (EUR/year)
    annual_elec_fixed       = 150_000.0
    annual_elec_var         =  50_000.0
    # Storage and VFD settings
    storage_volume_m3       = 500_000.0

    # Call the calculator
    results = calculate_pipeline_economics(
        compressors_per_station   = compressors_per_station,
        enroute_stations          = enroute_stations,
        mass_flow_kg_per_day      = avg_daily_flow_kg,
        annual_elec_cost_fixed    = annual_elec_fixed,
        annual_elec_cost_var      = annual_elec_var,
        storage_volume_m3         = storage_volume_m3,
    )

    # Print results
    print("\nFixed-Speed Pipeline Economics:\n" + "-"*40)
    for key, value in results['fixed_speed'].items():
        print(f"{key:22s}: {value:,.2f}")

    print("\nVariable-Speed Pipeline Economics:\n" + "-"*40)
    for key, value in results['variable_speed'].items():
        print(f"{key:22s}: {value:,.2f}")
