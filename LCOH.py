import pandas as pd
from config import (
    RATED_POWER_KW, CEPCI_2007,
    CEPCI_2019, CEPCI_2024, EXCHANGE_RATE,
    DISCOUNT_RATE, COMP_LIFE_YRS,
    BASE_ANNUAL_HOURS, LABOR_RATE,
    OM_REPAIRS_RATE, INSURANCE_TAX_RATE,
    COST_FROM_HDSAM_USD, SCALING_EXPONENT,
    INSTALLATION_FACTOR, INDIRECT_FRACTION,
    INDIRECT_LABOR_FACTOR,
    VFD_COST_PER_KW_USD,
    LOAD_CSV, STR_CAP
)

# Helper functions
def convert_USD2007_to_EUR2024(usd2007: float) -> float:
    return usd2007 * (CEPCI_2024 / CEPCI_2007) / EXCHANGE_RATE

def convert_USD2019_to_EUR2024(usd2019: float) -> float:
    return usd2019 * (CEPCI_2024 / CEPCI_2019) / EXCHANGE_RATE

def capital_recovery_factor(r: float, n: int) -> float:
    return r * (1 + r)**n / ((1 + r)**n - 1)

def storage_capex_eur(Vcav_m3: float) -> float:
    # HDSAM V5 eqns
    C_cav = ((2193.5*Vcav_m3/9e6) + 7e6)*Vcav_m3/9e6
    C_misc = (-0.0191*((Vcav_m3/1000)/(Vcav_m3/9e6))**2 + 343.88*((Vcav_m3/1000)/(Vcav_m3/9e6)) + 282603) * (Vcav_m3/9e6)
    print("Installed cost: ", C_cav)
    print("Misc cost: ", C_misc)
    usd19 = C_cav + C_misc
    return convert_USD2019_to_EUR2024(usd19)


def calculate_pipeline_economics(
    n_init, n_enroute, stations, m_dot_design,
    Ei_fix, Ee_fix, Es_fix,
    Ei_var, Ee_var, price, str_units
):
    fracs = pd.read_csv(LOAD_CSV)["NormalizedLoad"].to_numpy()
    factor = (24 * 365) / len(fracs)
    mass = fracs.sum() * m_dot_design * 3600.0 * factor

    # electricity OPEX
    elec_fix = (Ei_fix + Ee_fix + Es_fix) * price
    elec_var = (Ei_var + Ee_var) * price

    # CAPEX unit and CRF
    UC = convert_USD2007_to_EUR2024(COST_FROM_HDSAM_USD) * (RATED_POWER_KW**SCALING_EXPONENT)
    CRF = capital_recovery_factor(DISCOUNT_RATE, COMP_LIFE_YRS)

    # labor
    annh = BASE_ANNUAL_HOURS * ((mass / 1e5)**INDIRECT_LABOR_FACTOR)
    labor = annh * LABOR_RATE * (1 + INDIRECT_LABOR_FACTOR)

    def metrics(cap, oandm, lab, elec):
        total = cap + oandm + lab + elec
        return {
            'capex_per_kg': cap / mass,
            'energy_opex_per_kg': elec / mass,
            'nonenergy_opex_per_kg': (oandm + lab) / mass,
            'lcoh_per_kg': total / mass
        }

    # Fixed-speed: include storage capex & O&M
    cap_i = (n_init + str_units) * UC
    cap_e = n_enroute * stations * UC
    TIC_i = cap_i * INSTALLATION_FACTOR
    TIC_e = cap_e * INSTALLATION_FACTOR
    TCI_i = TIC_i * (1 + INDIRECT_FRACTION)
    TCI_e = TIC_e * (1 + INDIRECT_FRACTION)
    Ann_i = TCI_i * CRF
    Ann_e = TCI_e * CRF
    OM_i = OM_REPAIRS_RATE * TIC_i + INSURANCE_TAX_RATE * TCI_i
    OM_e = OM_REPAIRS_RATE * TIC_e + INSURANCE_TAX_RATE * TCI_e

    # storage capex & O&M
    TIC_st = storage_capex_eur(STR_CAP)
    TCI_st = TIC_st * (1 + INDIRECT_FRACTION)
    Ann_st = TCI_st * CRF
    OM_st = OM_REPAIRS_RATE * TIC_st + INSURANCE_TAX_RATE * TCI_st

    fixed = metrics(
        Ann_i + Ann_e + Ann_st,
        OM_i + OM_e + OM_st,
        labor * (1 + stations),
        elec_fix
    )

    # Variable-speed: add VFD cost
    vfd_u = convert_USD2007_to_EUR2024(RATED_POWER_KW * VFD_COST_PER_KW_USD)
    cap_vfd_i = n_init * vfd_u
    cap_vfd_e = n_enroute * stations * vfd_u
    cap_i_var = cap_i + cap_vfd_i
    cap_e_var = cap_e + cap_vfd_e
    TIC_vi = cap_i_var * INSTALLATION_FACTOR
    TIC_ve = cap_e_var * INSTALLATION_FACTOR
    TCI_vi = TIC_vi * (1 + INDIRECT_FRACTION)
    TCI_ve = TIC_ve * (1 + INDIRECT_FRACTION)
    Ann_vi = TCI_vi * CRF
    Ann_ve = TCI_ve * CRF
    OM_vi = OM_REPAIRS_RATE * TIC_vi + INSURANCE_TAX_RATE * TCI_vi
    OM_ve = OM_REPAIRS_RATE * TIC_ve + INSURANCE_TAX_RATE * TCI_ve

    variable = metrics(
        Ann_vi + Ann_ve,
        OM_vi + OM_ve,
        labor * (1 + stations),
        elec_var
    )

    return {'fixed': fixed, 'variable': variable}
