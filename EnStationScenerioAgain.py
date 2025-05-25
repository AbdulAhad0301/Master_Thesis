from rich.console import Console
from rich.table import Table
from config import TOTAL_LEN_KM, STORAGE_LOC, ELEC_TARIFF
from LCOH import calculate_pipeline_economics
from FlowCalc import simulate_energy

def main():

    DESIGN_FLOW = 4_278_788 / 86400.0
    INIT_UNITS = 13
    EN_PER_STATION = 13
    scenarios = [
        {'stations': 2, 'label': '2 stations'},
        {'stations': 3, 'label': '3 stations'},
        {'stations': 4, 'label': '4 stations'},
        {'stations': 5, 'label': '5 stations'}
    ]

    console = Console()

    # LCOH Table
    tbl = Table(title='En-route Station Sensitivity (€/kg H₂)')
    for col in ['Scenario', 'Mode', 'CapEx €/kg', 'Energy Opex €/kg', 'Non-energy Opex €/kg', 'LCOH €/kg']:
        tbl.add_column(col, justify='right' if col != 'Scenario' else 'left')

    # Raw Energy Breakdown Table
    tbl_energy = Table(title='Raw Energy Results (MWh)')
    for col in ['Scenario', 'Fix_Init', 'Fix_Enr', 'Var_Init', 'Var_Enr', 'Storage (MWh)', 'Mass (t)']:
        tbl_energy.add_column(col, justify='right' if col != 'Scenario' else 'left')

    for scen in scenarios:
        n_st = scen['stations']
        sec_len = TOTAL_LEN_KM / (n_st + 1)
        b_to_sto = STORAGE_LOC % sec_len

        # Energy simulation
        fi, fe, vi, ve, sk, max_store, mass, Z_initial, Z_enroute = simulate_energy(
            DESIGN_FLOW, INIT_UNITS, n_st, sec_len, b_to_sto
        )
        print ("Z for initial station: ", Z_initial)
        print ("Z for enroute stations: ", Z_enroute)
        # Energy breakdown table (convert kWh to MWh, mass to tonnes)
        tbl_energy.add_row(
            scen['label'],
            f"{fi/1000:.1f}",
            f"{fe/1000:.1f}",
            f"{vi/1000:.1f}",
            f"{ve/1000:.1f}",
            f"{sk/1000:.1f}",
            f"{mass/1000:.1f}"
        )

        # LCOH computation
        econ = calculate_pipeline_economics(
            INIT_UNITS, EN_PER_STATION, n_st,
            DESIGN_FLOW, fi, fe, sk, vi, ve, ELEC_TARIFF, max_store
        )
        for mode in ['fixed', 'variable']:
            vals = econ[mode]
            tbl.add_row(
                scen['label'], mode.capitalize(),
                f"{vals['capex_per_kg']:.2f}",
                f"{vals['energy_opex_per_kg']:.2f}",
                f"{vals['nonenergy_opex_per_kg']:.2f}",
                f"{vals['lcoh_per_kg']:.2f}"
            )

    console.print(tbl)
    console.print(tbl_energy)

if __name__ == '__main__':
    main()
