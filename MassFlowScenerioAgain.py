import pandas as pd
from rich.console import Console
from rich.table import Table
from config import (
  ELEC_TARIFF, TEST_MASS_FLOW_KG_PER_DAY, SEGMENT_LEN_KM, BOOSTER_TO_STORAGE_KM
)

from LCOH import calculate_pipeline_economics
from FlowCalc import simulate_energy

# --- Main scenarios ---
def main():
    scenarios = [
        {"label": "40% of capacity", "flow_kg_day": int(0.4 * TEST_MASS_FLOW_KG_PER_DAY), "init": 5,  "enroute": 5},
        {"label": "60% of capacity", "flow_kg_day": int(0.6 * TEST_MASS_FLOW_KG_PER_DAY), "init": 8,  "enroute": 8},
        {"label": "80% of capacity", "flow_kg_day": int(0.8 * TEST_MASS_FLOW_KG_PER_DAY), "init": 10, "enroute": 10},
        {"label": "100% of capacity","flow_kg_day": int(1.0 * TEST_MASS_FLOW_KG_PER_DAY), "init": 13, "enroute": 13}
    ]

    console = Console()
    tbl = Table(title='Scenario Comparison (€/kg H₂)')
    for col in ['Scenario', 'Mode', 'CapEx €/kg', 'Energy Opex €/kg', 'Non-energy Opex €/kg', 'LCOH €/kg']:
        tbl.add_column(col, justify='right' if 'Scenario' not in col else 'left')

    for s in scenarios:
        m_dot_design = s['flow_kg_day'] / 86400.0
        fi, fe, vi, ve, sk, max_u, mass, Z_initial, Z_enroute = simulate_energy(
            m_dot_design, s['init'], 2, SEGMENT_LEN_KM, BOOSTER_TO_STORAGE_KM
        )
        print("Max storae units: ", max_u)
        print ("Z for initial station: ", Z_initial)
        print ("Z for enroute stations: ", Z_enroute)
        # Compute economics for both modes
        for mode, (Ei, Ee, Es) in [('Fixed', (fi, fe, sk)), ('Variable', (vi, ve, 0.0))]:
            vals = calculate_pipeline_economics(
                s['init'], s['enroute'], 2, m_dot_design,
             Ei, Ee, Es, vi, ve, ELEC_TARIFF, max_u
            )
            metrics = vals[mode.lower()]
            tbl.add_row(
                s['label'], mode,
                f"{metrics['capex_per_kg']:.2f}",
                f"{metrics['energy_opex_per_kg']:.2f}",
                f"{metrics['nonenergy_opex_per_kg']:.2f}",
                f"{metrics['lcoh_per_kg']:.2f}"
            )

    console.print(tbl)

if __name__ == '__main__':
    main()
