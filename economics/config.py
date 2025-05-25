
# --- Load profile ---
LOAD_CSV = "load_profile.csv"  # path to normalized load profile

# --- Technical constants ---
PSUC_BAR       = 20.0            # Suction pressure (bar)
PDISC_BAR      = 70.0            # Discharge pressure (bar)
x              = 2.1             # pressure ratio per stage
NUMBER_OF_STAGES = 2
PER_STAGE_PR   = 2.234 
TSUC_K         = 293.15          # Suction temperature (K)
R_SPEC         = 4124            # Specific gas constant for H2 (J/kg·K)
R              = 8.314           # gas constant (J/mol.K)
K_GAS          = 1.4             # Heat capacity ratio
REF_RPM        = 4000            # Reference compressor speed (rpm)
CONV_FACTOR    = 0.0834          # kg H2 per standard m3
MOTOR_Eff      = 0.95
ETA            = 0.75             #isentropic eff

# Pipeline geometry & hydraulics
Q_STD                     = 51_255_602.23  # Standard m3/day
TEST_MASS_FLOW_KG_PER_DAY = 4_278_788     # kg/day      
DESIGN_MASS_FLOW_4kg      = 4              # kg/s
D_INCHES                  = 36.0           # Diameter in inches
MILES_TO_KM               = 1.60934        # Conversion factor
SEGMENT_LEN_KM            = 500.0          # Length between enroute stations (km)
TOTAL_LEN_KM              = 1500.0         # Total pipeline length (km)
STORAGE_LOC               = TOTAL_LEN_KM/2
BOOSTER_TO_STORAGE_KM     = STORAGE_LOC % SEGMENT_LEN_KM
DIAM_MM                   = 36 * 25.4      # Pipeline diameter (mm)
ROUGH_M                   = 0.0178e-3      # Pipeline roughness (m)
ELEV_H1                   = 0.0            # Inlet elevation (m)
ELEV_H2                   = 100.0          # Outlet elevation (m)
TFLOW_K                   = 288.15         # Flowing temperature (K)
MW_H2                     = 0.002          # kg/mol for H2
TARGET_VELOCITY           = 35             # Target velocity m/s 

# --- Compressor counts ---
INITIAL_UNITS  = 13              # Compressors at initial station
ENROUTE_UNITS  = 13               # Compressors at each enroute station
MAX_CAPACITY   = 16_000          # kW per compressor unit

# --- Economic constants ---
CEPCI_2007      = 525.4          # CEPCI index for 2007
CEPCI_2019      = 619.2          # CEPCI index for 2019
CEPCI_2024      = 800            # CEPCI index for 2024
EXCHANGE_RATE   = 0.924          # USD per EUR
DISCOUNT_RATE   = 0.08           # Discount rate (8%)
COMP_LIFE_YRS   = 15             # Compressor lifetime (years)
PIPE_LIFE_YRS   = 50             # Pipeline lifetime (years)
ELEC_TARIFF     = 0.1867         # Electricity tariff (EUR/kWh)
LABOR_RATE      = 33.9           # Labor rate (EUR/h)

# Fixed O&M breakdown pipeline
OM_REPAIRS_RATE     = 0.04       # 4% of TIC
INSURANCE_TAX_RATE  = 0.021      # 2.1% of TCI

# Pipeline throughput & conversions


# Compressor cost‐correlation from HDSAM
COST_FROM_HDSAM_USD = 3083.3#1962.2       # Base cost (2007 USD)
SCALING_EXPONENT    = 0.8335#0.8225       # Size exponent for cost correlation
RATED_POWER_KW      = 16000        # your standard compressor rating
VFD_COST_PER_KW_USD = 79.93        # latest vendor quote

# Installation & indirect factors
INSTALLATION_FACTOR = 2.0          # IF
INDIRECT_FRACTION   = 0.40         # 40% of TIC

# Labor scaling
BASE_ANNUAL_HOURS = 288           # Base annual labor hours
LABOR_EXPONENT    = 0.25           # Exponent for labor scaling
INDIRECT_LABOR_FACTOR = 0.5      # 50% ovehead

# Availability & calendar
AVAILABILITY     = 0.90            # System availability
DAYS_PER_YEAR    = 365             # Days per year
TIMESTEP_HOURS   = 1.0             # Hours per simulation step
