
LoadCSV = "industrial_load.csv" 

# Parameters
Psuc        = 20.0            # bar
Pdis        = 70.0            # bar
x           = 2.23     
N           = 2
PR          = 2.234 
Tsuc        = 293.15          # (K)
R           = 4124            # (J/kg·K)
K           = 1.4             
RPMRef      = 3600          
SM3tokg     = 0.0834          
Moteff      = 0.95
eff_isen    = 0.75             
VFD_Eff     = 0.985

# Pipeline geometry & hydraulics
Qsm3         = 51_255_602.23  # Standard m3/day
MassFlow     = 4_278_788   # kg/day      
DesignMassFlow   = 4   # kg/s
DiaInc       = 36.0       # inch
M2K          = 1.60934    
SegLen       = 500.0      # (km)
TotLength    = 1500.0     # (km)
StorLoc      = TotLength/2 #km
Booster2Str  = StorLoc % SegLen # km
DiaMM        = 36 * 25.4     # mm
Roughness    = 0.0178e-3     # m
H1           = 0.0            # m
H2           = 100.0    # m
TFLOW        = 288.15   # K
MW_H2        = 0.002    # kg/mol
Velocity     = 35       # m/s 

# Compressors
InStaUnits    = 13              
EnStaUnits    = 13               
MaxPower      = 16000       

# Economic constants 
CEPCI2007      = 525.4          
CEPCI2019      = 619.2          
CEPCI2024      = 800            
ExRate         = 0.924          
DisRate        = 0.08           
CompLife       = 15            
PipeLife       = 50             
ElecTerrif     = 0.1867         
Labrate        = 33.9         
StrCap         = 100e6

# Fixed O&M breakdown pipeline
OMR     = 0.04       
InTax  = 0.021      

# Compressor cost correlation from HDSAM
HDSAMCostUSD  = 1962.2       
ScalingExp    = 0.8225             
VFDCostUSD    = 79.93        

# Installation & indirect factors
InstFactor = 2.0          
IndirectFraction   = 0.40         

# Labor scaling
BASE_ANNUAL_HOURS = 288           
LABOR_EXPONENT    = 0.25           
INDIRECT_LABOR_FACTOR = 0.5      

# Availability & calendar
        
DAYS_PER_YEAR    = 365             
TIMESTEP_HOURS   = 1.0             
