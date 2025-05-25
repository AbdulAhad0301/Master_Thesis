import pandas as pd
import numpy as np
from scipy.interpolate import LinearNDInterpolator, interp1d, griddata
import CoolProp.CoolProp as CP

def cal_Z_Rho(Pressure, Temperature):
    Z = CP.PropsSI('Z', 'T', Temperature, 'P', Pressure, 'Hydrogen')
    rho = CP.PropsSI('D', 'T', Temperature, 'P', Pressure, 'Hydrogen')
    return Z, rho

def comp_power(P_suc, P_disc, T_suc, mfr, R, z, eff_isen, k, N):
    Power = (k / (k - 1)) * N *  (z / eff_isen) * T_suc * mfr * R * ((P_disc / P_suc)**((k - 1) / (N*k)) - 1)
    return Power

effFolder = "EffMap/"
effFiles = [
    "Eff_65.csv",
    "Eff_70.csv",
    "Eff_72.csv",
    "Eff_74.csv",
    "Eff_75.csv"
]

effPoints = []  
effValues = [] 
for file in effFiles:
    path = effFolder + file
    df = pd.read_csv(path)
    df.columns = ['mass_flow', 'pressure_ratio']
    try:
        effVal = float(file.replace("Eff", "").replace("_", "0.").replace(".csv", ""))
    except:
        continue
    for _, row in df.iterrows():
        if not pd.isna(row['mass_flow']) and not pd.isna(row['pressure_ratio']):
            effPoints.append([row['mass_flow'], row['pressure_ratio']])
            effValues.append(effVal)
effPoints = np.array(effPoints)
effValues = np.array(effValues)

effInterp = LinearNDInterpolator(effPoints, effValues)

speedFiles = [
    ("SpeedLine1.csv", 1500),
    ("SpeedLine2.csv", 2000),
    ("SpeedLine3.csv", 2500),
    ("SpeedLine4.csv", 3000),
    ("SpeedLine5.csv", 3500),
    ("SpeedLine6.csv", 4000),  
    ("SpeedLine7.csv", 4500)
]

speedFolder = "EffMap/"
speedInterpFuncs = {}
mass_flow_min = 2  #kg/s
mass_flow_max = 6 #kg/s

for filename, rpm in speedFiles:
    path = speedFolder + filename
    temp = pd.read_csv(path)
    temp.columns = ['mass_flow', 'pressure_ratio']
    fInterp = interp1d(temp['mass_flow'], temp['pressure_ratio'],kind='linear', fill_value="extrapolate")
    speedInterpFuncs[rpm] = fInterp


def compModel(mass_flow, PRQ):
    rpms = []
    prList = []
    for rpm in sorted(speedInterpFuncs.keys()):
        prVal = speedInterpFuncs[rpm](mass_flow)
        rpms.append(rpm)
        prList.append(prVal)
    rpms = np.array(rpms)
    prList = np.array(prList)
   
    if PRQ < prList.min() or PRQ > prList.max():
        if PRQ < prList.min():
            iLow, iHigh = 0, 1
        else:
            iLow, iHigh = -2, -1
    else:
        iLow = None
        for i in range(len(prList)-1):
            if prList[i] <= PRQ <= prList[i+1]:
                iLow = i
                break
        if iLow is None:
            for i in range(len(prList)-1):
                if prList[i] >= PRQ >= prList[i+1]:
                    iLow = i
                    break
        if iLow is None:
            raise ValueError("Cannot bracket the target pressure ratio with the given data.")
        iHigh = iLow + 1

    rpmLow = rpms[iLow]
    rpmHigh = rpms[iHigh]
    prLow = prList[iLow]
    prHigh = prList[iHigh]
    
    if prHigh == prLow:
        rpm_required = rpmLow
    else:
        rpm_required = rpmLow + (PRQ - prLow) * (rpmHigh - rpmLow) / (prHigh - prLow)
    

    query_point = np.array([mass_flow, PRQ])
    predicted_eff = griddata(effPoints, effValues, query_point, method='linear')
    if predicted_eff is None or np.isnan(predicted_eff):
        predicted_eff = griddata(effPoints, effValues, query_point, method='nearest')
    predicted_eff = float(predicted_eff.item())
    
    if predicted_eff < 0.5:
        predicted_eff = 0.5  
        
    return rpm_required, predicted_eff
