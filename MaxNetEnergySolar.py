import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot

def main():
    maximiseNetEnergyGrid()
    print "Hello"
    data = genfromtxt('../WindPotentialScala/opti_inputs_solar', delimiter='\t', dtype=None)[0,:]
    embodiedEPV = data[7]; operationEPV = data[8];embodiedECSP_fixed = data[9]; embodiedECSP_variable = data[10]; default_area_CSP = data[11]; operationECSP = data[12];
    irr = range(1,400)
    n = len(irr)
    netEPV = np.zeros(n)
    netECSP = np.zeros(n)
    optiSM = np.zeros(n)
    
    for i in range(0,n):
        netEPV[i] = maximizeNetEnergySolarPV(irr[i], 1, 1, embodiedEPV, operationEPV)[2]
        resCSP = maximizeNetEnergySolarCSP(irr[i], 1, 1, embodiedECSP_fixed, embodiedECSP_variable, default_area_CSP, operationECSP)
        netECSP[i] = resCSP[3] #*1E6/(365*24*25)
        if resCSP[2]: optiSM[i] = resCSP[1]
        
    plt.subplot(3, 1, 1)
    plot(irr,netEPV,irr,netECSP)
    plt.xlabel("Irradiance [W/m2]"); plt.ylabel("Net Power Produced [MW/km2]")
     
    plt.subplot(3, 1, 2)
    plot(irr, optiSM)
    plt.xlabel("Irradiance [W/m2]"); plt.ylabel("Optimal SM [MW/km2]")
     
    # Find the DNI needs for CSP to overcome PV given a fixed DNI
    min_DNI = np.zeros(n)
    for i in range(0,n):
        if netEPV[i]>0 :
            res = find_minimum_DNI(netEPV[i],embodiedECSP_fixed, embodiedECSP_variable, default_area_CSP, operationECSP)
            min_DNI[i] = res[0] 
        else: min_DNI[i] = 0
        
    plt.subplot(3, 1, 3)
    plot(irr, min_DNI)
    plt.xlabel("GHI [W/m2]");plt.ylabel("Minimum DNI [W/m2]")    
    
    plt.show()

def find_minimum_DNI(netEPV, embodiedECSP_fixed, embodiedECSP_variable, default_area_CSP, operationECSP):
    def f(x):
        return x[0]
    # Constraint : netEPV - netCSP > 0
    def constraint(x):
        return netEnergySolarCSP(1, x[1], x[0], 1, 1, embodiedECSP_fixed, embodiedECSP_variable, default_area_CSP, operationECSP) - netEPV
    cons = {'type':'ineq', 'fun': constraint}
    # x = [DNI, SM]
    res = minimize(f, x0=(100, 1.0), bounds=[(0.0, None), (0.5, 4.0)], constraints=cons)
    
    return (res.x[0], res.x[1], netEnergySolarCSP(1, res.x[1], res.x[0], 1, 1, embodiedECSP_fixed, embodiedECSP_variable, default_area_CSP, operationECSP))
    
def maximiseNetEnergyGrid():
     data = genfromtxt('../WindPotentialScala/opti_inputs_solar', delimiter='\t', dtype=None)
     
     lats = data[:, 0]; lon = data[:, 1]
     ghi = data[:, 2]; dni = data[:, 3];
     # Area in km2, suitability factor for PV and CSP
     area = data[:, 4]; sfPV = data[:, 5]; sfCSP = data[:, 6];
     # In MWh
     embodiedEPV = data[:, 7]; 
     # In % of outputs
     operationEPV = data[:, 8];
     # In MWh
     embodiedECSP_fixed = data[:, 9]; embodiedECSP_variable = data[:, 10]; default_area_CSP = data[:, 11]; 
     # In % of outputs
     operationECSP = data[:, 12];
     n = len(ghi)

     print "Size problem ", n
     
     output = open('solar', 'w')
     for i in range(0, n):
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
        if sfPV[i] > 0.0001:
            res = maximiseNetEnergySolar(ghi[i], dni[i], sfPV[i], sfPV[i], area[i], embodiedEPV[i], operationEPV[i], embodiedECSP_fixed[i], embodiedECSP_variable[i], default_area_CSP[i], operationECSP[i], 5.0, 5.0)
            csp_dni = find_minimum_DNI(netEnergySolarPV(1, ghi[i], 1, 1, embodiedEPV[i], operationEPV[i]), embodiedECSP_fixed[i], embodiedECSP_variable[i], default_area_CSP[i], operationECSP[i])[0]
            if(dni[i] >= csp_dni): has_csp = 1.0
            else: has_csp = 0.0
            output.write(str(res[0]+2*res[1]) +"\t" +str(res[0]) + "\t" + str(res[1]) + "\t" + str(res[2]) + "\t" + str(res[3]/(area[i]*sfPV[i])) +"\t" +str(has_csp))
           
        else: 
            output.write("0.0" + "\t" + "0.0 " + "\t" + "0.0 " + "\t" + "0.0"+ "\t" + "0.0" +"\t" + "0.0")
         
        output.write("\n")
     output.close()

# x is the fraction of the suitable area covered with solar panels
# Results in TWh !
def netEnergySolarPV(x, ghi, sf, area, embodiedEnergy1MW, operationE, gcr=5.0): 
    # In MWh / year
    outputs = x * sf * area / gcr * efficiencyPV() * ghi * 365 * 24 * (1 - operationE) * 25
    # Rated power in MW (240 MWp/km2)
    ratedPower = 240 * x * sf * area / gcr
    inputs = embodiedEnergy1MW * ratedPower 
    return (outputs - inputs) / 1E6

def netEnergySolarCSP(x, sm, dni, sf, area, embodiedEnergy1MW, embodiedEnergyArea, defaultArea, operationE, gcr=5.0):
    aperture_area = x * sf * area / gcr
    rated_power = aperture_area * 950 * 0.22 / sm  # In MW !
    outputs = aperture_area * dni * 365 * 24 * efficiencyCSP(dni, sm) * (1 - operationE) * 25
    inputs = embodiedEnergy1MW * rated_power + embodiedEnergyArea * aperture_area / defaultArea
    return (outputs - inputs) / 1E6 

def maximiseNetEnergySolar(ghi, dni, sfPV, sfCSP, area, embodiedEPV, operationEPV, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP, gcrPV, gcrCSP):
    # If net energy is < 0 it needs to choose to put x=0 !
    def net_energy(x):
        netEPV = netEnergySolarPV(x[0], ghi, sfPV, area, embodiedEPV, operationEPV, gcrPV)
        netECSP = netEnergySolarCSP(x[1], x[2], dni, sfCSP, area, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP, gcrCSP)
        return -(netEPV + netECSP)
        
    # Constraint : x[0] + x[1] <= 1
    def constraint(x):
        return -x[0] - x[1] + 1.0
    cons = {'type':'ineq', 'fun': constraint}
    
    # x[0] is the proportion of the suitable area covered with PV, x[1] the proportion covered with CSP, x[2] the solar multiple for the CSP plants
    res = minimize(net_energy, x0=(0.5, 0.5, 1.0), bounds=[(0.0, 1.0), (0.0, 1.0), (0.5, 4.0)], constraints=cons)
    # res = minimize(net_energy, x0=(0.5,0.5), bounds=[(0.0, 1), (0.0, 1)], constraints=cons)
    netEPV = netEnergySolarPV(res.x[0], ghi, sfPV, area, embodiedEPV, operationEPV)
    netECSP = netEnergySolarCSP(res.x[1], res.x[2], dni, sfCSP, area, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP)
    if netEPV > 0.0001: factorPV = res.x[0]
    else: factorPV = 0.0
    if netECSP > 0.0001: 
        factorCSP = res.x[1]
        n = res.x[2]
    else: 
        factorCSP = 0.0
        n = 0.0
    
    return (factorPV, factorCSP, n, -net_energy(res.x))

def maximizeTotalNetEnergySolarPV(ghi, suitableArea, embodiedEnergy1MW, operationE):
    size = len(ghi)
    def total_net_energy(x):
        res = 0
        for i in range(0, size):
            netE = netEnergySolarPV(x[i], ghi[i], suitableArea[i], embodiedEnergy1MW[i], operationE[i])
            if netE > 0:
                res = res - netE
        return res
    
    x0 = np.full((size, 1), 0)
    bounds = [(0, 1)] * size
    res = minimize(total_net_energy, x0, bounds=bounds)  # , options={'disp': True})
    return res.x

def maximizeNetEnergySolarPV(ghi, sf, area, embodiedEnergy1MW, operationE):
    def net_energy(x):
        return -netEnergySolarPV(x[0], ghi, sf, area, embodiedEnergy1MW, operationE)
        
    res = minimize(net_energy, 0, bounds=[(0, 1)])  # , options={'disp': True})
    return (res.x[0], -net_energy(res.x) > 0, -net_energy(res.x))

# Choose the optimal solar mutiple, and the fraction of the suitable area that should be occupied
def maximizeNetEnergySolarCSP(dni, sf, area, embodiedEnergy1MW_fixed, embodiedEnergyArea, defaultArea, operationE):
    def net_energy(x):
        if dni == 0: return 0
        return -netEnergySolarCSP(x[0], x[1], dni, sf, area, embodiedEnergy1MW_fixed, embodiedEnergyArea, defaultArea, operationE)
    
    res = minimize(net_energy, x0=(1, 2.7), bounds=[(0, 1), (0.5, 4.0)])  # , options={'disp': True})
    return (res.x[0], res.x[1], -net_energy(res.x) > 0, -net_energy(res.x))

def efficiencyPV(designEfficiency=0.24, performanceRatio=0.81, degradationRate=0.36 / 100, lifeTime=25):
    if (degradationRate == 0): return designEfficiency * performanceRatio
    else: return designEfficiency * performanceRatio * ((1.0 - math.pow(1.0 - degradationRate, lifeTime)) / degradationRate) / lifeTime
  
# DNI in W/m2
def efficiencyCSP(dni, sm):
    def a(sm):
        return -1.578 * sm + 11.17
    def b(sm):
        return 10.65 * sm - 66.33
    if dni == 0: return 0
    else: return (a(sm) * math.log(dni * 8.76) + b(sm)) / 100.0

if __name__ == "__main__":
    sys.exit(main())
