import sys
import numpy as np
from numpy import genfromtxt
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc

def main():
     data = genfromtxt('../WindPotentialScala/opti_inputs_solar', delimiter='\t', dtype=None)
    
     lats = data[:, 0]; lon = data[:, 1]
     ghi = data[:, 2]; dni = data[:, 3];
     totalArea = data[:, 4]; suitableAreaPV = data[:, 5]; suitableAreaCSP = data[:, 6];
     # In MWh
     embodiedEPV = data[:, 7]; 
     # In % of outputs
     operationEPV = data[:, 8];
     # In MWh
     embodiedECSP_fixed = data[:, 9]; embodiedECSP_variable = data[:, 10]; defautl_area_CSP = data[:, 11]; 
     # In % of outputs
     operationECSP = data[:, 12];
     i = 200
     # print ghi[i], dni[i], suitableAreaPV[i]
     res = maximiseNetEnergySolar(i, ghi, dni, suitableAreaPV, suitableAreaPV, embodiedEPV, operationEPV, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP)
     print maximizeNetEnergySolarCSP(i, dni, suitableAreaPV, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP)
     print maximizeNetEnergySolarPV(i, ghi, suitableAreaPV, embodiedEPV, operationEPV)
     print res
     print "---"
         
     n = len(ghi)
     print "Size problem ", n
     output = open('solar', 'w')
     for i in range(0, n):
         res = maximiseNetEnergySolar(i, ghi, dni, suitableAreaPV, suitableAreaPV, embodiedEPV, operationEPV, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP)
         output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t" + str(res[0]) + "\t" + str(res[1]) + "\t" + str(res[2]) + "\t" + str(res[3]) + "\t" + str(res[4]) + 
                      "\n")
     output.close()

# x is the fraction of the suitable area covered with solar panels
# Results in TWh !
def netEnergySolarPV(i, x, ghi, suitableArea, embodiedEnergy1MW, operationE): 
    # In MWh / year
    outputs = x * suitableArea[i] / 5.0 * efficiencyPV() * ghi[i] * 365 * 24 * (1 - operationE[i]) * 25
    # Rated power in MW (240 MWp/km2)
    ratedPower = 240 * x * suitableArea[i] / 5.0
    inputs = embodiedEnergy1MW[i] * ratedPower 
    # print "PV", x, outputs, inputs
    return (outputs - inputs)/1E6

def netEnergySolarCSP(i, x, sm, dni, suitableArea, embodiedEnergy1MW, embodiedEnergyArea, defaultArea, operationE):
    aperture_area = x * suitableArea[i] / 7.5
    rated_power = aperture_area * 950 * 0.22 / sm  # In MW !
    outputs = aperture_area * dni[i] * 365 * 24 * efficiencyCSP(dni[i], sm) * (1 - operationE[i]) * 25
    inputs = embodiedEnergy1MW[i] * rated_power + embodiedEnergyArea[i] * aperture_area / defaultArea[i]
    # print "CSP", x, sm, outputs, inputs
    return (outputs - inputs)/1E6

def maximiseNetEnergySolar(i, ghi, dni, suitableAreaPV, suitableAreaCSP, embodiedEPV, operationEPV, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP):
    # If net energy is < 0 it needs to choose to put x=0 !
    def net_energy(x):
        netEPV = netEnergySolarPV(i, x[0], ghi, suitableAreaPV, embodiedEPV, operationEPV)
        netECSP = netEnergySolarCSP(i, x[1], x[2], dni, suitableAreaCSP, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP)
        return -(netEPV + netECSP)
        
    # Constraint : x[0] + x[1] <= 1
    def constraint(x):
        return -x[0] - x[1] + 1.0
    cons = {'type':'ineq', 'fun': constraint}
   
    # x[0] is the proportion of the suitable area covered with PV, x[1] the proportion covered with CSP, x[2] the solar multiple for the CSP plants
    res = minimize(net_energy, x0=(0.5,0.5,1.0), bounds=[(0.0, None), (0.0, None), (0.5, 10.0)], constraints=cons)
    #res = minimize(net_energy, x0=(0.5,0.5), bounds=[(0.0, 1), (0.0, 1)], constraints=cons)
     
    netEPV = netEnergySolarPV(i, res.x[0], ghi, suitableAreaPV, embodiedEPV, operationEPV)
    netECSP = netEnergySolarCSP(i, res.x[1], 2.7, dni, suitableAreaCSP, embodiedECSP_fixed, embodiedECSP_variable, defautl_area_CSP, operationECSP)
    return (res.x[0], res.x[1], res.x[2], netEPV > 0, netECSP > 0, -net_energy(res.x))

def maximizeTotalNetEnergySolarPV(ghi, suitableArea, embodiedEnergy1MW, operationE):
    size = len(ghi)
    def total_net_energy(x):
        res = 0
        for i in range(0, size):
            netE = netEnergySolarPV(i, x[i], ghi, suitableArea, embodiedEnergy1MW, operationE)
            if netE > 0:
                res = res - netE
        return res
    
    x0 = np.full((size, 1), 0)
    bounds = [(0, 1)] * size
    res = minimize(total_net_energy, x0, bounds=bounds)  # , options={'disp': True})
    return res.x

def maximizeNetEnergySolarPV(i, ghi, suitableArea, embodiedEnergy1MW, operationE):
    def net_energy(x):
        return -netEnergySolarPV(i, x[0], ghi, suitableArea, embodiedEnergy1MW, operationE)
        
    res = minimize(net_energy, 0, bounds=[(0, 1)])  # , options={'disp': True})
    return (res.x[0], -net_energy(res.x)>0, -net_energy(res.x))

# Choose the optimal solar mutiple, and the fraction of the suitable area that should be occupied
def maximizeNetEnergySolarCSP(i, dni, suitableArea, embodiedEnergy1MW_fixed, embodiedEnergyArea, defaultArea, operationE):
    def net_energy(x):
        if dni[i] == 0: return 0
        return -netEnergySolarCSP(i, x[0], x[1], dni, suitableArea, embodiedEnergy1MW_fixed, embodiedEnergyArea, defaultArea, operationE)
    
    res = minimize(net_energy, x0=(1, 2.7), bounds=[(0, 1), (0.5, 10)])  # , options={'disp': True})
    return (res.x[0], res.x[1], -net_energy(res.x)>0, -net_energy(res.x))

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
