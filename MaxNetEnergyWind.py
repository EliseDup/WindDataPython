import sys
import numpy as np
from numpy import genfromtxt
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc

def main():
     data = genfromtxt('../WindPotentialScala/opti_inputs_wind', delimiter='\t', dtype=None)
     
     lats = data[:, 0]; lon = data[:, 1]
     c = data[:, 2]; k = data[:, 3];
     totalArea = data[:, 4]; suitableArea = data[:, 5];
     # In MWh !!
     embodiedE = data[:, 6]; 
     # Fraction of the energy outputs
     operationE = data[:, 7]; 
     availFactor = data[:, 8];
     n = len(c)
     
     print "Size problem ", n
     output = open('wind', 'w')
     #res = maximizeTotalNetEnergyWind(c, k, suitableArea, embodiedE, operationE, availFactor)
     for i in range(0, n):
         res = maximizeNetEnergyWind(c[i], k[i], suitableArea[i], embodiedE[i], operationE[i], availFactor[i])
         output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t" +
                      str(res[0]) + "\t"+str(res[1]) + "\t"+str(res[2]) + "\t" +
                      str(productionDensity(c[i], k[i], res[1],res[0], availFactor[i])) +"\n") 
     
     output.close()

# Maximise net energy in one cell
def maximizeNetEnergyWind(c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5):
    
    def net_energy(x):
        netE = netEnergyWind(x[0], x[1], c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp)
        if netE > 0: return -netE
        else: return 0
    
    res = minimize(net_energy, x0=(10,11), bounds=[(1,20),(10.0, 16.0)]) #, options={'disp': True})
    
    if -net_energy(res.x) > 0: 
        n = res.x[0]; vr = res.x[1]
    else: n = 0; vr = 0;
    return (n, vr, -net_energy(res.x))

def netEnergyWind(n, vr, c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5): 
    outputs = (suitableArea * installedCapacityDensity(vr, n)) * (energyPerYear1MW(c, k, vr, n, availFactor) * 25 * (1 - operationE))
    inputs = (suitableArea * installedCapacityDensity(vr, n)) *embodiedEnergy1MW          
    return  (outputs-inputs)/1E6

# c = installed capacity density; rp = rated power; d = rotor diameter
# Lambda = Area_rotor / Area_turbine = PI / 4n^2 = PI / (4 *  Rated Power / (CD * D2)) = PI * CD * D^2 / Rated Power 
# n = turbine spacing, depends on CD and Rated Power = Math.sqrt( Rated Power / ( CD * D^2 ) )      
def n(c, rp, d):
    return math.sqrt(rp / (c * d * d))
def lam(c, rp, d): 
    return math.pi / (4 * n(c, rp, d) * n(c, rp, d))
def arrayEfficiency(c, rp, d): 
     # Gustavson extrapolation for arrays of 50x50 wind turbines 
    a = 0.9838; b = 42.5681;
    return a * math.exp(-b * lam(c, rp, d))

# def arrayEfficiencyNInf(n): 
     # Gustavson extrapolation for arrays of 50x50 wind turbines 
#    a = 0.9619; b = 88.9204;
#    return a * math.exp(-b * math.pi / (4*n*n))

def arrayEfficiencyN(n):
    a = 0.9838; b = 42.5681;
    return a * math.exp(-b * math.pi / (4 * n * n))

def capacityFactor(c, k, vr):
    vf = 25.0; vc = 3.0;
    return -math.exp(-math.pow(vf / c, k)) + (3 * math.pow(c, 3) * math.gamma(3.0 / k) / (k * (math.pow(vr, 3) - math.pow(vc, 3)))) * (gammainc(3.0 / k, math.pow(vr / c, k)) - gammainc(3.0 / k, math.pow(vc / c, k)))

def installedCapacityDensity(vr, n, airDensity=1.225, cp=0.5):
    if n==0: return 0
    else: return (0.5 * cp * airDensity * math.pi / 4.0 * math.pow(vr, 3)) / math.pow(n, 2)

def eroi(c, k, vr, n, embodiedEnergy1MW, operationE, area, availFactor, airDensity):
    mw = installedCapacityDensity(vr, n, airDensity) * area
    out = mw * energyPerYear1MW(c, k, vr, n, availFactor) * 25
    return out / (mw * embodiedEnergy1MW + out * operationE)

def energyPerYear1MW(c, k, vr, n, availFactor):
    return capacityFactor(c, k, vr) * arrayEfficiencyN(n) * availFactor * 365 * 24

def productionDensity(c, k, vr, n, availFactor, airDensity=1.225):
    if(n==0): return 0.0
    else: return installedCapacityDensity(vr, n, airDensity) * capacityFactor(c, k, vr) * arrayEfficiencyN(n) * availFactor

if __name__ == "__main__":
    sys.exit(main())
