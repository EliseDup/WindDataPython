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
     c = data[1:100, 2]; k = data[:, 3];
     totalArea = data[:, 4]; suitableArea = data[:, 5];
     # In MWh
     embodiedE = data[:, 6]; operationE = data[:, 7]; availFactor = data[:, 8];
     n = len(c)
     
     print "Size problem ", n
     output = open('test2', 'w')
     res = maximizeTotalNetEnergyWind(c, k, suitableArea, embodiedE, operationE, availFactor)
     for i in range(0, n):
         output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t" + str(res[i]) + "\t"  + 
                      str(maximizeNetEnergyWind(i, c, k, suitableArea, embodiedE, operationE, availFactor)) + "\t"+
                      str(netEnergyWind(i,res[i],c,k,suitableArea,embodiedE,operationE,availFactor)>0) + "\t" + 
                      str(installedCapacityDensity(11,res[i])) +"\n") 
     
     output.close()

# Maximize the net energy of N cells
def maximizeTotalNetEnergyWind(c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5):
    size = len(c)
    def total_net_energy(x):
        res = 0
        for i in range(0, size):
            netE = netEnergyWind(i, x[i], c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp)
            if netE > 0:
                res = res - netE
        return res
    
    x0 = np.full((size, 1), 10)
    bounds = [(1, 20)] * size
    res = minimize(total_net_energy, x0, bounds=bounds) #, options={'disp': True})
    return res.x

# Maximise net energy in one cell
def maximizeNetEnergyWind(i, c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5):
    def net_energy(x):
        netE = netEnergyWind(i, x, c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp)
        if netE > 0: return -netE
        else: return 0
    res = minimize(net_energy, 10, bounds=[(1,20)]) #, options={'disp': True})
    return res.x[0]

def netEnergyWind(i, n, c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5): 
    return (suitableArea[i] * installedCapacityDensity(11, n) * (energyPerYear1MW(c[i], k[i], 11, n, availFactor[i]) * 25 * (1 - operationE[i]) - embodiedEnergy1MW[i]))  

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
    return (0.5 * cp * airDensity * math.pi / 4.0 * math.pow(vr, 3)) / math.pow(n, 2)

def eroi(c, k, vr, n, embodiedEnergy1MW, operationE, area, availFactor, airDensity):
    mw = installedCapacityDensity(vr, n, airDensity) * area
    out = mw * energyPerYear1MW(c, k, vr, n, availFactor) * 25
    return out / (mw * embodiedEnergy1MW + out * operationE)

def energyPerYear1MW(c, k, vr, n, availFactor):
    return capacityFactor(c, k, vr) * arrayEfficiencyN(n) * availFactor * 365 * 24

def productionDensity(c, k, vr, n, availFactor, airDensity):
    return installedCapacityDensity(vr, n, airDensity) * capacityFactor(c, k, vr) * arrayEfficiencyN(n) * availFactor

if __name__ == "__main__":
    sys.exit(main())
