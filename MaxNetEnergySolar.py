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
     c = data[1:100, 2]; k = data[:, 3];
     totalArea = data[:, 4]; suitableArea = data[:, 5];
     # In MWh
     embodiedE = data[:, 6]; operationE = data[:, 7]; availFactor = data[:, 8];
     n = len(c)
     
     print "Size problem ", n
     output = open('test2', 'w')
     res = maximizeTotalNetEnergySolar(c, k, suitableArea, embodiedE, operationE, availFactor)
     for i in range(0, n):
         output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t" + str(res[i]) + "\t"  + 
                      str(netEnergySolar(i,res,c,k,suitableArea,embodiedE,operationE,availFactor)>0)
                      + "\t" + str(installedCapacityDensity(11,res[i])) +"\n")
     output.close()

def maximizeTotalNetEnergySolar(c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5):
    size = len(c)
    def total_net_energy(x):
        res = 0
        for i in range(0, size):
            netE = netEnergyWind(i, x, c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp)
            if netE > 0:
                res = res - netE
        return res
    
    x0 = np.full((size, 1), 10)
    bounds = [(1, 20)] * size
    res = minimize(total_net_energy, x0, bounds=bounds) #, options={'disp': True})
    return res.x

def netEnergySolar(i, n, c, k, suitableArea, embodiedEnergy1MW, operationE, availFactor, cp=0.5): 
    return (suitableArea[i] * installedCapacityDensity(11, n[i]) * (energyPerYear1MW(c[i], k[i], 11, n[i], availFactor[i]) * 25 * (1 - operationE[i]) - embodiedEnergy1MW[i]))  

if __name__ == "__main__":
    sys.exit(main())
