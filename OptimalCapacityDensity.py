import sys
import numpy as np
from numpy import genfromtxt
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc

def main():
     # Load capacity factors and area of each cell
     data = genfromtxt('../WindPotential/world_c_k', delimiter='\t', dtype=None)
     
     lats = data[:, 0]; lon = data[:, 1]
     c = data[:, 2]; k = data[:, 3];
     totalArea = data[:, 4]; suitableArea = data[:, 5];
     # In MWh
     embodiedE = data[:, 6]
     n = len(lats)
     output = open('res_world_net_energy', 'w')
     for i in range(0, n):        
        if(suitableArea[i] > 0):
            output.write(str(lats[i]) + '\t' + str(lon[i]))
            if i % 1000 == 0: print i
            for e in np.linspace(5,12,8):
                res = maximizeNetEnergy(e,c[i], k[i], suitableArea[i], embodiedE[i])
                output.write('\t' + str(res[0]) + '\t'+  str(res[1]) + '\t' + str(res[3]))
            output.write('\n')
     
     output.close()

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

def arrayEfficiencyN(n):
    a = 0.9838; b = 42.5681;
    return a * math.exp(-b * math.pi / (4*n*n))

def capacityFactor(k, c, vr):
    vf = 25.0; vc = 3.0;
    return -math.exp(-math.pow(vf/c,k)) + (3*math.pow(c,3)*math.gamma(3.0/k) / (k*(math.pow(vr,3) - math.pow(vc,3)))) * (gammainc(3.0/k,math.pow(vr/c,k)) - gammainc(3.0/k,math.pow(vc/c,k)))

def maximizeNetEnergy(eroi_min, c, k, area, embodiedEnergy1MW, cp=0.5):
    # W/m^2 == MW/km^2
    def installedCapacityDensity(x):
        return (0.5*cp*1.225*math.pi/4.0*math.pow(x[0],3)) / math.pow(x[1],2)
    # x = (rated wind speed, turbine spacing n)
    def net_energy(x):
        if -eroi(x) >= eroi_min:
            return -(area * installedCapacityDensity(x) * (capacityFactor(k, c, x[0])*arrayEfficiencyN(x[1])*365*24*25 - embodiedEnergy1MW))
        else:
            return 1000
        
    def eroi(x):
        return -(capacityFactor(k, c, x[0])*365*24*25*arrayEfficiencyN(x[1]) / embodiedEnergy1MW)
    
    res = minimize(net_energy, x0=(8,20), bounds=[(8.0,14.0),(1, 50)])  # options={'disp': True})
    
    return (res.x[0], res.x[1], -eroi(res.x), net_energy(res.x)<1000)

# All in km^2 and MW ?
def maximizeOneDensity(eroi_min, cf, area, diameters, ratedPower, embodiedEnergy1MW):

    size = len(cf)
    print "Size = ", size
    
    def fun_EROI_min(x):
        power = 0
        for j in range(0, size):
            if(cf[j] * arrayEfficiency(x, ratedPower[j], diameters[j]) * 20 * 365 * 24 >= eroi_min * embodiedEnergy1MW[j]):
               power += cf[j] * area[j] * x * arrayEfficiency(x, ratedPower[j], diameters[j])
        return -power

    res = minimize(fun_EROI_min, x0=5, bounds=[(0.1, None)])  # , options={'disp': True})
    
    print eroi_min, '\t', 'Optimal CD', '\t', res.x, '\t', 'Power Output', '\t', -fun_EROI_min(res.x) / 1E6, '\t', 'TW'
    
    return res.x[0]

def maximizePowerInCell(eroi_min, cf, area, diameters, ratedPower, embodiedEnergy1MW):
 
    def power_output(x):
        if(cf * arrayEfficiency(x, ratedPower, diameters) * 20 * 365 * 24 >= eroi_min * embodiedEnergy1MW):
            return -cf * area * x * arrayEfficiency(x, ratedPower, diameters)
        return 1000
    
    res = minimize(power_output, x0=1, bounds=[(0.1, None)])  # options={'disp': True})
    
    return (res.x[0], (cf * arrayEfficiency(res.x[0], ratedPower, diameters) * 20 * 365 * 24 >= eroi_min * embodiedEnergy1MW))

def optimizeSpacing(eroi_min, cf, area, embodiedEnergy1MW):
 
    def power_output(x):
        if(cf * arrayEfficiencyN(x) * 25 * 365 * 24  >= eroi_min * embodiedEnergy1MW):
            return -cf * area * 320 / (x*x) * arrayEfficiencyN(x)
        return 1000
    
    res = minimize(power_output, x0=10, bounds=[(0.1, None)])  # options={'disp': True})
    
    return (res.x[0], (cf * arrayEfficiencyN(res.x[0]) * 20 * 365 * 24 >= eroi_min * embodiedEnergy1MW))

def maximizeEROIInCell(cf, area, diameters, ratedPower, embodiedEnergy1MW):
    
    def eroi(x): 
        return -cf * area * x * arrayEfficiency(x, ratedPower, diameters) * 20 * 365 * 24 / (x * area * embodiedEnergy1MW)
     
    res = minimize(eroi, x0=1, bounds=[(0.01, None)])  # options={'disp': True})
    
    return (res.x[0], -eroi(res.x[0]))
    
def maximizeDensityPerCell(eroi_min, start, cf, area, diameters, ratedPower, embodiedEnergy1MW):
    
    size = len(cf)
    print "Size = ", size
    
    def power_output(x):
        p = 0
        for i in range(0, size):
            # In order to reach min eroi we need :
            # CD * Area * CF * efficiency * 20 years >= EROI_min * Energy/MW * MW Installed(=CD*Area) ==> 1MW * CF * efficiency * 20 years >= EROI_min * Energy 1MW 
            if(cf[i] * arrayEfficiency(x[i], ratedPower[i], diameters[i]) * 20 * 365 * 24 >= eroi_min * embodiedEnergy1MW[i]):
                  p += cf[i] * area[i] * x[i] * arrayEfficiency(x[i], ratedPower[i], diameters[i])
        return -p
    
    bnds = list()
    for i in range(0, size): 
        bnds.append((0.1, None))
    res = minimize(power_output, x0=np.ones((size, 1)), bounds=bnds, options={'disp': True})
    
    print 'Max Wind Power', -power_output(res.x) / 1E6, 'TW'
    # Write optimization results (overwrite the file if we ran a new simulation)
    if(start == 0): output = open('output_' + str(int(eroi_min)), 'w')
    else: output = open('output_' + str(int(eroi_min)), 'a')
    for i in range(0, size):
        used = cf[i] * arrayEfficiency(res.x[i], ratedPower[i], diameters[i]) * 20 * 365 * 24 >= eroi_min * embodiedEnergy1MW[i]
        output.write(str(start + i) + "\t" + str(res.x[i]) + "\t" + str(used) + "\n")
    output.close()
   

if __name__ == "__main__":
    sys.exit(main())
