import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot
# from scipy.optimize import LinearConstraint
import time

# Fixed inputs !
efficiency_PV = 0.18622918619208484
installed_capacity_density_PV = 240
operationE_PV = 0.0097
operationE_CSP = 0.07300000000000001
embodiedE1y_PV = 180.00355555555558
embodiedE1y_CSP_fixed = 231.78966666666662
embodiedE1y_CSP_area = 109.73390740740737
defaultCSP_area = 13.533834586466165

def main():
    inputs = '../../WindPotentialScala/params'
    vr = genfromtxt('../vr', delimiter='\t', dtype=None)
    output = open('cf', 'w')
    n = len(vr)
    (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)=loadData(inputs)
    for i in range (0,n):
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t" + str(vr[i]) + "\t" + str(capacityFactor(c[i], k[i], vr[i])) + "\n" )
    output.close()
    #inputs_total = '../WindPotentialScala/simple_total'
    # results_maximiseNetEnergyCell(inputs_total, 'results_cell_total', True, 500)
    #results_maximiseNetEnergyCell(inputs, 'results_cell_params_all_end', True, 500)
    
def loadData(opti_inputs):
     data = genfromtxt(opti_inputs, delimiter='\t', dtype=None)
     lats = data[:, 0]; lon = data[:, 1]
     # Potential suitable area for each tecnology
     area = data[:, 2:5]
     c = data[:, 5]; k = data[:, 6]; ghi = data[:, 7]; dni = data[:, 8]
     embodiedE1y_wind = data[:, 9]
     operationE_wind = data[:,10]
     avail_wind = data[:,11]
     return (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)
        
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind) = loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
    print "# Cells:", n
    output = open(output_file, 'w')
    for i in range(24100, n):
        if(i % 1000 == 0):
            print "Progress ", round(float(i) / float(n) * 100), "%"
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
        if sum(area[i]) > 0:
           res = maximiseNetEnergyCell(area[i], c[i], k[i], ghi[i], dni[i], embodiedE1y_wind[i],operationE_wind[i], avail_wind[i])
           output.write(str(res[0] / 1000) + "\t")
           output.write(str(res[1].x[0]) + "\t" + str(res[1].x[1]) + "\t" + str(res[1].x[2]) + "\t" + str(res[1].x[3])+ "\t" + str(res[1].x[4])+ "\t" + str(res[1].x[5]))
           # Write Cf as it is difficult to calculate afterwards
           output.write("\t" + str(capacityFactor(c[i], k[i], res[1].x[3])))
        else:
           output.write("0.0" + "\t" + "0.0" + "\t" + "0.0" + "\t" + "0.0" + "\t" + "0.0")
        output.write("\n")
    output.close()
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"
 
# Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
# x is the vector corresponding to the % of each cells covered by a given technology
def netEnergy(x, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind):
    return netEnergyWind(x[0], area[0], x[3], x[4], c, k, embodiedE1y_wind, operationE_wind, avail_wind) + netEnergyPV(x[1], area[1], ghi) + netEnergyCSP(x[2], area[2], dni, x[5])

def netEnergyWind(x, area, vr, n, c, k, embodiedE1y_wind, operationE_wind, avail_wind):
    return x * area * installedCapacityDensityWind(vr, n) * ((capacityFactor(c, k, vr) * arrayEfficiency(n) * avail_wind * (365 * 24) * (1 - operationE_wind)) - embodiedE1y_wind)
def netEnergyPV(x, area, ghi):
    return x * area * (efficiency_PV * ghi * (365 * 24) * (1 - operationE_PV) - installed_capacity_density_PV * embodiedE1y_PV)
def netEnergyCSP(x, area, dni, sm):
    eeCSP = embodiedE1y_CSP_fixed * ratedPowerCSP(sm, x * area) + embodiedE1y_CSP_area * x * area / defaultCSP_area
    return x * area * (efficiencyCSP(dni,sm) * dni * (365 * 24) * (1 - operationE_CSP)) - eeCSP
   
# WIND
def capacityFactor(c, k, vr):
    vf = 25.0; vc = 3.0;
    return -math.exp(-math.pow(vf / c, k)) + (3 * math.pow(c, 3) * math.gamma(3.0 / k) / (k * (math.pow(vr, 3) - math.pow(vc, 3)))) * (gammainc(3.0 / k, math.pow(vr / c, k)) - gammainc(3.0 / k, math.pow(vc / c, k)))
def arrayEfficiency(n):
    a = 0.9838; b = 42.5681;
    return a * math.exp(-b * math.pi / (4 * n * n))
# W/m2
def installedCapacityDensityWind(vr, n, airDensity=1.225, cp=0.5):
    if n == 0: return 0
    else: return (0.5 * cp * airDensity * math.pi / 4.0 * math.pow(vr, 3)) / math.pow(n, 2)

# SOLAR
def lifeTimeEfficiency(eff, pr, degradationRate, lifetime):
    if degradationRate == 0: return eff * pr
    else: return  eff * pr * ((1.0 - math.pow(1.0 - degradationRate, lifetime)) / degradationRate) / lifetime
# CSP
def efficiencyCSP(dni, sm):
    def a(sm):
        return -1.578 * sm + 11.17
    def b(sm):
        return 10.65 * sm - 66.33
    if dni == 0: return 0
    else: return lifeTimeEfficiency((a(sm) * math.log(dni * 8.76) + b(sm)) / 100.0, 1.0, 0.2 / 100, 30)
def ratedPowerCSP(sm, area):
    return area * 950 * 0.22 / sm 

# Maximise the net energy produced on one cell: 3 variables x_ij + vr_i + n_i + SM_i
def maximiseNetEnergyCell(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    def obj(x): 
        return -netEnergy(x, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)
    
    cons = ({'type': 'ineq', 'fun' : lambda x:-x[1] - x[2] + 1})
    res = minimize(obj, x0=(1.0, 1.0, 1.0, 11.0, 10.0, 2.7), bounds=[(0, 1.0), (0, 1.0), (0, 1.0), (10.0, 16.0), (1.0, 20.0), (1.0, 4.0)], constraints=[cons], method='trust-constr')
    return (-obj(res.x), res)

if __name__ == "__main__":
    sys.exit(main())
