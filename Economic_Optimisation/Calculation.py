import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot

# Generic function for the inequality constraint x[1] + x[2] <= 1
def non_complementary_constraint(i, j):
    def g(x):
        return np.array([-x[i] - x[j] + 1])
    return g
def binary_bounds(n):
    bnds = []
    for i in range(0, n):
        bnds.append((0.0, 1.0))
    return bnds

def loadData(opti_inputs):
     data = genfromtxt(opti_inputs, delimiter='\t', dtype=None)
     lats = data[:, 0]; lon = data[:, 1]
     # Potential suitable area for each tecnology
     area = data[:, 2:5]
     c = data[:, 5]; k = data[:, 6]; ghi = data[:, 7]; dni = data[:, 8]
     embodiedE1y_wind = data[:, 9]
     operationE_wind = data[:, 10]
     avail_wind = data[:, 11]
     return (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)

def loadDataSimpleModel(opti_inputs):
     data = genfromtxt(opti_inputs, delimiter='\t', dtype=None)
     lats = data[:, 0]; lon = data[:, 1]
     # Potential suitable area for each tecnology
     area = data[:, 2:5]
     eff = data[:, 5:8]
     ressources = data[:, 8:11]
     installed_capaciy_density = data[:, 11:14]
     embodiedE1y = data[:, 14:17]
     operationE = data[:, 17:20]
     return (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE)
        
   
# Fixed inputs !
efficiency_PV = 0.18622918619208484
installed_capacity_density_PV = 240
operationE_PV = 0.0097
operationE_CSP = 0.07300000000000001
embodiedE1y_PV = 180.00355555555558
embodiedE1y_CSP_fixed = 231.78966666666662
embodiedE1y_CSP_area = 109.73390740740737
defaultCSP_area = 13.533834586466165

# Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
# x is the vector corresponding to the % of each cells covered by a given technology
def netEnergySimpleModel(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y): 
    if len(area) == 1:
        one = 1
    else:
        one = np.ones(len(x))
    production = x * area * eff * ressources * (365 * 24) * (one - operationE)
    ee = x * area * installed_capaciy_density * embodiedE1y 
    return sum(production - ee)

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
    return x * area * (efficiencyCSP(dni, sm) * dni * (365 * 24) * (1 - operationE_CSP)) - eeCSP
   
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

# Results grid with 3 decision variable per cell (x_wind, x_pv, x_csp)
def writeResultsGrid(output_file, n, lats, lon, res):
    output = open(output_file, 'w')
    for i in range(0, n):
       output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
       resIndex = i * 3; x = np.zeros(3);
       for j in range(0, 3):
           x[j] = res[1][resIndex + j]
       # netE = netEnergy(x, area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
       # output.write(str(netE / 1000) + "\t")
       output.write(str(x[0]) + "\t" + str(x[1]) + "\t" + str(x[2]))
       output.write("\n")
    output.close()
    return

def writeResultsCell(output, lat, lon, res):
    output.write(str(lat) + "\t" + str(lon) + "\t")
    output.write(str(res[0] / 1000))
    for i in range(0, len(res[1])):
        output.write("\t" + str(res[1][i]))
    output.write("\n")
    return