import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc

import time

def main():
    inputs = 'simple_sf' # '../WindPotentialScala/simple_sf'
    #inputs_total = '../WindPotentialScala/simple_total'
    #results_maximiseNetEnergyCell(inputs, 'test', False, 500)
    results_maximiseNetEnergyGrid(inputs, 'results_grid', False, 1000)
    
def loadData(opti_inputs):
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
        
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE) = loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
    print "# Cells:", n
    output = open(output_file, 'w')
    for i in range(0, n):
        if(i%10000==0):
            print "Progress ", round(float(i)/float(n)*100), "%"
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
        if sum(area[i]) > 0:
           res = maximiseNetEnergyCell(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
           output.write(str(res[0]/1000) + "\t")
           output.write(str(res[1].x[0]) + "\t" + str(res[1].x[1]) + "\t" + str(res[1].x[2]))
        else:
           output.write("0.0" + "\t" +"0.0" + "\t" +"0.0" + "\t" + "0.0")
        output.write("\n")
    output.close()
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"
    
def results_maximiseNetEnergyGrid(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE) = loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
    print "# Cells:", n
    res = maximiseNetEnergyGrid(area[0:n, :], eff[0:n, :], ressources[0:n, :], installed_capaciy_density[0:n, :], operationE[0:n, :], embodiedE1y[0:n, :])
    print "Results Grid ", res[0]/1E6, " TWh "
    output = open(output_file, 'w')
    for i in range(0, n):
       output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
       resIndex = i*3; x = np.zeros(3);
       for j in range(0, 3):
           x[j] = res[1].x[resIndex + j]
       netE = netEnergy(x, area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
       output.write(str(netE/1000) + "\t")
       output.write(str(x[0]) + "\t" + str(x[1]) + "\t" + str(x[2]))
       output.write("\n")
    output.close()
    print "Optimization for the whole grid ends in ", (time.time() - t0), " seconds"

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

# Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
# x is the vector corresponding to the % of each cells covered by a given technology
def netEnergy(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y): 
    production = x * area * eff * ressources * (365 * 24) * (np.ones(len(x)) - operationE)
    ee = x * area * installed_capaciy_density * embodiedE1y 
    return sum(production - ee)
       
# Maximise the net energy produced on one cell: only 3 variables
def maximiseNetEnergyCell(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    def obj(x): 
        return -netEnergy(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y)
    # linear_constraint = LinearConstraint([[0, 1, 1]], [0], [1])
    cons = ({'type': 'ineq', 'fun' : non_complementary_constraint(1, 2)})
    res = minimize(obj, x0=np.ones(3), bounds=binary_bounds(3), constraints=[cons], method='trust-constr')
    return (-obj(res.x), res)

# Maximise the net energy produced for the whole grid: 3 variables / grid cell !
def maximiseNetEnergyGrid(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    def obj(x): 
        return -netEnergy(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y)
    
    cons = []
    for i in range(0, (n / 3)):
        cons.append({'type': 'ineq', 'fun' : non_complementary_constraint(i * 3 + 1, i * 3 + 2)})
    res = minimize(obj, x0=np.ones(n), bounds=binary_bounds(n), constraints=cons, method='trust-constr')
    return  (-obj(res.x), res)

if __name__ == "__main__":
    sys.exit(main())
