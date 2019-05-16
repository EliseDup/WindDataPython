import sys
import numpy as np
from numpy import genfromtxt, ndarray, array
import math
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
import Calculation
import time
import pulp

def main():
    inputs = 'inputs/simple_sf'
    inputs_total = 'inputs/simple_total'
   
    results_maximiseNetEnergyCell(inputs, 'outputs/test_simple_pulp', True, 10, True)
    results_maximiseNetEnergyCell(inputs, 'outputs/test_simple_scipy', True, 10, False)
    
    #results_maximiseNetEnergyGrid(inputs, 'outputs/test_simple_grid', False, 100)
    
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
        
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, size, pulp):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE) = loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
    print "# Cells:", n
    output = open(output_file, 'w')
    for i in range(9620, n):
        if(i%10000==0):
            print "Progress ", round(float(i)/float(n)*100), "%"
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
        if sum(area[i]) > 1E-3:
           if pulp:
               res = maximiseNetEnergyCell_Pulp(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
           else:
               res = maximiseNetEnergyCell(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
           output.write(str(res[0]/1000) + "\t")
           output.write(str(res[1][0]) + "\t" + str(res[1][1]) + "\t" + str(res[1][2]))
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

# Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
# x is the vector corresponding to the % of each cells covered by a given technology
def netEnergy(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y): 
    if len(area)==1:
        one = 1
    else:
        one = np.ones(len(x))
    production = x * area * eff * ressources * (365 * 24) * (one - operationE)
    ee = x * area * installed_capaciy_density * embodiedE1y 
    return sum(production - ee)
 # Only take the indexes where there is a potential (area > 0) for the optimization problem
def getIndexes(area):
    indexes = []
    for i in range(0,len(area)):
        if area[i] > 0:
            indexes.append(i) 
    return np.array(indexes)  
# Maximise the net energy produced on one cell: only 3 variables
def maximiseNetEnergyCell(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    def obj(x): 
        return -netEnergy(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y)
    # linear_constraint = LinearConstraint([[0, 1, 1]], [0], [1])
    cons = ({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(1, 2)})
    res = minimize(obj, x0=np.ones(3), bounds=Calculation.binary_bounds(3), constraints=[cons], method='trust-constr')
    return (-obj(res.x), res.x)

def maximiseNetEnergyCell_Pulp(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    my_lp_problem = pulp.LpProblem("NE Maximisation Cell", pulp.LpMaximize)
    indexes = getIndexes(area)
    x = []
    for i in indexes:
        x.append(pulp.LpVariable('x'+str(i), lowBound=0, upBound=1, cat='Continuous'))
    # Objective function
    my_lp_problem += netEnergy(x, area[indexes], eff[indexes], ressources[indexes], installed_capaciy_density[indexes], operationE[indexes], embodiedE1y[indexes]), "Z"
    # Constraints
    if (indexes == 1).sum() + (indexes == 2).sum() == 2:
        if (indexes == 0).sum()==1:
            my_lp_problem += x[1] + x[2] <= 1
        else:
            my_lp_problem += x[0] + x[1] <= 1
    
    my_lp_problem.solve()
    
    x_res = np.zeros(3); j = 0;
    for i in indexes:
        x_res[i] = (my_lp_problem.variables()[j].varValue); j = j + 1;
    #for variable in my_lp_problem.variables():
    #    x_res.append(variable.varValue)
    
    return (pulp.value(my_lp_problem.objective), x_res)

# Maximise the net energy produced for the whole grid: 3 variables / grid cell !
def maximiseNetEnergyGrid(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    def obj(x): 
        return -netEnergy(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y)
    
    cons = []
    for i in range(0, (n / 3)):
        cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(i * 3 + 1, i * 3 + 2)})
    res = minimize(obj, x0=np.ones(n), bounds=Calculation.binary_bounds(n), constraints=cons, method='trust-constr')
    return  (-obj(res.x), res)

if __name__ == "__main__":
    sys.exit(main())
