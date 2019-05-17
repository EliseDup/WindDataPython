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

# Pulp is faster than scipy.optimize.minimize !
def main():
    #results_maximiseNetEnergyCell(Calculation.inputs_simple, 'outputs/netE_Simple_Cell', False, 1000, True)
    #results_maximiseNetEnergyCell(inputs, 'outputs/test2', False, 100, False)
    results_maximiseNetEnergyGrid(Calculation.inputs_simple, 'outputs/netE_Simple_Grid', False, 10000, True)
    #results_maximiseNetEnergyGrid(inputs, 'outputs/test4', False, 100, False)
       
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, size, pulp):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModel(opti_inputs)
    if total: 
        n = len(lats) 
    else:
        n = size
    print "# Cells:", n
    
    output = open(output_file, 'w')
    total = 0
    for i in range(0, n):
        if(i % (n/10) == 0):
            print "Progress ", round(float(i)/float(n)*100), "%"
        if pulp:
            res = maximiseNetEnergy_Pulp(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i], keMax[i])
        else:
            res = maximiseNetEnergyCell(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
        total += res[0]
        Calculation.writeResultsCell(output, lats[i], lon[i], res)
    output.close()
    print "Results Grid ", total / 1E6, " TWh "
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"
    
def results_maximiseNetEnergyGrid(opti_inputs, output_file, total, size, pulp):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModel(opti_inputs)
    if total: 
        n = len(lats)
    else: 
        n = size
    print "# Cells:", n
    if pulp:
        res = maximiseNetEnergy_Pulp(area[0:n, :], eff[0:n, :], ressources[0:n, :], installed_capaciy_density[0:n, :], operationE[0:n, :], embodiedE1y[0:n, :], keMax[0:n])
    else:
        res = maximiseNetEnergyGrid(area[0:n, :], eff[0:n, :], ressources[0:n, :], installed_capaciy_density[0:n, :], operationE[0:n, :], embodiedE1y[0:n, :])
    print "Results Grid ", res[0] / 1E6, " TWh "
    Calculation.writeResultsGrid(output_file, n, lats, lon, res)
    print "Optimization for the whole grid ends in ", (time.time() - t0), " seconds"

def maximiseNetEnergy_Pulp(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    
    if(sum(area) < 0.001):
        return(0.0, np.zeros(n))
    else:
        my_lp_problem = pulp.LpProblem("NE Maximisation Cell", pulp.LpMaximize)
        indexes = Calculation.getIndexesSimpleModel(area)
        x = []
        for i in indexes[0]:
            x.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Continuous'))
        # Objective function
        my_lp_problem += Calculation.netEnergySimpleModel(x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]]), "Z"
        # Constraints: non complementary constraints for PV and CSP
        for i in indexes[1]:
            my_lp_problem += x[i[0]] + x[i[1]] <= 1
        # Constraints: max KE generation for wind 
        for i in indexes[2]:
            my_lp_problem += x[i]*area[i]*ressources[i]*eff[i] <= keMax[i]
            
        my_lp_problem.solve()
        
        x_res = np.zeros(n); j = 0;
        for i in indexes[0]:
            x_res[i] = (my_lp_problem.variables()[j].varValue); j = j + 1;
            
        return (pulp.value(my_lp_problem.objective), x_res)

# Maximise the net energy produced on one cell: only 3 variables
def maximiseNetEnergyCell(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    if(sum(area) < 0.001):
        return (0, np.zeros(3))
    else:
        # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
        indexes = Calculation.getIndexesSimpleModel(area)
        n = len(indexes[0])
        def obj(x):
            return -Calculation.netEnergySimpleModel(x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]])
        cons = []
        for i in indexes[1]:
            cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(i[0], i[1])})
        res = minimize(obj, x0=np.ones(n), bounds=Calculation.binary_bounds(n), constraints=cons, method='trust-constr')
        x_res = np.zeros(3); j = 0;
        for i in indexes[0]:
            x_res[i] = res.x[j]; j = j + 1;
        
        return (-obj(res.x), x_res)

# Maximise the net energy produced for the whole grid: 3 variables / grid cell !
def maximiseNetEnergyGrid(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    def obj(x): 
        return -Calculation.netEnergySimpleModel(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y)
    cons = []
    for i in range(0, (n / 3)):
        cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(i * 3 + 1, i * 3 + 2)})
    res = minimize(obj, x0=np.ones(n), bounds=Calculation.binary_bounds(n), constraints=cons, method='trust-constr')
    return  (-obj(res.x), res.x)

if __name__ == "__main__":
    sys.exit(main())
