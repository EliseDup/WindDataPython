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
    inputs = Calculation.inputs_simple
    
    #results_maximiseNetEnergyGrid(Calculation.inputs_simple_onshore_wind, 'neteSimple_onshorewind', True, 0, 100, True, "onshore_wind", True)
    #results_maximiseNetEnergyGrid(Calculation.inputs_simple_offshore_wind, 'neteSimple_offshorewind', True, 0, 100, True, "offshore_wind", True)
    #results_maximiseNetEnergyGrid(Calculation.inputs_simple_pv, 'neteSimple_pv', True, 0, 100, True, "pv", False)
    #results_maximiseNetEnergyGrid(Calculation.inputs_simple_csp, 'neteSimple_csp', True, 0, 100, True, "csp", False)
    
    results_maximiseNetEnergyGrid(inputs, 'netESimple_grid', False, 0, 10000)
      
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, start, size, singleTechnology=False, nameTechnology="", wind=False):
    t0 = time.time()
    if singleTechnology:
        (lats, lon, total_area, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModelOneTech(opti_inputs)
    else:
        (lats, lon, total_area, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModel(opti_inputs)
 
    if total: 
        n = len(lats); start = 0;
    else:
        n = size
    print "# Cells:", n
    
    output = open(output_file, 'w')
    total = 0
    for i in range(start, n+start):
        if(n > 10 and (i-start) % (n/10) == 0):
            print "Progress ", round(float(i-start)/float(n)*100), "%"
            if singleTechnology:
                res = maximiseNetEnergyTechnology(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i], keMax[i], wind, nameTechnology)
            else:
                res = maximiseNetEnergy(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i], keMax[i])
        total += res[0]
        Calculation.writeResultsCell(output, lats[i], lon[i], res)
    output.close()
    print "Results Grid ", total / 1E6, " TWh "
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"
    
def results_maximiseNetEnergyGrid(opti_inputs, output_file, total, start, size, singleTechnology=False, nameTechnology="", wind=False):
    t0 = time.time()
    if singleTechnology:
        (lats, lon, total_area, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModelOneTech(opti_inputs)
    else:
        (lats, lon, total_area, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModel(opti_inputs)
    if total: 
        n = len(lats); start = 0;
    else: 
        n = size
    print "# Cells:", n
    if singleTechnology:
        res = maximiseNetEnergyTechnology(area[start:n+start], eff[start:n+start], ressources[start:n+start], installed_capaciy_density[start:n+start], operationE[start:n+start], embodiedE1y[start:n+start], keMax[start:n+start], wind,nameTechnology)
    else:
        res = maximiseNetEnergy(area[start:n+start, :], eff[start:n+start, :], ressources[start:n+start, :], installed_capaciy_density[start:n+start, :], operationE[start:n+start, :], embodiedE1y[start:n+start, :], keMax[start:n+start])
    print "Results Grid ", res[0] / 1E6, " TWh "
    #Calculation.writeDetailedResultsGrid(output_file, start, n, lats, lon, res, area, total_area)
    if singleTechnology:
        nX = 1
    else: 
        nX = 3
    Calculation.writeResultsGrid(output_file, start, n, lats, lon, res, nX)
    
    print "Optimization for ",n,"cells in ", (time.time() - t0), " seconds"

def maximiseNetEnergy(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    if(sum(area) == 0):
        return(0.0, np.zeros(n))
    else:
        my_lp_problem = pulp.LpProblem("Net Energy Maximisation Simple", pulp.LpMaximize)
        indexes = Calculation.getIndexesSimpleModel(area)
        nX = len(indexes[0])
        if(nX > 3):
            print nX, " decision variables"
        x = []
        
        for i in indexes[0]:
            x.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Continuous'))
        # Objective function
        my_lp_problem += Calculation.netEnergySimpleModel(x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]]), "Z"
        # Constraints: non complementary constraints for PV and CSP
        for i in indexes[1]:
            my_lp_problem += x[i[0]] + x[i[1]] <= 1
        # Constraints: max KE generation for wind (ke in W/m^2 !!)
        for i in indexes[2]:
            if len(area) == 3:
                ke = keMax
            else:
                ke = keMax[i[1]/3]
            my_lp_problem += x[i[0]]*ressources[i[1]]*eff[i[1]] <= ke
             
        my_lp_problem.solve()
       
        x_res = np.zeros(n); j = 0;
        for i in indexes[0]:
            x_res[i] = (x[j].varValue); j = j + 1;
  
        return (pulp.value(my_lp_problem.objective), x_res)
    
# Maximise the net energy produced for a single technology
def maximiseNetEnergyTechnology(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax, wind, name):
    n = len(area)
    my_lp_problem = pulp.LpProblem("Net Energy Maximisation "+name+" Simple", pulp.LpMaximize)
    
    x = []
    for i in range (0, n):
        x.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Continuous'))
        if wind: my_lp_problem += x[i]*ressources[i]*eff[i] <= keMax[i]
            
    my_lp_problem += Calculation.netEnergySimpleModel(x, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y), "Z"
    my_lp_problem.solve()
    
    x_res = np.zeros(n); 
    for i in range (0, n):
        x_res[i] = (x[i].varValue); 
  
    return (pulp.value(my_lp_problem.objective), x_res)      

# Maximise the net energy produced on one cell: only 3 variables
def maximiseNetEnergyCell(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax):
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
        for i in indexes[2]:
            cons.append({'type': 'ineq', 'fun' : ke_constraint(i[0], ressources[i[1]], eff[i[1]], keMax)})
            
        res = minimize(obj, x0=np.ones(n), bounds=Calculation.binary_bounds(n), constraints=cons, method='trust-constr')
        x_res = np.zeros(3); j = 0;
        for i in indexes[0]:
            x_res[i] = res.x[j]; j = j + 1;
        
        return (-obj(res.x), x_res)
def ke_constraint(i, ressources, eff, ke):
    def g(x):
        return np.array([-x[i]*ressources*eff+ke])
    return g    
# Maximise the net energy produced for the whole grid: 3 variables / grid cell !
def maximiseNetEnergyGrid(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    
    indexes = Calculation.getIndexesSimpleModel(area)
    nX = len(indexes[0])
    print nX, " decision variables"
    def obj(x): 
        return -Calculation.netEnergySimpleModel(x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]])
    
    cons = []          
    for i in indexes[1]:
        cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(i[0], i[1])})
    for i in indexes[2]:
        cons.append({'type': 'ineq', 'fun' : ke_constraint(i[0], ressources[i[1]], eff[i[1]], keMax[i[1]/3])})
            
    res = minimize(obj, x0=np.ones(nX), bounds=Calculation.binary_bounds(nX), constraints=cons, method='trust-constr')
    x_res = np.zeros(n); j = 0;
    for i in indexes[0]:
        x_res[i] = res.x[j]; j = j + 1;
        
    return (-obj(res.x), x_res)

if __name__ == "__main__":
    sys.exit(main())
