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
    results_maximiseEROIGrid('inputs/inputs_simple_sf', 'outputs/test_EROI', False, 0, 100)
      
def results_maximiseEROICell(opti_inputs, output_file, total, start, size):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModel(opti_inputs)
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
        res = maximiseEROI(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i], keMax[i])
        total += res[0]
        Calculation.writeResultsCell(output, lats[i], lon[i], res)
    output.close()
    print "Results Grid ", total / 1E6, " TWh "
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"
    
def results_maximiseEROIGrid(opti_inputs, output_file, total, start, size):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) = Calculation.loadDataSimpleModel(opti_inputs)
    if total: 
        n = len(lats); start = 0;
    else: 
        n = size
    print "# Cells:", n
    res = maximiseEROI(area[start:n+start, :], eff[start:n+start, :], ressources[start:n+start, :], installed_capaciy_density[start:n+start, :], operationE[start:n+start, :], embodiedE1y[start:n+start, :], keMax[start:n+start])
    print "Results Grid ", res[0] / 1E6, " TWh "
    Calculation.writeResultsGrid(output_file, start, n, lats, lon, res)
    print "Optimization for ",n,"cells in ", (time.time() - t0), " seconds"

def maximiseEROI(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    
    if(sum(area) == 0):
        return(0.0, np.zeros(n))
    else:
        my_lp_problem = pulp.LpProblem("EROI Maximisation", pulp.LpMaximize)
        indexes = Calculation.getIndexesSimpleModel(area)
        nX = len(indexes[0])
        if(nX > 3):
            print nX, " decision variables"
        x = []
        for i in indexes[0]:
            x.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Continuous'))
        # Objective function
        my_lp_problem += Calculation.eroi(x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]]), "Z"
        # Constraints: non complementary constraints for PV and CSP
        for i in indexes[1]:
            my_lp_problem += x[i[0]] + x[i[1]] <= 1
        # Constraints: max KE generation for wind 
        for i in indexes[2]:
            if len(area) == 3:
                ke = keMax
            else:
                ke = keMax[i[1]/3]
            my_lp_problem += x[i[0]]*area[i[1]]*ressources[i[1]]*eff[i[1]] <= ke
             
        my_lp_problem.solve()
       
        x_res = np.zeros(n); j = 0;
        for i in indexes[0]:
            x_res[i] = (x[j].varValue); j = j + 1;
           
        return (pulp.value(my_lp_problem.objective), x_res)

def ke_constraint(i, area, ressources, eff, ke):
    def g(x):
        return np.array([-x[i]*area*ressources*eff+ke])
    return g    

if __name__ == "__main__":
    sys.exit(main())
