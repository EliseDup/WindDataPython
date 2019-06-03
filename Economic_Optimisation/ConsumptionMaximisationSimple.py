import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot
# from scipy.optimize import LinearConstraint
import time
import Calculation
import pulp

# Parameters of the economy
# Y = E / qF = K / vF [$]
# C = Y - delta_E K_E - delta_F K_F = (1 - delta_F) / vF * E - delta_E / q_F * \tilde{K}_E
# E is the energy available for the economy ( = Gross E - Operation E)
# \tilde{K}_E is the energy embodied in the energy sector capital stock
# delta are depreciation rates
deltaE = 1.0/25
deltaF = 1.0/20
# Estimated based on statistical data : PIB and secondary energy supply [MWh/US$]
qF = 0.002
# Considering the value of 4.3 from Picketty as an upper bound for the world capital ratio
vF = 4

# Consumption in Dollars !
def main():
    res = results_maximiseConsumptionGrid('inputs_simple_sf', 'consumptionMax', True, 0, 10000, qF, vF)
    
# Start is the first index used to compute results, size is the size of the cells where the optimization is completed
def results_maximiseConsumptionGrid(opti_inputs, output_file, total, start, size, qF = qF, vF = vF):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE, keMax) =  Calculation.loadDataSimpleModel(opti_inputs)
    if total: 
        n = len(lats) 
        start = 0
    else: 
        n = size
          
    print "# Cells:", n
    output = open(output_file, 'w')
    res = maximiseConsumptionGrid(qF,vF,start,area[start:n+start, :], eff[start:n+start, :], ressources[start:n+start, :], installed_capaciy_density[start:n+start, :], operationE[start:n+start, :], embodiedE1y[start:n+start, :], keMax[start:n+start])
    print res[0] / 1E9
    print "Total consumption ", round(res[0] / 1E9), "GUS$ ", round((res[0]/1E9/80935)*100), "% of 2017 data"
    print "Net Energy Grid ", round(res[2]/1E6), " TWh "
    Calculation.writeResultsGrid(output_file, start, n, lats, lon, res)
    print "Consumption Maximisation grid completed in ", (time.time() - t0), " seconds" 
    
    return (res[0]/1E6, res[2]/1E6)

# Maximise the consumption on the whole grid
def maximiseConsumptionGrid(qF,vF, start, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y, keMax):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    
    if(sum(area) < 0.001):
        return(0.0, np.zeros(n), 0.0)
    else:
        my_lp_problem = pulp.LpProblem("Consumption Maximisation", pulp.LpMaximize)
        indexes = Calculation.getIndexesSimpleModel(area)
        
        print len(indexes[0]), " decision variables"
        
        x = []
        for i in indexes[0]:
            x.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Continuous'))
        # Objective function
        my_lp_problem += Calculation.consumptionSimple2(qF, vF, deltaE, deltaF, x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]]), "Z"
        
        # Constraints
        for i in indexes[1]:
            my_lp_problem += x[i[0]] + x[i[1]] <= 1
        for i in indexes[2]:
            my_lp_problem += x[i[0]]*area[i[1]]*ressources[i[1]]*eff[i[1]] <= keMax[i[1]/3]
            
        my_lp_problem.solve()
        
        x_res = np.zeros(n); j = 0;
        for i in indexes[0]:
            x_res[i] = (x[j].varValue); j = j + 1;
        
        netE = Calculation.netEnergySimpleModel(x_res, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y)
        return (pulp.value(my_lp_problem.objective), x_res, netE)

if __name__ == "__main__":
    sys.exit(main())
