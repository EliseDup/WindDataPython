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
# Y = E / qF = K / vF
# C = Y - delta_E K_E - delta_F K_F = (1 - delta_F) / vF * E - delta_E / q_F * \tilde{K}_E
# E is the energy available for the economy ( = Gross E - Operation E)
# \tilde{K}_E is the energy embodied in the energy sector capital stock
# delta are depreciation rates
deltaE = 1.0/25
deltaF = 1.0/20
qF = 1.4
vF = 1.4

def main():
    results_maximiseConsumptionGrid(Calculation.inputs_simple, 'outputs/test_cons', False, 1000)
       
def results_maximiseConsumptionGrid(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, eff, ressources, installed_capaciy_density, embodiedE1y, operationE) =  Calculation.loadDataSimpleModel(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
        
    print "# Cells:", n
    output = open(output_file, 'w')
    res = maximiseConsumptionGrid(area[0:n, :], eff[0:n, :], ressources[0:n, :], installed_capaciy_density[0:n, :], operationE[0:n, :], embodiedE1y[0:n, :])
    print "Total consumption ", res[0] / 1E6, " TWh "
    
    Calculation.writeResultsGrid(output_file, n, lats, lon, res)
    print "Consumption Maximisation grid completed in ", (time.time() - t0), " seconds"
 
# Maximise the consumption on the whole grid
def maximiseConsumptionGrid(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    area = area.flatten();eff = eff.flatten();ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten();operationE = operationE.flatten();embodiedE1y = embodiedE1y.flatten()
    n = len(area)
    
    if(sum(area) < 0.001):
        return(0.0, np.zeros(n))
    else:
        my_lp_problem = pulp.LpProblem("NE Maximisation Consumption", pulp.LpMaximize)
        indexes = Calculation.getIndexesSimpleModel(area)
        x = []
        for i in indexes[0]:
            x.append(pulp.LpVariable('x' + str(i), lowBound=0, upBound=1, cat='Continuous'))
        # Objective function
        my_lp_problem += Calculation.consumptionSimple(qF, vF, deltaE, deltaF, x, area[indexes[0]], eff[indexes[0]], ressources[indexes[0]], installed_capaciy_density[indexes[0]], operationE[indexes[0]], embodiedE1y[indexes[0]]), "Z"
        # Constraints
        for i in indexes[1]:
            my_lp_problem += x[i[0]] + x[i[1]] <= 1
        my_lp_problem.solve()
        x_res = np.zeros(n); j = 0;
        for i in indexes[0]:
            x_res[i] = (my_lp_problem.variables()[j].varValue); j = j + 1;
        print "Net Energy Grid ", Calculation.netEnergySimpleModel(x_res, area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y) / 1E6, " TWh "
        return (pulp.value(my_lp_problem.objective), x_res)

if __name__ == "__main__":
    sys.exit(main())
