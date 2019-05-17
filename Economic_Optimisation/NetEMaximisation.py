import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot
import Calculation

# from scipy.optimize import LinearConstraint
import time

def main():
    results_maximiseNetEnergyCell(Calculation.inputs_params, 'outputs/test_netE', False, 100)
           
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind) = Calculation.loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
    print "# Cells:", n
    total = 0
    output = open(output_file, 'w')
    for i in range(0, n):
        if(i % (n/10) == 0):
            print "Progress ", round(float(i) / float(n) * 100), "%"
        res = maximiseNetEnergyCell(area[i], c[i], k[i], ghi[i], dni[i], embodiedE1y_wind[i],operationE_wind[i], avail_wind[i])
        total += res[0]
        Calculation.writeResultsCell(output, lats[i], lon[i], res)
    output.close()
    print "Results Grid ", total / 1E6, " TWh "
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"

# Maximise the net energy produced on one cell: 3 variables x_ij + vr_i + n_i + SM_i
def maximiseNetEnergyCell(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind):
    if(sum(area) < 1E-3):
        return(0.0, np.zeros(7))
    else:
        # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
        def obj(x): 
            return -Calculation.netEnergy(x, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)
        cons = ({'type': 'ineq', 'fun' : lambda x:-x[1] - x[2] + 1})
        res = minimize(obj, x0=(1.0, 1.0, 1.0, 11.0, 10.0, 2.7), bounds=[(0, 1.0), (0, 1.0), (0, 1.0), (10.0, 16.0), (1.0, 20.0), (1.0, 4.0)], constraints=[cons], method='trust-constr')
        return (-obj(res.x), np.append(res.x,Calculation.capacityFactor(c, k, res.x[3])))


# To do -> Not possible with Pulp  as the objective is non linear ?
def maximiseNetEnergy(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind):
    area = area.flatten();
    (ind_x, ind_wind, ind_csp, ind_cons) = Calculation.getIndexes(area)
    
    my_lp_problem = pulp.LpProblem("NE Maximisation", pulp.LpMaximize)
    
    return

if __name__ == "__main__":
    sys.exit(main())
