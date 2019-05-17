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

def main():
    inputs = 'inputs/params'
    results_maximiseConsumptionGrid(inputs, 'outputs/test_cons', False, 10)
       
def results_maximiseConsumptionGrid(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind) = Calculation.loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
        
    print "# Cells:", n
    output = open(output_file, 'w')
    res = maximiseConsumptionGrid(area[0:n, :], c[0:n], k[0:n], ghi[0:n], dni[0:n], embodiedE1y_wind[0:n], operationE_wind[0:n], avail_wind[0:n])
    print "Results Grid ", res[0] / 1E6, " TWh "
    Calculation.writeResultsGrid(output_file, lats, lon, res)
    print "Consumption Maximisation grid completed in ", (time.time() - t0), " seconds"
 
# Maximise the consumption on the whole grid
def maximiseConsumptionGrid(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind):
    n = 3 * len(c)
    def obj(x): 
        return -sum(x)
    
    cons = []
    for i in range(0, (n / 3)):
        cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(i * 3 + 1, i * 3 + 2)})
    res = minimize(obj, x0=np.ones(n), bounds=Calculation.binary_bounds(n), constraints=cons, method='trust-constr')
    return  (-obj(res.x), res.x)

if __name__ == "__main__":
    sys.exit(main())
