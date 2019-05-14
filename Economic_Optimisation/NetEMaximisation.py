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
    inputs = 'inputs/params'
    results_maximiseNetEnergyCell(inputs, 'outputs/test_netE', False, 10)
           
def results_maximiseNetEnergyCell(opti_inputs, output_file, total, size):
    t0 = time.time()
    (lats, lon, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind) = Calculation.loadData(opti_inputs)
    if total: 
        n = len(lats) 
    else: 
        n = size
    print "# Cells:", n
    output = open(output_file, 'w')
    for i in range(0, n):
        if(i % 1000 == 0):
            print "Progress ", round(float(i) / float(n) * 100), "%"
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
        if sum(area[i]) > 0:
           res = maximiseNetEnergyCell(area[i], c[i], k[i], ghi[i], dni[i], embodiedE1y_wind[i],operationE_wind[i], avail_wind[i])
           output.write(str(res[0] / 1000) + "\t")
           output.write(str(res[1].x[0]) + "\t" + str(res[1].x[1]) + "\t" + str(res[1].x[2]) + "\t" + str(res[1].x[3])+ "\t" + str(res[1].x[4])+ "\t" + str(res[1].x[5]))
           # Write Cf as it is difficult to calculate afterwards
           output.write("\t" + str(Calculation.capacityFactor(c[i], k[i], res[1].x[3])))
        else:
           output.write("0.0" + "\t" + "0.0" + "\t" + "0.0" + "\t" + "0.0" + "\t" + "0.0")
        output.write("\n")
    output.close()
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"

# Maximise the net energy produced on one cell: 3 variables x_ij + vr_i + n_i + SM_i
def maximiseNetEnergyCell(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    def obj(x): 
        return -Calculation.netEnergy(x, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)
    
    cons = ({'type': 'ineq', 'fun' : lambda x:-x[1] - x[2] + 1})
    res = minimize(obj, x0=(1.0, 1.0, 1.0, 11.0, 10.0, 2.7), bounds=[(0, 1.0), (0, 1.0), (0, 1.0), (10.0, 16.0), (1.0, 20.0), (1.0, 4.0)], constraints=[cons], method='trust-constr')
    return (-obj(res.x), res)

if __name__ == "__main__":
    sys.exit(main())
