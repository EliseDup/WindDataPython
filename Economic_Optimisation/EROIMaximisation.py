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
    #results_maximiseEROICell(Calculation.inputs_params, 'test_eroi_cell', False, 0, 1000)
    results_maximiseEROIGrid(Calculation.inputs_params, 'test_eroi_grid', False, 0, 100)
      
def results_maximiseEROICell(opti_inputs, output_file, total, start, size):
    t0 = time.time()
    (lats, lon, total_area, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind, keMax) = Calculation.loadData(opti_inputs)
     
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
        res = maximiseEROIScalar(area[i], c[i], k[i], ghi[i], dni[i], embodiedE1y_wind[i], operationE_wind[i], avail_wind[i], keMax[i])
        total += res[0]
        Calculation.writeDetailedResultsCell(output, lats[i], lon[i], res, area[i], total_area[i])
        
    output.close()
    print "Results Grid: mean EROI = ", total / n
    print "Optimization cell per cell completed in ", (time.time() - t0), " seconds"
    
def results_maximiseEROIGrid(opti_inputs, output_file, total, start, size):
    t0 = time.time()
    (lats, lon, total_area, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind, keMax) = Calculation.loadData(opti_inputs)
    if total: 
        n = len(lats); start = 0;
    else: 
        n = size
    print "# Cells:", n
    res = maximiseEROI(area[start:n+start, :], c[start:n+start], k[start:n+start], ghi[start:n+start], dni[start:n+start], embodiedE1y_wind[start:n+start],operationE_wind[start:n+start], avail_wind[start:n+start], keMax[start:n+start])
    print "Results Grid: mean EROI = ", res[0]
    Calculation.writeDetailedResultsGrid(output_file, start, n, lats, lon, res, area, total_area, c, k, 6)
    print "Optimization for ",n,"cells in ", (time.time() - t0), " seconds"

def maximiseEROI(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind, keMax):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    # indexes = Calculation.getIndexes(area)
    area = area.flatten()
    def obj(x): 
        return -Calculation.eroi(x, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)
        
    cons = []; xstart = []; bnds=[];
    for i in range(0, len(c)):
        cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(i*6+1,i*6+2)})
        cons.append({'type': 'ineq', 'fun' : ke_constraint(i*6, c[i], k[i], avail_wind[i], keMax[i])})
        xstart.append(1.0);xstart.append(1.0);xstart.append(1.0);xstart.append(11.0);xstart.append(10.0);xstart.append(2.7);
        bnds.append((0, 1.0));bnds.append((0, 1.0));bnds.append((0, 1.0));
        bnds.append((10.0, 16.0)); bnds.append((1.0, 20.0));bnds.append((1.0, 4.0))
    
    res = minimize(obj, x0=xstart, bounds=bnds, constraints=cons, method='trust-constr')
        
    return (-obj(res.x), res.x)

def maximiseEROIScalar(area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind, keMax):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    # indexes = Calculation.getIndexes(area)
    def obj(x): 
        return -Calculation.eroiScalar(x, area, c, k, ghi, dni, embodiedE1y_wind, operationE_wind, avail_wind)
        
    cons = []; xstart = []; bnds=[];
    cons.append({'type': 'ineq', 'fun' : Calculation.non_complementary_constraint(1, 2)})
    cons.append({'type': 'ineq', 'fun' : ke_constraint(0, c, k, avail_wind, keMax)})
    xstart.append(1.0);xstart.append(1.0);xstart.append(1.0);xstart.append(11.0);xstart.append(10.0);xstart.append(2.7);
    bnds.append((0, 1.0));bnds.append((0, 1.0));bnds.append((0, 1.0));
    bnds.append((10.0, 16.0)); bnds.append((1.0, 20.0));bnds.append((1.0, 4.0))
    
    res = minimize(obj, x0=xstart, bounds=bnds, constraints=cons, method='trust-constr')
        
    return (-obj(res.x), np.append(res.x,Calculation.capacityFactor(c, k, res.x[3])))

def ke_constraint(i, c, k, avail_wind, keMax):
    def g(x):
        return np.array([- x[i] * Calculation.installedCapacityDensityWind(x[i+3], x[i+4]) * Calculation.capacityFactor(c, k, x[i+3]) * Calculation.arrayEfficiency(x[i+4]) * avail_wind + keMax])
    return g

if __name__ == "__main__":
    sys.exit(main())
