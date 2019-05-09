import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot
#from scipy.optimize import LinearConstraint

def main():
    maximiseNetEnergy('../WindPotentialScala/simple')
    print "Hello"
    
def maximiseNetEnergy(opti_inputs):
     import scipy
     print scipy.__version__
    
     data = genfromtxt(opti_inputs, delimiter='\t', dtype=None)
     
     lats = data[:, 0]; lon = data[:, 1]
     # Potential suitable area for each tecnology
     area = data[:, 2:5]
     eff = data[:, 5:8]
     ressources = data[:, 8:11]
     installed_capaciy_density = data[:, 11:14]
     embodiedE1y = data[:, 14:17]
     operationE = data[:, 17:20]
     n = len(lats)
     print "Size problem ", n

     m = 200
     res = maximiseNetEnergyGrid(area[0:m,:], eff[0:m,:], ressources[0:m,:], installed_capaciy_density[0:m,:], operationE[0:m,:], embodiedE1y[0:m,:])
     print res[1].x
     output = open('test', 'w')
     for i in range(0, m):
       output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
       resIndex = i*3
       # res = maximiseNetEnergyCell(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
       #output.write(str(res[i]) + "\t")
       output.write(str(area[i][0]) + "\t" + str(area[i][1]) + "\t" +str(area[i][2]))
       
       output.write(str(res[1].x[resIndex]) + "\t" + str(res[1].x[resIndex+1]) + "\t" + str(res[1].x[resIndex+2]))
       output.write("\n")
      
     output.close()
#      output = open('test', 'w')
#      for i in range(1, 10):
#         output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
#         res = maximiseNetEnergyCell(area[i], eff[i], ressources[i], installed_capaciy_density[i], operationE[i], embodiedE1y[i])
#         output.write(str(res[0]) + "\t")
#         output.write(str(res[1].x[0]) + "\t" + str(res[1].x[1]) + "\t" + str(res[1].x[2]))
#         output.write("\n")
#      
#      output.close()

# Maximise the net energy produced on one cell
def maximiseNetEnergyCell(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    def netEnergy(x): 
        production = x * area * eff * ressources * (365 * 24) * (np.ones(3) - operationE)
        ee = x * area * installed_capaciy_density * embodiedE1y 
        return -sum(production - ee)
    
    #linear_constraint = LinearConstraint([[0, 1, 1]], [0], [1])
    cons = ({'type': 'ineq', 'fun' : lambda x: np.array([-x[1]-x[2]+1])})
    res = minimize(netEnergy, x0=(1.0, 1.0, 1.0), bounds=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)], constraints=[cons], method='trust-constr')
    return (-netEnergy(res.x), res)

# Maximise the net energy produced for the whole grid
def maximiseNetEnergyGrid(area, eff, ressources, installed_capaciy_density, operationE, embodiedE1y):
    area = area.flatten()
    eff = eff.flatten()
    ressources = ressources.flatten()
    installed_capaciy_density = installed_capaciy_density.flatten()
    operationE = operationE.flatten()
    embodiedE1y = embodiedE1y.flatten()
    
    n = len(area)
    
    # Net Energy = energy produced [MWh/an] - embodied energy in installed capacity
    def netEnergy(x): 
        production = x * area * eff * ressources * (365 * 24) * (np.ones(n) - operationE)
        ee = x * area * installed_capaciy_density * embodiedE1y 
        return -sum(production - ee)

    x0 = np.ones(n)
    
    bnds = []
    for i in range(0, n):
        bnds.append((0.0,1.0))

    # In each cell pv and csp are not compatible = x[1]+x[2] <= 1   
    cons = []
    for i in range(0, (n/3)):
        cons.append({'type': 'ineq', 'fun' : lambda x: np.array([-x[i*3+1]-x[i*3+2]+1])})
   
    res = minimize(netEnergy,x0, bounds=bnds, constraints=cons, method='trust-constr')
    return (-netEnergy(res.x), res)

if __name__ == "__main__":
    sys.exit(main())
