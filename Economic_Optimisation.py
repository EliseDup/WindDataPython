import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot

def main():
    maximiseNetEnergyGrid('../WindPotentialScala/simple')
    print "Hello"
    
def maximiseNetEnergyGrid(opti_inputs):
     data = genfromtxt(opti_inputs, delimiter='\t', dtype=None)
     
     lats = data[:, 0]; lon = data[:, 1]
     # Potential suitable area for each tecnology
     area = data[:, 2:5]
     eff = data[:, 5:8]
     ressources = data[:, 8:11]
     installed_capaciy_density = data[:, 11:14]
     operationE = data[:, 14:17]
     embodiedE = data[:, 17:20]
     
     n = len(lats)

     print "Size problem ", n
     
     output = open('solar', 'w')
     for i in range(0, n):
        output.write(str(lats[i]) + "\t" + str(lon[i]) + "\t")
         
        output.write("\n")
     output.close()

def maximiseNetEnergyCell(area, eff, ressources, installed_capaciy_density, operationE, embodiedE):
    
    #def netEnergy(x):
       
    # Constraint : x[0] + x[1] <= 1
    def constraint(x):
        return -x[0] - x[1] + 1.0
    cons = {'type':'ineq', 'fun': constraint}
    res = minimize(net_energy, x0=(1.0, 1.0, 1.0), bounds=[(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)], constraints=cons)
    
if __name__ == "__main__":
    sys.exit(main())
