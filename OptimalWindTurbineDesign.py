import sys
import numpy as np
from numpy import genfromtxt
import math
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root
from datetime import datetime
from scipy.special import gammainc

def main():
    
    data = genfromtxt('../WindPotential/world_c_k', delimiter='\t', dtype=None)
    output = open('res_world_vr', 'w')
    lats = data[:, 0]; lon = data[:, 1]
    c = data[:, 2]; k = data[:, 3]; energyInputs = data[:, 4]
    n = len(lats)
    for i in range(0, n):
        if i%1000==0: print 'Progress :',i,'/',n
        res = optimalRatedSpeed(k[i], c[i], energyInputs[i])
        output.write(str(lats[i]) + '\t' + str(lon[i]) + '\t' +  str(res[0]) + '\t' + str(rotorDiameter(1E6, res[0])) + '\n')
    output.close()
         
def capacityFactor(k, c, vr):
    vf = 25.0; vc = 3.0;
    return -math.exp(-math.pow(vf/c,k)) + (3*math.pow(c,3)*math.gamma(3.0/k) / (k*(math.pow(vr,3) - math.pow(vc,3)))) * (gammainc(3.0/k,math.pow(vr/c,k)) - gammainc(3.0/k,math.pow(vc/c,k)))

# Rated power pr in MW, rated wind speed in m/s, returns the rotor diameter in m
def rotorDiameter(pr, vr, cp = 0.45):
    return math.sqrt(pr / (cp*0.5*1.225*math.pi/4*math.pow(vr,3) ) )
# Diameter in meters, rated wind speed in m/s, returns the corresponding rated power in MW
def ratedPower(d, vr, cp = 0.45):
    return (cp*0.5*1.225*math.pi*math.pow(d,2)/4*math.pow(vr,3))

# Determine performance coefficient from wind turbine specifications
def cp(rp, d, vr):
    return rp*1E6/(0.5*1.225*math.pi*d*d/4.0*vr*vr*vr)

def optimalRatedSpeed(k, c, energyInputs, n = 10, cp = 0.45):
    # Returns - power production density (We / m^2)
    def productionDensity(x):
        return -( math.pow(x,3) * capacityFactor(k, c, x) *cp*0.5*1.225*math.pi/4/(n*n) * 25*365*24 - energyInputs*1E6 / (n*n*rotorDiameter(1E6, x, cp)*rotorDiameter(1E6, x, cp)))

    res = minimize(productionDensity, x0=8, bounds=[(10.0, 16.0)]) #, options={'disp': True})
    return (res.x[0], -productionDensity(res.x[0]))
    
    
if __name__ == "__main__":
    sys.exit(main())
