import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
import scipy as sp
import scipy.interpolate as interp
from Tkconstants import BOTTOM
from scipy.interpolate.interpolate import interp1d

def main():
    print "Hello"
    efficiency_dni(2000, "../sam_data/latitude_dni_eff_tower12h_sm2_7")
    efficiency_dni(2000, "../sam_data/latitude_dni_eff_trough0h_sm1_615")
    
    #efficiency_dni(2000, "../sam_data/latitude_dni_eff_trough12h_sm2_7")
    #efficiency_dni(2000, "../sam_data/latitude_dni_eff_trough0h_sm1_3")

    print "End"
    plt.show()
    
def efficiency_dni(dni, file):
    data_dni = genfromtxt(file, delimiter='\t', dtype=None)
    dni_eff = ln_interpolation(data_dni[:, 1], data_dni[:, 2], 1)
    print dni_eff
    return dni_eff(np.log(dni))

def efficiency_trough0h(dni, sm):  
    data_dni = genfromtxt("../sam_data/latitude_dni_eff_trough0h_sm1_3", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough0h", delimiter='\t', dtype=None)
    dni_eff = ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    sm_eff = interpolation(data_sm[:, 0], data_sm[:, 1], 4)
    return dni_eff(np.log(dni)) * sm_eff(sm) / 100.0

def efficiency_trough12h(dni, sm):  
    data_dni = genfromtxt("../sam_data/latitude_dni_eff_trough12h_sm2_7", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough12h", delimiter='\t', dtype=None)
    dni_eff = ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    sm_eff = interpolation(data_sm[:, 0], data_sm[:, 1], 4)
    return dni_eff(np.log(dni)) * sm_eff(sm) / 100.0

def ln_interpolation(x, y, plot=False):
    logx = np.log(x)
    z = np.poly1d(np.polyfit(logx, y, 1))
    if(plot):
        xp = np.linspace(min(x), max(x), 1000)
        plt.figure; plt.plot(x, y, '.', xp, z(np.log(xp)), '-'); 
        plt.xlabel("DNI [kWh/m2/year]"); plt.ylabel("Efficiency [%]"); plt.draw();
    return z

def interpolation(x, y, degree=1, plot=False):
    z = np.poly1d(np.polyfit(x, y, degree))
    if (plot):
        xp = np.linspace(min(x), max(x), 1000)
        plt.figure; plt.plot(x, y, '.', xp, z(xp)); 
        plt.xlabel("DNI [kWh/m2/year]"); plt.ylabel("Efficiency [%]"); plt.draw();
  
    return z

if __name__ == "__main__":
    sys.exit(main())
 
