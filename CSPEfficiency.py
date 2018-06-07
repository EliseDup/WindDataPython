import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt
import scipy as sp
import scipy.interpolate

def main():
    
    sm = np.linspace(0.1, 5, 100)
    plt.figure; plt.plot(sm, efficiency_trough0h(2000,sm)); plt.show()
#     data_dni = genfromtxt("../sam_data/dni_eff_trough0h", delimiter='\t', dtype=None)
#     data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough0h", delimiter='\t', dtype=None)
#     dni_eff = ln_interp1d(data_dni[:,0], data_dni[:,1],  np.linspace(500, 5000, 100),False)
#     sm_eff = interp1d(data_sm[:,0], data_sm[:,1],  np.linspace(0, 5, 100), 5,False)
#     x = np.linspace(0.1, 5, 100)
#     y1000 = dni_eff(np.log(1000))*sm_eff(x)
#     y2000 = dni_eff(np.log(2000))*sm_eff(x)
#     y3000 = dni_eff(np.log(3000))*sm_eff(x)
#     plt.figure; plt.plot(x,y1000,x,y2000,x,y3000); plt.show()
    
    
def efficiency_trough0h(dni, sm):  
    data_dni = genfromtxt("../sam_data/dni_eff_trough0h", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough0h", delimiter='\t', dtype=None)
    dni_eff = np.poly1d(np.polyfit(np.log(data_dni[:,0]), data_dni[:,1], 1))
    sm_eff = np.poly1d(np.polyfit(data_sm[:,0], data_sm[:,1],4))
    return dni_eff(np.log(dni))*sm_eff(sm)/100.0

def efficiency_trough12h(dni, sm):  
    data_dni = genfromtxt("../sam_data/dni_eff_trough12h", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough12h", delimiter='\t', dtype=None)
    dni_eff = np.poly1d(np.polyfit(np.log(data_dni[:,0]), data_dni[:,1], 1))
    sm_eff = np.poly1d(np.polyfit(data_sm[:,0], data_sm[:,1],4))
    return dni_eff(np.log(dni))*sm_eff(sm)/100.0

def ln_interp1d(x, y, xp, plot):
    logx = np.log(x)
    z = np.poly1d(np.polyfit(logx, y, 1))
    print 'interpolation a ln(x) + b:', z
    if (plot):
        plt.figure; plt.plot(x,y,'.',xp,z(np.log(xp))); plt.show()
    return z

def interp1d(x, y, xp, degree, plot):
    z = np.poly1d(np.polyfit(x, y, degree))
    print 'interpolation', degree, ':', z
    if (plot):
        plt.figure; plt.plot(x,y,'.',xp,z(xp)); plt.show()
    return z

if __name__ == "__main__":
    sys.exit(main())
 