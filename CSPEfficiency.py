import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import genfromtxt
import scipy as sp
import scipy.interpolate as interp
from Tkconstants import BOTTOM
from scipy.interpolate.interpolate import interp1d
import math
import Interpolation

def main():
    print "Hello"
    #efficiency_dni("../resources/data_solar/sam_data/latitude_dni_eff_trough0h_sm2") 
    #sm = ("2_7","3","3_3","3_6")
    #for s in sm:
    #    print s
    #    efficiency_dni("../resources/data_solar/sam_data/latitude_dni_eff_tower12h_sm"+s)
     
#===============================================================================
#     fig = plt.figure()
#     ax = plt.axes(projection='3d')
#     def f(x, y):
#         return (-1.585*y + 11.19)*np.log(x*8.76) + (10.58*y - 66.09)
# 
#     x = np.linspace(50, 500, 1000)
#     y = np.linspace(2, 4, 1000)
# 
#     X, Y = np.meshgrid(x, y)
#     Z = f(X, Y)
#    
#     fig = plt.figure()
#     ax = plt.axes(projection='3d')
#     ax.contour3D(Y, X, Z, 1000, rstride=1, cstride=1, cmap='viridis', edgecolor='none');
#     ax.set_ylabel('DNI [W/m2]')
#     ax.set_xlabel('SM')
#     ax.set_zlabel('Efficiency [%]');
#     
#===============================================================================
    #
  
    #ptpp = genfromtxt("../resources/data_solar/sam_data/PTPP", delimiter='\t', dtype=None)
    #interpolation(ptpp[:,0],ptpp[:,1],1)
    #interpolation(ptpp[:,0],ptpp[:,2],1)
    
    ptpptes = genfromtxt("../resources/data_solar/sam_data/PTPP-TES", delimiter='\t', dtype=None)
    Interpolation.interpolation(ptpptes[:,0],ptpptes[:,1],1)
    Interpolation.interpolation(ptpptes[:,0],ptpptes[:,2],1)
    
    
    #stpptes = genfromtxt("../resources/data_solar/sam_data/STPP-TES", delimiter='\t', dtype=None)
    #interpolation(stpptes[:,0],stpptes[:,1],1)
    #interpolation(stpptes[:,0],stpptes[:,2],1)
    plt.show()
    print "End"
    
    
def efficiency_dni(file):
    data_dni = genfromtxt(file, delimiter='\t', dtype=None)
    dni_eff = Interpolation.ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    print dni_eff
    return 

def efficiency_trough0h(sm):  
    data_dni = genfromtxt("../sam_data/latitude_dni_eff_trough0h_sm1_3", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough0h", delimiter='\t', dtype=None)
    dni_eff = Interpolation.ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    sm_eff = interpolation(data_sm[:, 0], data_sm[:, 1], 4)
    return 

def efficiency_trough12h(sm):  
    data_dni = genfromtxt("../sam_data/latitude_dni_eff_trough12h_sm2_7", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough12h", delimiter='\t', dtype=None)
    dni_eff = Interpolation.ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    sm_eff = Interpolation.interpolation(data_sm[:, 0], data_sm[:, 1], 4,"DNI [kWh/m2/year]","Efficiency [%]")
    return 

if __name__ == "__main__":
    sys.exit(main())
 
