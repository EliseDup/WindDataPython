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
def main():
    print "Hello"
    sm = ("0_9","1","1_1","1_2","1_3","1_4","1_5","1_6","1_7","1_8")
    #for s in sm:
    #    print s
    #    efficiency_dni("../resources/data_solar/sam_data/latitude_dni_eff_trough0h_sm"+s)
     
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
    interpolation(ptpptes[:,0],ptpptes[:,1],1)
    interpolation(ptpptes[:,0],ptpptes[:,2],1)
    
    plt.show()
    #stpptes = genfromtxt("../resources/data_solar/sam_data/STPP-TES", delimiter='\t', dtype=None)
    
    #interpolation(stpptes[:,0],stpptes[:,1],1)
    #interpolation(stpptes[:,0],stpptes[:,2],1)
   
    print "End"
    
    
def efficiency_dni(file):
    data_dni = genfromtxt(file, delimiter='\t', dtype=None)
    dni_eff = ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    print dni_eff
    return 

def efficiency_trough0h(sm):  
    data_dni = genfromtxt("../sam_data/latitude_dni_eff_trough0h_sm1_3", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough0h", delimiter='\t', dtype=None)
    dni_eff = ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    sm_eff = interpolation(data_sm[:, 0], data_sm[:, 1], 4)
    return 

def efficiency_trough12h(sm):  
    data_dni = genfromtxt("../sam_data/latitude_dni_eff_trough12h_sm2_7", delimiter='\t', dtype=None)
    data_sm = genfromtxt("../sam_data/sm_norm_efficiency_trough12h", delimiter='\t', dtype=None)
    dni_eff = ln_interpolation(data_dni[:, 1], data_dni[:, 2])
    sm_eff = interpolation(data_sm[:, 0], data_sm[:, 1], 4)
    return 

def ln_interpolation(x, y, plot=True):
    print "Interpolation with ", x.size, " measurements"
    logx = np.log(x)
    z = np.poly1d(np.polyfit(logx, y, 1))
    
    num = 0; den = 0; mean = sum(y) / y.size
    for i in range(0, x.size):
        xi = np.log(x[i]); 
        yi = y[i]
        num += (yi - z(xi)) * (yi - z(xi))
        den += (yi - mean) * (yi - mean)
    rsquare = 1 - num / den
    print 'R2 = ', rsquare
    
    if(plot):
        xp = np.linspace(min(x), max(x), 1000)
        plt.subplot(211) 
        #plt.plot(np.log(x), y, '.', np.log(xp), z(np.log(xp)), '-'); 
        # plt.xlabel("ln(DNI [kWh/m2/year])"); plt.ylabel("Efficiency [%]"); # plt.draw();
        plt.plot(x, y, '.', xp, z(np.log(xp)), '-'); 
        plt.xlabel("DNI [kWh/m2/year]"); plt.ylabel("Efficiency [%]");  # plt.draw();
        
        plt.subplot(212) 
        plt.plot(np.log(x), y - z(np.log(x)),'.'); 
        
        plt.xlabel('R2='+str(rsquare)); plt.ylabel("Residuals");
        plt.draw();
    
    return z

def sqrt_interpolation(x, y, plot=True):
    print "Interpolation with ", x.size, " measurements"
    sqrtx = np.sqrt(x)
    z = np.poly1d(np.polyfit(sqrtx, y, 1))
    
    num = 0; den = 0; mean = sum(y) / y.size
    for i in range(0, x.size):
        xi = np.sqrt(x[i]); 
        yi = y[i]
        num += (yi - z(xi)) * (yi - z(xi))
        den += (yi - mean) * (yi - mean)
    rsquare = 1 - num / den
    print 'R2 = ', rsquare
    
    if(plot):
        xp = np.linspace(min(x), max(x), 1000)
        # plt.subplot(211) 
        plt.plot(np.sqrt(x), y, '.', np.sqrt(xp), z(np.sqrt(xp)), '-'); 
        plt.xlabel("ln(DNI [kWh/m2/year])"); plt.ylabel("Efficiency [%]");  # plt.draw();
        # plt.subplot(212) 
        # plt.plot(np.log(x), y - z(np.log(x)),'.'); 
        
        # plt.xlabel('R2='+str(rsquare)); plt.ylabel("Residuals");
        plt.draw();
    
    return z


def interpolation(x, y, degree, plot=True):
    z = np.poly1d(np.polyfit(x, y, degree))
    print z
    num = 0; den = 0; mean = sum(y) / y.size
    for i in range(0, x.size):
        xi = x[i] 
        yi = y[i]
        num += (yi - z(xi)) * (yi - z(xi))
        den += (yi - mean) * (yi - mean)
    rsquare = 1 - num / den
    print 'R2 = ', rsquare
    
    if (plot):
        xp = np.linspace(min(x), max(x), 1000)
        # plt.subplot(211) 
        plt.plot(x, y, '.', xp, z(xp)), '-'; 
        plt.xlabel("DNI [kWh/m2/year]"); plt.ylabel("Efficiency [%]"); plt.draw();
        # plt.subplot(212) 
        # plt.plot(x, y - z(x),'.'); 
        
        # plt.xlabel('R2='+str(rsquare)); plt.ylabel("Residuals");
        # plt.draw();
    
    return z

if __name__ == "__main__":
    sys.exit(main())
 
