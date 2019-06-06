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
    techs = ("Wind-onshore","Wind-offshore","mono-Si-PV","ST-salt-TES")
    for t in techs:
        data = genfromtxt("../WindPotentialScala/"+t, delimiter='\t', dtype=None)
        interpolation(data[:,0],data[:,1],2,xname="Embodied Energy [EJ/year]",yname=t+"Production [EJ/year]")
    
    plt.show()
    print "---"
    
def interpolation(x, y, degree, plot=True, plotResidual=False, xname="",yname=""):
    print "Interpolation with ", x.size, " measurements"
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
        plt.figure()
        xp = np.linspace(min(x), max(x), 1000)
        if plotResidual: plt.subplot(211) 
        plt.plot(x, y, '.', xp, z(xp)), '-'; 
        plt.xlabel(xname); plt.ylabel(yname); plt.draw();
        if plotResidual:
            plt.subplot(212) 
            plt.plot(x, y - z(x),'.'); 
            plt.xlabel('R2='+str(rsquare)); plt.ylabel("Residuals");
        plt.draw();
    
    return z

def ln_interpolation(x, y, plot=True, xname="",yname=""):
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
        plt.figure
        xp = np.linspace(min(x), max(x), 1000)
        plt.subplot(211) 
        
        plt.plot(x, y, '.', xp, z(np.log(xp)), '-'); 
        plt.xlabel(xname); plt.ylabel(yname);
        
        plt.subplot(212) 
        plt.plot(np.log(x), y - z(np.log(x)),'.'); 
        
        plt.xlabel('R2='+str(rsquare)); plt.ylabel("Residuals");
        plt.draw();
    
    return z

def sqrt_interpolation(x, y, plot=True, xname="",yname=""):
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
        plt.plot(np.sqrt(x), y, '.', np.sqrt(xp), z(np.sqrt(xp)), '-'); 
        plt.xlabel(xname); plt.ylabel(yname);
        plt.draw();
    
    return z

if __name__ == "__main__":
    sys.exit(main())
