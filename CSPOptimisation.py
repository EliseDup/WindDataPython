import sys
import numpy as np
from scipy.optimize import minimize, root

def main():
    dni = 400*8.76
    print optimizeSolarMultiple(dni, 2922368.611, 4487437.377, 0.1, default_area(1))
    print "Hello"

def default_area(powerMW):
    return powerMW*1E6/(950*0.15)

def fullLoadHours(dni, sm):
    return (2.5715*dni-694)*(-0.0371*sm*sm+0.4171*sm-0.0744)

def efficiency(sm):
    return (-0.8522*sm*sm*sm*sm+8.115*sm*sm*sm-27.56*sm*sm+36.01*sm-1.9)/100.0

def efficiency(dni):
    return (4.318e-10*dni*dni*dni - 5.991e-06*dni*dni + 0.01982*dni - 4.291)/100.0
# Calculate Solar Multiple that maximizes the EROI
# With : 
# - DNI : direct normal irradiance in kWh/m2/year
# - eff : solar to electricity efficiency (15%)
# - ee_fixed : fixed part of the embodied energy (linked to the power block [kWh]
# - ee_area : part of the embodied energy varying with the aperture area [kWh]
# - default_area : default aperture area [m2]
     
def optimizeSolarMultiple(dni, ee_fixed, ee_area, ee_prod, default_area):
    
    def eroi(x):
        prod = 30.0 * dni * efficiency(dni) * x * default_area 
        return - prod / (ee_fixed + ee_area * x + prod * ee_prod)
        #return - 30.0 * fullLoadHours(dni, x) * 1E3 / (ee_fixed + ee_area * x)
    res = minimize(eroi, x0=8, bounds=[(1.0, 4.0)]) #, options={'disp': True})
    return (res.x[0], -eroi(res.x[0]))
 
    
if __name__ == "__main__":
    sys.exit(main())
 