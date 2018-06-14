import sys
import numpy as np
from scipy.optimize import minimize, root
from CSPEfficiency import efficiency_trough0h, efficiency_trough12h
from numpy import genfromtxt
import matplotlib.pyplot as plt

# Calculate Solar Multiple that maximizes the EROI
# With : 
# - DNI : direct normal irradiance in kWh/m2/year
# - eff : solar to electricity efficiency (15%)
# - ee_fixed : fixed part of the embodied energy (linked to the power block [kWh]
# - ee_area : part of the embodied energy varying with the aperture area [kWh]
# - default_area : default aperture area [m2]

# For the 3 technologies tested :
#EE_Fixed[kWh/MW]    EE_Variable[kWh/area for 1MW]    EE_output[GJ/GJ]
#2921710.2777777775    2664706.018702423    0.1
#5455141.666666667    1193483.7107424117    0.07300000000000001
#7208229.166666666    1996433.3375287117    0.07300000000000001
#
_12h = False

ee_fixed_trough0h = 2921710.2777777775   
ee_area_trough0h = 2664706.018702423
ee_output_trough0h = 0.1
ee_fixed_trough12h = 5455141.666666667  
ee_area_trough12h = 1193483.7107424117
ee_output_trough12h = 0.073
ee_fixed_power12h = 7208229.166666666
ee_area_power12h = 1996433.3375287117
ee_output_power12h = 0.073

def efficiency(dni, sm):
    if _12h: return efficiency_trough12h(dni,sm)
    else: return efficiency_trough0h(dni,sm)

def eroi(dni, sm):
    eff = efficiency(dni,sm)
    prod = 30 * dni * eff * sm * default_area(1)
    if _12h:
        ee_fixed = ee_fixed_trough12h; ee_area = ee_area_trough12h; ee_prod = ee_output_trough12h;
    else:
        ee_fixed = ee_fixed_trough0h; ee_area = ee_area_trough0h; ee_prod = ee_output_trough0h;
    return prod / (ee_fixed + ee_area * sm + prod * ee_prod)

def default_area(powerMW):
    return powerMW*1E6/(950*0.22)

def main():
#     x = np.linspace(0.1, 6, 100)
#     plt.plot(x, eroi(1500,x), label='1500 kWh/m2/y')
#     plt.plot(x, eroi(2000,x), label='2000 kWh/m2/y')
#     plt.plot(x, eroi(2500,x), label='2500 kWh/m2/y')
#     plt.plot(x, eroi(3000,x), label='3000 kWh/m2/y')
#     plt.plot(x, eroi(3500,x), label='3500 kWh/m2/y')
#    
#     plt.legend()
#     plt.show()
    
    data = genfromtxt('../WindPotentialScala/data_csp_optimisation', delimiter='\t', dtype=None)
    if(_12h): output = open('res_through_12h', 'w')
    else: output = open('res_through_0h', 'w')
    lats = data[:, 0]; lon = data[:, 1]; dni = data[:, 2];
   
    n = len(lats)
    for i in range(0, n):
        if i%1000==0: print 'Progress :',i,'/',n
        res = maximiseEROI(dni[i])
        output.write(str(lats[i]) + '\t' + str(lon[i]) + '\t' +  str(dni[i]) + '\t'+ str(res[0]) + '\t' + str(res[1]) + '\t' + str(efficiency_trough12h(dni[i],res[0])) + '\n')
     
    output.close()
     
    print "End"

def maximiseEROI(dni):
    def inv_eroi(x):
        return -eroi(dni, x)
    res = minimize(inv_eroi, x0=1, bounds=[(1, 6)]) #, options={'disp': True})
    
    return (res.x[0], -inv_eroi(res.x[0]))

def maximiseEfficiency(dni, ee_fixed, ee_area, ee_prod, default_area):
    def inv_eff(x):
        return -efficiency(dni,x)
    
    res = minimize(inv_eff, x0=1, bounds=[(1, 6)]) #, options={'disp': True})
    return (res.x[0], -eroi(res.x[0]))

if __name__ == "__main__":
    sys.exit(main())
 
# def fullLoadHours(dni, sm):
#     return (2.5715*dni-694)*(-0.0371*sm*sm+0.4171*sm-0.0744)
# 
# def efficiency(sm):
#     return (-0.8522*sm*sm*sm*sm+8.115*sm*sm*sm-27.56*sm*sm+36.01*sm-1.9)/100.0
# 
# def efficiency(dni):
#     return (4.318e-10*dni*dni*dni - 5.991e-06*dni*dni + 0.01982*dni - 4.291)/100.0