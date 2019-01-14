from numpy import*
import matplotlib.pyplot as plt

t0 = 1800
t = t0
nYears =300

# Capital stock of the three sectors [EJ].
# We make the assumption that at its incept date, each energy sub-sector begins 
# with an installed capital of 0.1 [EJ].  This value is arbitrary and can of
# course be changed.  We assume that the main economy begins with 0 capital
# stock in 1800.
Cap1 = 0.1 * ones((6,nYears))
Cap2 = 0.1 * ones((8,nYears))
Cap3 = zeros(nYears)

# Annual energy production of the two energy sectors [EJ/yr]
P1 = zeros((6,nYears))
P2 = zeros((8,nYears))

# Total annual energy production [EJ/yr]
P = zeros(nYears)

# Net annual energy yield of the two energy sectors [EJ/yr]
Net1 = zeros((6,nYears))
Net2 = zeros((8,nYears))

# Total net annual energy yield [EJ/yr]
Net = zeros(nYears)

# Total industrial output [EJ/yr]
Out = zeros(nYears)

# Energy demand for each energy subsector [EJ/yr]
D1 = zeros((6,nYears))
D2 = zeros((8,nYears))

# Total energy demand [EJ/yr]
D = zeros(nYears)
 
# Exploited resource for each energy subsector [/]
# We make the assumption that at its incept date, each energy sub-sector begins 
# with an exploited resource of 0.05.  This value is arbitrary and can of
# course be changed.
rho1 = 0.05 * ones((6,nYears))
rho2 = 0.05 * ones((8,nYears))

# EROI of each energy source [/]
EROI1 = zeros((6,nYears))
EROI2 = zeros((8,nYears))

# Accumulation of capital stock for the three sectors [EJ/yr]
Acc1 = zeros((6,nYears))
Acc2 = zeros((8,nYears))
Acc3 = zeros(nYears)

# Favourability of the energy sources [/]
gamma1 = zeros((6,nYears))
gamma2 = zeros((8,nYears))

# Pre-allocation of vectors and matrices ----- END -----


# Obtained with calibration ----- START -----

URR = array([31500, 14000, 2500, 9050, 1500, 2500]) # [EJ]
TP = array([65, 35, 100, 175, 750, 1, 5, 15])  # [EJ/yr]

Peak_EROI1 = array([71, 400, 60, 350, 10, 15]) # [/]
Peak_EROI2 = array([24, 60, 10, 20, 10, 5, 15, 10]) # [/]

Incept_date1 = array([1800, 1915, 1950, 1925, 1990, 1970]) # [/]
Incept_date2 = array([1800, 1904, 1980, 2005, 2010, 2000, 2020, 2040]) # [/]

XI1 = array([0.904, 0.993, 0.98, 0.996, 0.8, 1.0]) # [/]  The value for unconventional oil is normally 1.0
XI2 = array([0.405, 0.98, 0.979, 0.8, 0.8, 0.75, 0.8, 0.8]) # [/]  Normally hydro=0.984, tidal=1.0

xi1 = array([7.877, 3.570, 0.997, 1.409, 4.0, 8.271]) # [/]
xi2 = array([25.0, 0.066, 0.0, 25.0, 25.0, 25.0, 25.0, 25.0]) # [/]

phi1 = array([6.570, 8.284, 0.001, 6.653, 1.75, 20.0]) # [/]
phi2 = array([0.575, 0.0, 20.0, 18.0, 18.0, 1.561, 18.0, 18.0]) # [/]

eps   = 1.67   # Energy requirement ratio [/]
kappa = 0.0933 # Industrial capital effectiveness [1/yr]

# Obtained with calibration ----- END -----


# The following values are assumed ----- START -----

PHI1 = ones(6) # [/]
PHI2 = ones(8) # [/]

# Capital factor of the energy subsectors [/]
chi1 = 0.1 * ones(6) # [/]
chi2 = 0.9 * ones(8) # [/]

# Lifetime of capital stock for the three sectors [yr]
L1 = 20 * ones(6)
L2 = 20 * ones(8)
L3 = 20

# The following values are assumed ----- END -----


# We store the technological learning component for the EROI of renewable 
# energy sources, so that their technology cannot regress when the exploited 
# resource diminishes
G2 = zeros(8)

Depleted = zeros(6)
# Depleted[i]==1 iff the resources of fossil fuel i have been depleted

for iYear in arange(nYears):
    i1 = logical_and(Incept_date1<=t, logical_not(Depleted.tolist()))      
    i2 = Incept_date2 <= t
    
    # We only consider biomass and tidal energy
    i1[:] = 0
    i2[1:5] = 0
    i2[6:] = 0
        
    # Technological component of the EROI [/]
    G1 = 1 - XI1[i1] * exp(- xi1[i1] * rho1[i1, iYear-1])
    G2[i2] = maximum(G2[i2], 1 - XI2[i2] * exp(- xi2[i2] * rho2[i2, iYear-1]))
        
    # Resource quality component of the EROI [/]
    H1 = PHI1[i1] * exp(- phi1[i1] * rho1[i1, iYear-1])
    H2 = PHI2[i2] * exp(- phi2[i2] * rho2[i2, iYear-1])
    
    EROI1[i1, iYear] = Peak_EROI1[i1] * G1 * H1
    EROI2[i2, iYear] = Peak_EROI2[i2] * G2[i2] * H2
    
    # If the EROI of an energy resource falls below 1, we consider that this resource is depleted
    Depleted[logical_and(Incept_date1<=t, EROI1[:, iYear]<=1)] = 1
    Depleted[logical_and(Incept_date1<=t, rho1[:, iYear-1]>0.99999)] = 1
    
    # We ignore depleted resources
    i1 = logical_and(Incept_date1<=t, logical_not(Depleted.tolist())) 
    
    P1[i1, iYear] = Cap1[i1, iYear-1] * EROI1[i1, iYear] / chi1[i1] / L1[i1]
    P2[i2, iYear] = Cap2[i2, iYear-1] * EROI2[i2, iYear] / chi2[i2] / L2[i2]

    rho1[i1, iYear] = rho1[i1, iYear-1] + P1[i1, iYear] / URR[i1]
    rho2[i2, iYear] = P2[i2, iYear] / TP[i2]
    if(iYear > 1):
        print iYear, " ", P2[0, iYear], " ", EROI2[0, iYear], " ", rho2[0, iYear], " ", TP[0], " ",Cap2[0, iYear-1]
    #print iYear, " ", P2[5, iYear], " ", EROI2[5, iYear], " ", rho2[5, iYear], " ", TP[5]
   
    # We handel the cases when the exploited resources go over 1
    iExceedLim1 = rho1[:, iYear] > 1
    iExceedLim2 = rho2[:, iYear] > 1
    
    P1[iExceedLim1, iYear] -= (rho1[iExceedLim1, iYear] - 1) * URR[iExceedLim1]
    P2[iExceedLim2, iYear] = TP[iExceedLim2]
    rho1[iExceedLim1, iYear] = 1
    rho2[iExceedLim2, iYear] = 1
    
    
    P[iYear] = sum(P1[:,iYear]) + sum(P2[:,iYear])
    
    D[iYear] = eps * kappa * Cap3[iYear-1]
    
    # Proportion of energy demand supplied by each energy subsector [/]
    gamma1[i1, iYear] = rho1[i1, iYear] * EROI1[i1, iYear] / (sum(rho1[i1, iYear] * EROI1[i1, iYear]) + sum(rho2[i2, iYear] * EROI2[i2, iYear]))
    gamma2[i2, iYear] = rho2[i2, iYear] * EROI2[i2, iYear] / (sum(rho1[i1, iYear] * EROI1[i1, iYear]) + sum(rho2[i2, iYear] * EROI2[i2, iYear]))
    
    D1[i1, iYear] = gamma1[i1, iYear] * D[iYear]
    D2[i2, iYear] = gamma2[i2, iYear] * D[iYear]
    
    # Required capital stock for the two energy sectors [EJ/yr]
    Req1 = chi1[i1] * D1[i1, iYear] * L1[i1] / EROI1[i1, iYear]
    Req2 = chi2[i2] * D2[i2, iYear] * L2[i2] / EROI2[i2, iYear]
    
    # Fuel subsidy for the two energy sectors [EJ/yr] 
    S11 = (1 - chi1[i1]) * P1[i1, iYear] / EROI1[i1, iYear]
    S12 = (1 - chi2[i2]) * P2[i2, iYear] / EROI2[i2, iYear]
    
    # Capital subsidy for the two energy sectors [EJ/yr] 
    S21 = maximum(Req1 - Cap1[i1, iYear-1], 0)
    S22 = maximum(Req2 - Cap2[i2, iYear-1], 0)
    
    Net1[i1, iYear] = P1[i1, iYear] - S11
    Net2[i2, iYear] = P2[i2, iYear] - S12
    
    Net[iYear] = sum(Net1[:,iYear]) + sum(Net2[:,iYear])
    
    # We warn the user if the net energy yield becomes negative
    if(Net[iYear]<-1e-08):
        print("year = %d : " %(iYear+t0))
        print("Net energy yield = %8.2f\n" %(Net[iYear]))
    
    Out[iYear] = Net[iYear] / eps
    
    # We handel the case when the required capital subsidies are greater than
    # the available industrial output
    Delta = sum(S21) + sum(S22) - Out[iYear]
    if(Delta>0):
        S21 -= Delta * gamma1[i1, iYear]
        S22 -= Delta * gamma2[i2, iYear]
    
    Acc1[i1, iYear] = S21
    Acc2[i2, iYear] = S22
    Acc3[iYear] = Out[iYear] - sum(S21) - sum(S22)
    
    
    # Rate of depreciation of capital stock for the three sectors [EJ/yr] 
    Dep1 = Cap1[i1, iYear-1] / L1[i1]
    Dep2 = Cap2[i2, iYear-1] / L2[i2]
    Dep3 = Cap3[iYear-1] / L3
    
    Cap1[i1, iYear] = Cap1[i1, iYear-1] + Acc1[i1, iYear] - Dep1
    Cap2[i2, iYear] = Cap2[i2, iYear-1] + Acc2[i2, iYear] - Dep2
    Cap3[iYear] = Cap3[iYear-1] + Acc3[iYear] - Dep3
    t = t + 1



plt.figure()
plt.plot(arange(nYears)+t0, P, label="total")
plt.plot(arange(nYears)+t0, P2[0], label="biomass")
plt.plot(arange(nYears)+t0, P2[5], label="tidal")
plt.legend(loc='upper left', fancybox=True, shadow=True)
plt.xlabel("Year"); plt.ylabel("Production [EJ/year]")

plt.figure()
plt.plot(arange(nYears)+t0, EROI2[0], label="biomass")
plt.plot(arange(nYears)+t0, EROI2[5], label="tidal")
plt.xlabel("Year"); plt.ylabel("EROI [/]")
plt.legend(loc='upper right', fancybox=True, shadow=True)
# 
# plt.figure()
# plt.plot(arange(nYears)+t0, rho2[0], label="biomass")
# plt.plot(arange(nYears)+t0, rho2[5], label="tidal")
# plt.xlabel("Year"); plt.ylabel("Exploited resource")
# plt.legend(loc='upper left', fancybox=True, shadow=True)
# 
# plt.figure()
# plt.plot(arange(nYears)+t0, Cap3, label="industrial")
# plt.legend(loc='upper left', fancybox=True, shadow=True)
# plt.xlabel("Year"); plt.ylabel("Installed capital [EJ]")
# 
# plt.figure()
# plt.plot(arange(nYears)+t0, gamma2[0], label="biomass")
# plt.plot(arange(nYears)+t0, gamma2[5], label="tidal")
# plt.legend(loc='upper left', fancybox=True, shadow=True)
# plt.xlabel("Year"); plt.ylabel("Favourability [/]")
# 
# plt.figure()
# plt.plot(arange(nYears)+t0, Acc3, label="industrial")
# plt.legend(loc='upper left', fancybox=True, shadow=True)
# plt.xlabel("Year"); plt.ylabel("Accumulation of capital [EJ/year]")


'''
plt.figure()
plt.plot(arange(nYears)+t0, P, label="total")
plt.plot(arange(nYears)+t0, P1[0], label="coal")
plt.plot(arange(nYears)+t0, P1[1], label="conventional oil")
plt.plot(arange(nYears)+t0, P1[2], label="unconventional oil")
plt.plot(arange(nYears)+t0, P1[3], label="conventional gas")
plt.plot(arange(nYears)+t0, P1[4], label="unconventional gas")
plt.plot(arange(nYears)+t0, P1[5], label="nuclear fission")
plt.legend(loc='upper left', fancybox=True, shadow=True)
plt.xlabel("Year"); plt.ylabel("Production [EJ/year]")

plt.figure()
plt.plot(arange(nYears)+t0, EROI1[0], label="coal")
plt.plot(arange(nYears)+t0, EROI1[1], label="conventional oil")
plt.plot(arange(nYears)+t0, EROI1[2], label="unconventional oil")
plt.plot(arange(nYears)+t0, EROI1[3], label="conventional gas")
plt.plot(arange(nYears)+t0, EROI1[4], label="unconventional gas")
plt.plot(arange(nYears)+t0, EROI1[5], label="nuclear fission")
plt.xlabel("Year"); plt.ylabel("EROI [/]")
plt.legend(loc='upper right', fancybox=True, shadow=True)

plt.figure()
plt.plot(arange(nYears)+t0, rho1[0], label="coal")
plt.plot(arange(nYears)+t0, rho1[1], label="conventional oil")
plt.plot(arange(nYears)+t0, rho1[2], label="unconventional oil")
plt.plot(arange(nYears)+t0, rho1[3], label="conventional gas")
plt.plot(arange(nYears)+t0, rho1[4], label="unconventional gas")
plt.plot(arange(nYears)+t0, rho1[5], label="nuclear fission")
plt.xlabel("Year"); plt.ylabel("Exploited resource")
plt.legend(loc='upper left', fancybox=True, shadow=True)

plt.figure()
plt.plot(arange(nYears)+t0, gamma1[0], label="coal")
plt.plot(arange(nYears)+t0, gamma1[1], label="conventional oil")
plt.plot(arange(nYears)+t0, gamma1[2], label="unconventional oil")
plt.plot(arange(nYears)+t0, gamma1[3], label="conventional gas")
plt.plot(arange(nYears)+t0, gamma1[4], label="unconventional gas")
plt.plot(arange(nYears)+t0, gamma1[5], label="nuclear fission")
plt.legend(loc='upper left', fancybox=True, shadow=True)
plt.xlabel("Year"); plt.ylabel("Favourability [/]")

plt.figure()
plt.plot(arange(nYears)+t0, P, label="total")
plt.plot(arange(nYears)+t0, P2[0], label="biomass")
plt.plot(arange(nYears)+t0, P2[1], label="hydro")
plt.plot(arange(nYears)+t0, P2[2], label="geothermal")
plt.plot(arange(nYears)+t0, P2[3], label="wind")
plt.plot(arange(nYears)+t0, P2[4], label="solar")
#plt.plot(arange(nYears)+t0, P2[5], label="tidal")
#plt.plot(arange(nYears)+t0, P2[6], label="wave")
plt.plot(arange(nYears)+t0, P2[7], label="otec")
plt.legend(loc='upper left', fancybox=True, shadow=True)
plt.xlabel("Year"); plt.ylabel("Production [EJ/year]")

plt.figure()
plt.plot(arange(nYears)+t0, EROI2[0], label="biomass")
plt.plot(arange(nYears)+t0, EROI2[1], label="hydro")
plt.plot(arange(nYears)+t0, EROI2[2], label="geothermal")
plt.plot(arange(nYears)+t0, EROI2[3], label="wind")
plt.plot(arange(nYears)+t0, EROI2[4], label="solar")
#plt.plot(arange(nYears)+t0, EROI2[5], label="tidal")
#plt.plot(arange(nYears)+t0, EROI2[6], label="wave")
plt.plot(arange(nYears)+t0, EROI2[7], label="otec")
plt.xlabel("Year"); plt.ylabel("EROI [/]")
plt.legend(loc='upper right', fancybox=True, shadow=True)

plt.figure()
plt.plot(arange(nYears)+t0, rho2[0], label="biomass")
plt.plot(arange(nYears)+t0, rho2[1], label="hydro")
plt.plot(arange(nYears)+t0, rho2[2], label="geothermal")
plt.plot(arange(nYears)+t0, rho2[3], label="wind")
plt.plot(arange(nYears)+t0, rho2[4], label="solar")
#plt.plot(arange(nYears)+t0, rho2[5], label="tidal")
#plt.plot(arange(nYears)+t0, rho2[6], label="wave")
plt.plot(arange(nYears)+t0, rho2[7], label="otec")
plt.xlabel("Year"); plt.ylabel("Exploited resource")
plt.legend(loc='upper left', fancybox=True, shadow=True)

plt.figure()
plt.plot(arange(nYears)+t0, gamma2[0], label="biomass")
plt.plot(arange(nYears)+t0, gamma2[1], label="hydro")
plt.plot(arange(nYears)+t0, gamma2[2], label="geothermal")
plt.plot(arange(nYears)+t0, gamma2[3], label="wind")
plt.plot(arange(nYears)+t0, gamma2[4], label="solar")
#plt.plot(arange(nYears)+t0, gamma2[5], label="tidal")
#plt.plot(arange(nYears)+t0, gamma2[6], label="wave")
plt.plot(arange(nYears)+t0, gamma2[7], label="otec")
plt.legend(loc='upper left', fancybox=True, shadow=True)
plt.xlabel("Year"); plt.ylabel("Favourability [/]")
'''


plt.show()




