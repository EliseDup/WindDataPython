'''

Sector 1 is the non-renewable part of the energy sector.
It is composed of:
    0: coal
    1: conventional oil
    2: unconventional oil
    3: conventional gas
    4: unconventional gas
    5: nuclear fission

Sector 2 is the renewable part of the energy sector.
It is composed of:
    0: biomass
    1: hydro
    2: geothermal
    3: wind
    4: solar
    5: tidal
    6: wave
    7: otec
    
Sector 3 is the industrial sector.
It is not broken into subsectors.

'''

from numpy import*
import matplotlib.pyplot as plt
import sys

non_renewable_sectors = array(["Coal", "Conventional Oil", "Unconventional Oil", "Conventional Gas", "Unconventional Gas", "Nuclear Fission" ])
renewable_sectors = array(["Biomass", "Hydro", "Geothermal", "Wind", "Solar", "Tidal", "Wave", "Otec"])
n1 = non_renewable_sectors.size
n2 = renewable_sectors.size

# Obtained with calibration ----- START -----
URR = array([31500, 14000, 2500, 9050, 1500, 2500])  # [EJ]
TP = array([65, 35, 100, 175, 750, 1, 5, 15])  # [EJ/yr]

Peak_EROI1 = array([71, 400, 60, 350, 10, 15])  # [/]
Peak_EROI2 = array([24, 60, 10, 20, 10, 5, 15, 10])  # [/]

Incept_date1 = array([1800, 1915, 1950, 1925, 1990, 1970])  # [/]
Incept_date2 = array([1800, 1904, 1980, 2005, 2010, 2000, 2020, 2040])  # [/]

XI1 = array([0.904, 0.993, 1.0, 0.996, 0.8, 1.0])  # [/]
XI2 = array([0.405, 0.984, 0.979, 0.8, 0.8, 1.0, 0.8, 0.8])  # [/]

xi1 = array([7.877, 3.570, 0.997, 1.409, 4.0, 8.271])  # [/]
xi2 = array([25.0, 0.066, 0.0, 25.0, 25.0, 25.0, 25.0, 25.0])  # [/]

phi1 = array([6.570, 8.284, 0.001, 6.653, 1.75, 20.0])  # [/]
phi2 = array([0.575, 0.0, 20.0, 18.0, 18.0, 1.561, 18.0, 18.0])  # [/]

eps = 1.67  # Energy requirement ratio [/]
kappa = 0.0933  # Industrial capital effectiveness [1/yr]


# Obtained with calibration ----- END -----

# The following values are assumed ----- END -----
PHI1 = ones(6)  # [/]
PHI2 = ones(8)  # [/]

# We store the technological learning component for the EROI of renewable 
# energy sources, so that their technology cannot regress when the exploited 
# resource diminishes
G2 = zeros(8)

def main():
    #plotEROI(False)
    simulate(151)
    print("End Simulation")

def EROI_renewable(i, rho):
    # Technological component of the EROI [/]
    G = maximum(G2[i], 1 - XI2[i] * exp(-xi2[i] * rho))
   # Resource quality component of the EROI [/]
    H = PHI2[i] * exp(-phi2[i] * rho)
    return Peak_EROI2[i] * G * H
   
def EROI_non_renewable(i, rho):
    # Technological component of the EROI [/]
    G = 1 - XI1[i] * exp(-xi1[i] * rho)
    # Resource quality component of the EROI [/]
    H = PHI1[i] * exp(-phi1[i] * rho)
    return Peak_EROI1[i] * G * H 

def plotEROI(save):
    rho_test = linspace(0, 1, 100)
    plt.figure()
    for i in range(0, renewable_sectors.size):
        plt.subplot(2, 4, i + 1)
        plt.plot(rho_test * TP[i], EROI_renewable(i, rho_test))
        plt.title(renewable_sectors[i])
        plt.draw()
        #plt.show()
    if(save):
        plt.savefig('EROI_renewables.png', dpi=250, bbox_inches='tight')
        plt.close()
        
    plt.figure()
    for i in range(0, non_renewable_sectors.size):
        plt.subplot(2, 3, i + 1)
        plt.plot(rho_test * URR[i], EROI_non_renewable(i, rho_test))
        plt.title(non_renewable_sectors[i])
        plt.draw()
    if(save):
        plt.savefig('EROI_non_renewables.png', dpi=250, bbox_inches='tight')
        plt.close()
    else: 
        plt.show()
    return

# nYears = Total number of time steps of the simulation
# One time step = one year, we go from the year 1800 to the year 2200
def simulate(nYears):
    # Pre-allocation of vectors and matrices ----- START -----

    # Capital stock of the three sectors [EJ].
    # We make the assumption that at its incept date, each energy sub-sector begins 
    # with an installed capital of 1 [EJ].  This value is arbitrary and can of
    # course be changed.  We assume that the main economy begins with 0 capital
    # stock in 1800.
    Cap1 = ones((6, nYears))
    Cap2 = ones((8, nYears))
    Cap3 = zeros(nYears)

    # Annual energy production of the two energy sectors [EJ/yr]
    P1 = zeros((6, nYears))
    P2 = zeros((8, nYears))

    # Total annual energy production [EJ/yr]
    P = zeros(nYears)

    # Net annual energy yield of the two energy sectors [EJ/yr]
    Net1 = zeros((6, nYears))
    Net2 = zeros((8, nYears))

    # Total net annual energy yield [EJ/yr]
    Net = zeros(nYears)

    # Total industrial output [EJ/yr]
    Out = zeros(nYears)

    # Energy demand for each energy subsector [EJ/yr]
    D1 = zeros((6, nYears))
    D2 = zeros((8, nYears))

    # Total energy demand [EJ/yr]
    D = zeros(nYears)
 
    # Exploited resource for each energy subsector [/]
    rho1 = zeros((6, nYears))
    rho2 = zeros((8, nYears))

    # EROI of each energy source [/]
    EROI1 = zeros((6, nYears))
    EROI2 = zeros((8, nYears))

    # Capital factor of the energy subsectors [/]
    chi1 = 0.1 * ones(6)  # [/]
    chi2 = 0.9 * ones(8)  # [/]

    # Lifetime of capital stock for the three sectors [yr]
    L1 = 20 * ones(6)
    L2 = 20 * ones(8)
    L3 = 20 
    t0 = 1800
    
    for iYear in arange(nYears):    
    
        t = t0 + iYear
        i1 = Incept_date1 <= t
        i2 = Incept_date2 <= t
    
        EROI1[i1, iYear] = EROI_non_renewable(i1, rho1[i1, iYear - 1])
        EROI2[i2, iYear] = EROI_renewable(i2, rho2[i2, iYear - 1])
        
        # Elise -> Add a minimum to make sure you don't produce more than the potential ?
        for i in arange(n1):
            if(i1[i]):
                P1[i, iYear] = min(URR[i],Cap1[i, iYear - 1] * EROI1[i, iYear] / chi1[i] / L1[i])
                URR[i] = URR[i] - P1[i, iYear]
        for i in arange(n2):
            if(i2[i]):
                P2[i, iYear] = min(TP[i], Cap2[i, iYear - 1] * EROI2[i, iYear] / chi2[i] / L2[i])
        
        P[iYear] = sum(P1[:, iYear]) + sum(P2[:, iYear])
    
        rho1[i1, iYear] = rho1[i1, iYear - 1] + P1[i1, iYear] / URR[i1]
        rho2[i2, iYear] = P2[i2, iYear] / TP[i2]
     
        D[iYear] = eps * kappa * Cap3[iYear - 1]
    
        # Proportion of energy demand supplied by each energy subsector [/]
        gamma1 = rho1[i1, iYear] * EROI1[i1, iYear] / sum(rho1[i1, iYear] * EROI1[i1, iYear])
        gamma2 = rho2[i2, iYear] * EROI2[i2, iYear] / sum(rho2[i2, iYear] * EROI2[i2, iYear])
    
        D1[i1, iYear] = gamma1 * D[iYear]
        D2[i2, iYear] = gamma2 * D[iYear]
    
        # Required capital stock for the two energy sectors [EJ/yr]
        Req1 = chi1[i1] * D1[i1, iYear] * L1[i1] / EROI1[i1, iYear]
        Req2 = chi2[i2] * D2[i2, iYear] * L2[i2] / EROI2[i2, iYear]
    
        # Fuel subsidy for the two energy sectors [EJ/yr] 
        S11 = (1 - chi1[i1]) * P1[i1, iYear] / EROI1[i1, iYear]
        S12 = (1 - chi2[i2]) * P2[i2, iYear] / EROI2[i2, iYear]
    
        # Capital subsidy for the two energy sectors [EJ/yr] 
        S21 = Req1 - Cap1[i1, iYear - 1]
        S22 = Req2 - Cap2[i2, iYear - 1]
    
        Net1[i1, iYear] = P1[i1, iYear] - S11
        Net2[i2, iYear] = P2[i2, iYear] - S12
    
        Net[iYear] = sum(Net1[:, iYear]) + sum(Net2[:, iYear])
    
        Out[iYear] = Net[iYear] / eps
    
        # Rate of accumulation of capital stock for the three sectors [EJ/yr]
        Acc1 = maximum(S21, 0)
        Acc2 = maximum(S22, 0)
        Acc3 = Out[iYear] - sum(S21) - sum(S22)
    
        # Rate of depreciation of capital stock for the three sectors [EJ/yr] 
        Dep1 = Cap1[i1, iYear - 1] / L1[i1]
        Dep2 = Cap2[i2, iYear - 1] / L2[i2]
        Dep3 = Cap3[iYear - 1] / L3
    
        Cap1[i1, iYear] = Cap1[i1, iYear - 1] + Acc1 - Dep1
        Cap2[i2, iYear] = Cap2[i2, iYear - 1] + Acc2 - Dep2
        Cap3[iYear] = Cap3[iYear - 1] + Acc3 - Dep3
        
        
    # Plot the annual production of energy
    plt.figure()
    plt.plot(arange(nYears)+t0,P);
    for i in arange(n1):
        plt.plot(arange(nYears)+t0,P1[i,:], label=non_renewable_sectors[i]);
    for i in arange(n2):    
        plt.plot(arange(nYears)+t0,P2[i,:], label=renewable_sectors[i])
    plt.legend(loc='upper right', bbox_to_anchor=(0.5, 1.05),
          ncol=3, fancybox=True, shadow=True)
    plt.xlabel("Year"); plt.ylabel("Production [EJ/year]")
    plt.show()
    return

if __name__ == "__main__":
    sys.exit(main())
 