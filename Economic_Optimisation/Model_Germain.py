# -*- coding: utf-8 -*-
import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot
import matplotlib.pyplot as plt

# Variables = EROI et qe (autocosommation d'energie) évoluent au cours du temps (ou ve et EROI, ou ve et qe)

# Données "observables"
qy = 1.4  # 1.4 # énergie finale par unité de PIB[ue/ub]
# vy = 0.1 # intensité capitalistique de l'économie [ub/ub]
m = 0.1  # partage du capital entre énergie et le reste de l'économie

delta = 1 / 20.0  # durée de vie du capital
s = 0.1  # taux d'épargne = Y/I

# Calibration
# vy, capital par unité de PIB produit [ub/ub]
def _vy(m, eroi, qe):
    return _ve(eroi, qe) / m * qy / (1 - qe)

# ve, capital par unité d'énergie finale [ub/ue], défini par le EROI et qe (auto consommation d'énergie), et par qy l'intensité énergétique de l'économie
def _ve(eroi, qe):
    return (1.0 - eroi * qe) / (eroi * delta * qy)
   
# m, partage du capital entre énergie et économie (=Ke/Ky)
def _m(eroi, qe, vy):
    return ve(eroi, qe) / vy * qy / (1 - qe)
# v, intensité capitalistique de toute l'économie 
def _v(eroi, qe):
    return vy * (1 + m(eroi, qe))
def _delta_v(eroi0, qe0, eroi1, qe1):
    return v(eroi1, qe1) - v(eroi0, qe0)
def _g(eroi0, qe0, eroi1, qe1):
    return (s - delta_v(eroi0, qe0, eroi1, qe1)) / (vy + 1 / delta * (1 / eroi1 - qe1) / (1 - qe1)) - delta

# Exogenous data
tf = 100
eroi = np.linspace(20, 5, tf)
qe = np.ones(tf) * 0.01
ve = _ve(eroi, qe)

def main():
    calibrationMarc()
    y = np.zeros(tf); c = np.zeros(tf);
    e = np.zeros(tf); u = np.zeros(tf); a = np.zeros(tf);
    # k = np.zeros(tf);
    ky = np.zeros(tf);
    ke = np.zeros(tf);
    m = np.zeros(tf);
    s = np.zeros(tf);
    vy = np.zeros(tf); 
    # t = 0
    y[0] = 1; s[0] = 0.1; m[0] = 0.2;
    
    for t in range(1, tf):
         # begin with update value in t-1
        e[t - 1] = qy * y[t - 1]; 
        u[t - 1] = (1 + qe[t - 1]) * e[t - 1];
        vy[t - 1] = _vy(m[t - 1], eroi[t - 1], qe[t - 1])
        ky[t - 1] = y[t - 1] * vy[t - 1]
        ke[t - 1] = u[t - 1] * ve[t - 1]
        m[t - 1] = ke[t - 1] / ky[t - 1]
        
        # estimate the new values for ? Y ? S ? M ?
        y[t] = 1; s[t] = 0.1; m[t] = 0.2;

    #plot(e)
    #plt.show()
    print "-- END --"
    
def calibrationMarc():
# donnees (2017)
    U = 10537.0;
    A = 817.0;
    E = 9717.0;
    Ey = 6316.0;
    Ce = 3401.0;
    PIB = 80250.0;
    PIBm1 = 77796.0;  # PIB annee precedente
    gPT = .00534;  # taux de progrès technique
    s = .25;
    T = 15.0;
    al = .08;

    # calibration
    ga = Ey / E;
    qe = A / U;
    g = PIB / PIBm1 - 1;
    gk = g - gPT;
    delta = 1 - math.pow(1 - .9,1 / T);
    v = s / (gk + delta);
    m = al / (1.0 - al);
    rho = 1.0 - al + al * ga;
    p = al * PIB / E;
    r = 1.0 / v;
    qy = al * ga / p / rho;
    vy = v * (1 - al) / rho;
    ve = al * ga * (1 - qe) * v / rho / qy;
    eps = 1 / (qe + delta * qy * ve);
    # contributions a l'eroi
    qe, delta * qy * ve
    # verif : on doit verifier l'egalite suivante
    'v,vy+ve*qy/(1-qe),rho*(vy+ve*qy/(1-qe)/ga)', [v, vy + ve * qy / (1 - qe), rho * (vy + ve * qy / (1 - qe) / ga)]

    # calcul variables macro
    Y = rho * PIB;
    Ky = vy * Y;
    Ke = ve * U;
    K = Ke + Ky;
    VAy = (1 - p * qy) * Y;
    VAe = p * E;
    I = s * PIB;
    Cy = Y - I;

    print "Calibration Marc", p, ve, vy, qy, eps
    # verif : on doit verifier les egalites suivantes
    'qy*Y,Ey', [qy * Y, Ey]
    'm,Ke/Ky', [m, Ke / Ky]
    'K/v,PIB', [K / v, PIB]
    'Y,p*Ey+r*Ky', [Y, p * Ey + r * Ky]
    'p*E,r*Ke', [p * E, r * Ke]
    'PIB,VAy+VAe,r*K', [PIB, VAy + VAe, r * K]
    'PIB,Y+p*Ce', [PIB, Y + p * Ce]
    'U,A+E,A+Ey+Ce', [U, A + E, A + Ey + Ce]
    'Cy+p*Ce,(1-s)*PIB', [Cy + p * Ce, (1 - s) * PIB]
    return

if __name__ == "__main__":
    sys.exit(main())
