# -*- coding: utf-8 -*-
import sys
import numpy as np
from numpy import genfromtxt, ndarray
import math
from datetime import datetime
from scipy.special import gammainc
from matplotlib.pyplot import plot
import matplotlib.pyplot as plt

# model macro-energie : calibration

#donnees (2017)
A=817.0;
E=9717.0;
U=A+E;
Ey=6316.0;
Ce=E-Ey;
PIB=80250.0;
PIBm1=77796.0; #PIB annee precedente
L=3422.0; # Pop active (10^6 personnes) 
L=2871.0; # Pop employee (10^6 personnes)
gPT=0.001; #.00534; #taux de progrès technique
s=.25;
T=25.0;
al=.05;#VAe/PIB
th=.4;#rK/PIB
m=.05;#.05/.95;#Ke/Ky

#calibration
ga=Ey/E;
qe=A/U;
g=PIB/PIBm1-1;
gk=g-gPT;
delta=1.0-math.pow(1.0-0.9,1.0/T);
v=s/(gk+delta);
K=v*PIB;
Ke=m*K/(1+m);
Ky=K/(1+m);
ve=Ke/U;
p=al*PIB/E;
Y=PIB-p*Ce;
vy=Ky/Y;
qy=Ey/Y;
eps=1/(qe+delta*qy*ve)

print v, qe, delta, vy, qy, ve, eps


