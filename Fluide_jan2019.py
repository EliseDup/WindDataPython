import numpy as np
import matplotlib.pyplot as plt
import sys

# Graphes question fluide janvier 2019
muA = 1.0
muB = muA / 5.0
p0 = 10.0
pL = 1.0
L = 1.0
b = 1.0

def p(x1):
    return (pL-p0)/L * x1 +p0

def v1(x2, mu):
    return ((pL-p0)/L) * (x2*x2/(2*mu) - b/(2*mu)*(muA-muB)/(muA+muB)*x2 - b*b/(muA+muB))

x1 = np.linspace(0, L, 100)
x2A = np.linspace(-b,0,100)
x2B = np.linspace(0,b,100)

plt.plot(v1(x2A,muA),x2A,v1(x2B,muB),x2B)
plt.ylabel("x2/b")
plt.xlabel("v1")
plt.show()