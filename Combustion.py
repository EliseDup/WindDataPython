# Masses molaires en [kg/mol]
Mm_O2 = 0.032;
Mm_N2 = 0.028;
Mm_CH4 = 0.016;
Mm_CO2 = 0.044;
Mm_H2O = 0.018;

# Cp moyen entre 0 et 1000C [kJ/kgK]
cp_O2  = 1.035 # 0.03312 / Mm_O2 
cp_N2  = 1.12 # 0.03131 / Mm_N2 
cp_CH4 = 3.78 # 0.06050 / Mm_CH4
cp_CO2 = 1.12 # 0.04939 / Mm_CO2
cp_H2O = 2.13 # 0.03847 / Mm_H2O

print cp_O2, cp_N2, cp_CH4, cp_CO2, cp_H2O

PCI_CH4 = 50e3; # [kJ/kg]
T0 = 0; #[C]
T1 = 15; #[C]

eN2 = cp_N2*Mm_N2
eO2 = cp_O2*Mm_O2
eCH4 = cp_CH4*Mm_CH4
eCO2 = cp_CO2*Mm_CO2
eH2O = cp_H2O*Mm_H2O

# CH4 + 2 lam O2 + 2 w lam N2 => CO2 + 2H20 + 2 w lam N2 + 2 (lam-1) 02
#CAS 1
w = 0 # Pas de N2
lam = 1 #Pas d'exces d'air
q = PCI_CH4*Mm_CH4
def T2(w, lam):
    return T0 + (Mm_CH4*(PCI_CH4 + cp_CH4*(T1-T0))+2*lam*Mm_O2*cp_O2*(T1-T0)+2*lam*w*Mm_N2*cp_N2*(T1-T0))/(Mm_CO2*cp_CO2 + 2*Mm_H2O *cp_H2O + 2*(lam-1)*cp_O2*Mm_O2 + 2*lam*w*Mm_N2*cp_N2)
def T_f(w,lam):
    return (q + T1*(eCH4 + 2*lam*eO2 + 2*lam*w*eN2))/(eCO2+2*eH2O+2*lam*w*eN2+2*(lam-1)*eO2)

def lam(Tf, w):
    return (q - Tf*(eCO2 + 2*eH2O - 2*eO2) + T1*(eCH4)) / (Tf*(2*w*eN2+2*eO2)-T1*(2*eO2+2*w*eN2))

print "Only O2 ", T2(0,1), T_f(0,1)
print "With N2, no excess", T2(0.79/0.21,1),T_f(0.79/0.21,1)
print "Excess", T2(0.79/0.21,2),T_f(0.79/0.21,2)
print "Tf = 1400", lam(1400,0.79/0.21)
print "lambda = 17.4", T_f(0.79/0.21,1.71385973206)
# Lambda such as Tf = 1400 ?

#figure;
#xx = [1:0.01:8];
#yy = (T0 + (Mm_CH4*(PCI_CH4 + cp_CH4*(T1-T0))+2*xx*Mm_O2*cp_O2*(T1-T0)+2*xx*w*Mm_N2*cp_N2*(T1-T0))./(Mm_CO2*cp_CO2 + 2* Mm_H2O * cp_H2O + 2*(xx-1)*cp_O2*Mm_O2 + 2*xx*w*Mm_N2*cp_N2));

#plot(xx,yy)
