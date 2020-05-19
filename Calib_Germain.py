# modelta macro-energie : calibration
import math
#donnees (2017)
U=10537.0;
A=820.0;
E=9717.0;
Ey=6316.0;
Ce=3401.0;
PIB=80250.0;
PIBm1=77796.0; #PIB annee precedente
L=3422.0; # Pop active (10^6 personnes) 
L=2871.0; # Pop employee (10^6 personnes)
gPT=0.0; #.00534; #taux de progres technique
s=.25;
T=15;
al=.06;#VAe/PIB
th=.4;#rK/PIB
m=.05;#.05/.95;#Ke/Ky

#calibration
ga=Ey/E;
qe=A/U;
g=PIB/PIBm1-1;
gk=g-gPT;
delta=1-math.pow(0.1,1.0/T);
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
#contributions a l'eroi
qe,delta*qy*ve

r=th*(gk+delta)/s;
w=(1-th)*PIB/L;
ly=(1-p*qy-r*vy)/w;
le=((1-qe)*p-r*ve)/w;
Ly=ly*Y;
Le=le*U;
rho=Y/PIB;

#autre variables 
VAy=(1-p*qy)*Y;
VAe=p*E;
I=s*PIB;
Cy=Y-I;

n=p*Ce/Cy; #coefficient de partage du budget de C
sy=s*(1+n)/(1+s*n); #I/Y
e=E/PIB;
c=Ce/Cy;
k=K/U;
#kmax=(1-qe)*vy/qy+ve; # Ce>0 => k<kmax
#kmax=(1-qe)*(1+s*n)*vy/(1+n)/qy+ve; # p<1/qy => k<kmax
#kmin=ve;
#intervalle de k <- w,r + (cas ou ve/vy<p<le/ly)
n2=(1-s)*n/(1+s*n);
n1=(1+n)/(1+s*n); #var aux
ki=ve+vy/(n2*vy/ve+n1*qy/(1-qe)); #borne inf
ks=ve+vy/(n2*ly/le+n1*qy/(1-qe)); #borne sup
k0=ve/(1-vy*(delta+gPT)*(1+s*n)/(1+n)/s);
gi=sy*(1-ve/ki)/vy-delta;
gs=sy*(1-ve/ks)/vy-delta;
'gi,s/ve/((1-s)*n*vy/(1+n)/ve+qy/(1-qe)+(1+s*n)*vy/(1+n)/ve)-delta-gPT',gi,s/ve/((1-s)*n*vy/(1+n)/ve+qy/(1-qe)+(1+s*n)*vy/(1+n)/ve)-delta-gPT
'gs,s/ve/((1-s)*n*ly/(1+n)/le+qy/(1-qe)+(1+s*n)*vy/(1+n)/ve)-delta-gPT',gs,s/ve/((1-s)*n*ly/(1+n)/le+qy/(1-qe)+(1+s*n)*vy/(1+n)/ve)-delta-gPT

print gi, g ,gs
#verif : on doit verifier les egalites suivantes
'qy*Y,Ey',[qy*Y,Ey]
'm,Ke/Ky',[m,Ke/Ky]
'K/v,PIB',[K/v,PIB]
'K,th*PIB/r',[K,th*PIB/r]
'Y,p*Ey+r*Ky+w*Ly',[Y,p*Ey+r*Ky+w*Ly]
'p*E,r*Ke+w*Le',[p*E,r*Ke+w*Le]
'PIB,VAy+VAe,r*K+w*L',[PIB,VAy+VAe,r*K+w*L]
'PIB,Y+p*Ce',[PIB,Y+p*Ce]
'U,A+E,A+Ey+Ce',[U,A+E,A+Ey+Ce]
'Cy+p*Ce,(1-s)*PIB',[Cy+p*Ce,(1-s)*PIB]
'K,Ke+Ky',[K,Ke+Ky]
'L,Le+Ly',[L,Le+Ly]
'p,(r*ve+w*le)/(1-qe)',[p,(r*ve+w*le)/(1-qe)]
'1-p*qy,r*vy+w*ly',[1-p*qy,r*vy+w*ly]
'K/E,ga*vy/qy+ve/(1-qe)',[K/E,ga*vy/qy+ve/(1-qe)]
'L/E,ga*ly/qy+le/(1-qe)',[L/E,ga*ly/qy+le/(1-qe)]
'g,gk+gPT,(1+n)*s*(1-ve/k)/vy/(1+s*n)-delta-gPT',g,gk+gPT,(1+n)*s*(1-ve/k)/vy/(1+s*n)-delta-gPT

#statique

#1) U,c,k donnes 
# ER=10#
Ur=1097.8861130882
Ar=5.0638975531
Unr=8644.4360709112
epsp=10.8985177298067;
# ER=20#
#Ur=2057.39822586319
#Ar=9.86316885054747
#Unr=7684.9239581362
#epsp=9.987632329
# ER=50#
Ur=4941.656940314
Ar=25.0128328287
Unr=4800.6652436853
epsp=7.755835605

qer=Ar/Ur
pr=Ur/(Ur+Unr)
qep=qer*pr+qe*(1-qe)
deltap=delta
vep=(1/epsp-qep)/deltap/qy;
np=n;
Up=1*U;
kp=2*k; 
syp=s*(1+np)/(1+s*np);
n2p=(1-s)*np/(1+s*np); #var aux
n1p=(1+np)/(1+s*np); #var aux
kip=vep+vy/(n2p*vy/vep+n1p*qy/(1-qep)); 
ksp=vep+vy/(n2p*ly/le+n1p*qy/(1-qep)); #borne sup
k0p=vep/(1-vy*(deltap+gPT)*(1+s*np)/(1+np)/s);
Kp=kp*Up;
Ep=(1-qep)*Up;
Kep=vep*Up;
Kyp=Kp-Kep;
mp=Kep/Kyp;
Yp=Kyp/vy;
Eyp=qy*Yp;
Cep=Ep-Eyp;
gap=Eyp/Ep;
Cyp=(1-s)*Yp/(1+s*n);
Ip=Yp-Cyp;
PIBp=Ip/s;
pp= np*Cyp/Cep;
'pp,((1-s)*PIBp-Cyp)/Cep',pp,((1-s)*PIBp-Cyp)/Cep
'PIBp,Cyp+pp*Cep+Ip',PIBp,Cyp+pp*Cep+Ip
#if pp>1/qy
#  'probleme'
#  end
ep=Ep/PIBp;
cp=Cep/Cyp;
gp=s*PIBp/Kp-delta;
gip=syp*(1-vep/kip)/vy-deltap;
gsp=syp*(1-vep/ksp)/vy-deltap;
print gip, gp, gsp
#matp=[vy,ly;vep,le]
#inv(matp)*[1-pp*qy;pp*(1-qep)]











