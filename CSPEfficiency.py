import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import genfromtxt

def main():
    dni_efficiency("../trough",3)
  
def dni_efficiency(file, degree):
    data = genfromtxt(file, delimiter='\t', dtype=None)
    
    x = data[:, 0]; y = data[:, 1] 
    z = np.poly1d(np.polyfit(x, y, degree))
    print z
    xp = np.linspace(500, 3000, 100)
   
    plt.figure; plt.plot(x,y,'.',xp,z(xp))
    plt.xlabel("DNI")
    plt.ylabel("Efficiency [%]")
    plt.title(file)
    
    plt.show()
    
    return z
    
def troughNoStorage(degree):
    x = np.array([0.5, 0.75, 1.0, 1.25, 1.3, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 0.5, 0.75, 1.0, 1.25, 1.3, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3, 0.5, 0.75, 1.0, 1.25, 1.3, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3]) 
    y = np.array([8.89,11.88,13.04,13.44,13.44,13.00,12.17,11.25,10.45,9.74,9.11,8.51,10.73,13.35,14.36,14.76,14.62,13.90,12.77,11.65 ,10.71,9.88,9.18 ,8.54,11.06,12.85,13.61,13.77,13.52,12.62,11.49,10.42,9.53,8.74,8.09,7.50])
    xp = np.linspace(0, 3, 100)
    z = np.poly1d(np.polyfit(x, y, degree))
    plt.figure; plt.plot(x,y,'.',xp,z(xp))
    plt.xlabel("Solar Multiple")
    plt.ylabel("Efficiency [%]")
    plt.title("Parabolic Trough No Storage")
    plt.show()
    
    return z

def trough12HoursStorage(degree):
    x = np.array([1,1.5,2,2.25,2.5,2.75,3,3.25,3.5,4,4.5,5,5.5,6,1,1.5,2,2.25,2.5,2.75,3,3.25,3.5,4,4.5,5,5.5,6,1,1.5,2,2.25,2.5,2.75,3,3.25,3.5,4,4.5,5,5.5,6]) 
    y = np.array([13.6738966164829,    14.187516819295,    14.6871832728589,    14.8352319519263,    14.7666836309427,    14.2837019403469,    13.7254867918645,    13.1756947422504,    12.5990100894301,    11.5542739186743,    10.6244989911952,    9.76724234365377,    9.00941347400859,    8.32275935700974,    13.0278914038241,    13.6708568985666,    14.1222468262367,    14.2942642005532,    14.455143282423,    14.5774844217982,    14.3309882000552,    13.9908434771544,    13.5790083108125,    12.6702348500374,    11.8096011241008,    10.9584950164815,    10.216532855903,    9.50302033008429,    14.4628398198745,    15.0858219533938,    15.6627371688348,    15.856220514364,    16.0088182696013,    15.9564400965044,    15.6754217674983,    15.2634606050314,    14.7132499857261,    13.5619262694737,    12.4350487255124,    11.3993005612062,    10.4916338698931,    9.67573876127865,])
    xp = np.linspace(0, 6, 100)
    z = np.poly1d(np.polyfit(x, y, degree))
    plt.figure; plt.plot(x,y,'.',xp,z(xp))
    plt.xlabel("Solar Multiple")
    plt.ylabel("Efficiency [%]")
    plt.title("Parabolic Trough 12 Hours Storage")
    plt.show()
    
    return z


if __name__ == "__main__":
    sys.exit(main())
 