import sys
import pygrib
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
from ncepgrib2 import dump

VERBOSE = 1  # verbose error reporting

def main():
    speedToFile("../ressources/wind/europe5years",50,50)
    
def gribToFile(name, data_name):
    grib_data = pygrib.open(name +'.grib')
    data, lats, lons = grib_data.select(name=data_name)
    res = open(name, 'w')
    for i in range(0, len(lats)):
        for j in range(0, len(lats[0])):
            res.write(str(lats[i][j]));res.write("\t");
            res.write(str(lons[i][j]));res.write("\t");
            res.write(str(data[i][j]));res.write("\n");
            
    res.close()
    return

def speedToFile(name, i, j):
    grib_data = pygrib.open(name +'.grib')
    uComponents = grib_data.select(name='10 metre U wind component') # To the EAST
    vComponents = grib_data.select(name='10 metre V wind component') # To the NORTH
    res = open('test', 'w')
    print len(uComponents),"-",len(vComponents)
    for m in range(0, min(len(vComponents),len(uComponents))):
        dataU, latsU, lonsU = uComponents[m].data(); dataV, latsV, lonsV = vComponents[m].data()
        if(m%100 == 0): print uComponents[m] , " ", vComponents[m], " ", latsV[i][j], " ", lonsV[i][j]
        
        res.write(str(dataU[i][j]));res.write("\t");
        res.write(str(dataV[i][j]));res.write("\t");
        res.write(str(math.sqrt(dataU[i][j]*dataU[i][j]+dataV[i][j]*dataV[i][j])));res.write("\n");
            
    res.close()
    return

if __name__ == "__main__":
    sys.exit(main())