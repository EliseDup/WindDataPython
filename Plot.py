import sys
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import scipy.interpolate
from matplotlib.mlab import griddata

from numpy import genfromtxt
solarFolder = '../resources/era_interim_data/solar/'

name = '../WindPotential/clearnessIndex'#solarFolder+'net40years' #../WindPotential/res'

index = 2
def main():
    plotData(name,index, output="clearnessIndex", xLabel="Clearness Index")
    print "Hello"
    
def plotData(csvFile, index, output, xLabel="", show=False):
    data = genfromtxt(csvFile, delimiter='\t', dtype=None)
    # data = genfromtxt(csvFile, delimiter='\t', dtype=None)
    #### data preparation 
    lats = data[:, 0] 
    # # lon => x 
    lons = data[:, 1] 
    # # values => z 
    values = data[:, index]
    print "Size :", len(lats)
    print "Lat Min :", min(lats), "Lat Max :", max(lats)
    print "Lon Min :", min(lons), "Lon Max :", max(lons)
    #### later in the defined map 
    map = Basemap(projection='cyl', llcrnrlat=min(lats), urcrnrlat=max(lats), llcrnrlon=min(lons), urcrnrlon=max(lons), resolution='l')
    
    map.drawcoastlines()
    map.drawstates()
    map.drawcountries()
    # map.drawlsmask(land_color='coral',ocean_color='blue') 
    print "Maximum ", max(values), " - Minimum ", min(values)
    cs = map.contourf(lons, lats, values, np.linspace(0.01, max(values)+0.01
                                                      , 250, endpoint=True), tri=True ) #, latlon=True)
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label(xLabel)
    cbar.set_ticks(np.linspace(0, 1, 11)) #0.0,math.ceil(max(values)), math.ceil(max(values))*2+1))
    
    if show: 
        plt.show()
    else :
        plt.savefig(output + '.png', dpi=250, bbox_inches='tight')
        plt.close()
        
    return

if __name__ == "__main__":
    sys.exit(main())
