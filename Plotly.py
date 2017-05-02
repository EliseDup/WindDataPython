import sys
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import scipy.interpolate
from matplotlib.mlab import griddata
import plotly.plotly as py

from numpy import genfromtxt

name = '../WindPotential/eea'

index = 2

def main():
    #plotData(name, index)
    #plotData(name, 2, "eroi", "EROI")
    eroi = (2,5,8,10,12)
    for i in range(0,5):
        print i
        plotData(name, i+2, "eu_gwi"+str(eroi[i]), "Installed Capacity [GW]")
        plotData(name, i+7, "eu_wi"+str(eroi[i]), "Installed Capacity Density [Wi/m2]")
        plotData(name, i+12, "eu_eroi"+str(eroi[i]), "EROI")
    
    print "Hello"
    
def plotData(csvFile, index, output, xLabel, show=False):
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
    mpl_fig = plt.figure()
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
    cbar.set_ticks(np.linspace(0, math.ceil(max(values)), math.ceil(max(values))*2+1))
    
    plotly_fig = tls.mpl_to_plotly(mpl_fig)

    unique_url = py.plot(plotly_fig)
    
    #if show: 
    #    plt.show()
    #else :
    #    plt.savefig(output + '.png', dpi=250, bbox_inches='tight')
    #    plt.close()
        
    return

if __name__ == "__main__":
    sys.exit(main())
