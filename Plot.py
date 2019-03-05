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
file = ''
name = '../WindPotentialScala/'+file
# sf_wind, wind100m, cf_wind_100m, wi_eroi5, wi_eroi12
index = 7

def main():
    #plotData('../WindPotentialScala/eu28',7,'eu28meanSlope',show=True)
    #plotData('wind_solar_grid/technology_12.0', 2, "technology_12",show=False)
    #plotDataPerTech(2, 'eroi', "EROI",('1'))
    eroi = ('4','6','9','12')
    plotDataPerTech(4, 'power', "Power Density [We/m2]",eroi)
    #plotDataPerTech(7, 'installed_power', "Installed Capacity Density [We/m2]",eroi)

def plotDataPerTech(index, label_short, label, eroi_list):
    #techs = ('Wind-onshore','Wind-offshore','PT-oil','PT-salt-TES','ST-salt-TES','Poly-Si-PV','Mono-Si-PV')
    techs = 'Wind' # Wind-onshore','Wind-offshore','ST-salt-TES','Mono-Si-PV')
    #for t in techs:
    for e in eroi_list:
        plotData('wind_solar_grid/'+techs+"_"+e,index, output='wind_solar_images/'+techs+"_"+e+'_'+label_short, xLabel=label, show=False)
    
def plotData(csvFile, index, output, xLabel="", show=True):
    
    data = genfromtxt(csvFile, delimiter='\t', dtype=None)
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
    if(max(values)==0):
        return
    cs = map.contourf(lons, lats, values,
                      np.linspace(min(values)+0.01, max(values)+0.01, 250, endpoint=True), tri=True ) #, latlon=True)
                      #np.linspace(0.01, 3, 250, endpoint=True), tri=True ) #, latlon=True
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label(xLabel)
    cbar.set_ticks(np.arange(0,1.1,0.1))
    
    if show: 
        plt.show()
    else :
        # High resolution
        #plt.savefig(output + '.pdf', dpi=250, bbox_inches='tight')
        plt.savefig(output + '.png', dpi=250, bbox_inches='tight')
        plt.close()
        
    return

if __name__ == "__main__":
    sys.exit(main())
