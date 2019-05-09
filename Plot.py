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

def main():
     plotData("results_simple_model", 2, "out")
    #plotData("../WindPotentialScala/EROI", 2,'EROI_PV',xLabel="EROI Mono-Si PV")
    # plotData("../WindPotentialScala/EROI", 3,'EROI_CSP',xLabel="EROI CSPPT-12h TES")
    # plotData('solar',3,'',xLabel="Solar PV",subplot_size=2,subplot_index=1)
    # plotData('solar',4,'solar_csp_pv',xLabel="Solar CSP",subplot_size=2,subplot_index=2,save=True)
    
     plt.show()

def plotData(csvFile, index, output, xLabel="", save=False, subplot=False, subplot_size = 1, subplot_index = 1):
    
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
    
    if subplot:
        plt.subplot(subplot_size,1,subplot_index)
       
    #### later in the defined map 
    map = Basemap(projection='cyl', llcrnrlat=min(lats), urcrnrlat=max(lats), llcrnrlon=min(lons), urcrnrlon=max(lons), resolution='l')
    
    map.drawcoastlines()
    map.drawstates()
    plt.xlabel(xLabel)
    map.drawcountries()
    # map.drawlsmask(land_color='coral',ocean_color='blue') 
    print "Maximum ", max(values), " - Minimum ", min(values)
    if(max(values)==0):
        return
    
    cs = map.contourf(lons, lats, values,
                      np.linspace(min(values)+0.01, max(values)+0.01, 
                                  250, endpoint=True), tri=True ) #, latlon=True)
                      #np.linspace(0.01, 3, 250, endpoint=True), tri=True ) #, latlon=True
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label(xLabel)
    cbar.set_ticks(np.arange(0,10,1))
    
    if save :
        # High resolution
        #plt.savefig(output + '.pdf', dpi=250, bbox_inches='tight')
        plt.savefig(output + '.png', dpi=250, bbox_inches='tight')
        plt.close()
        
    return

    
def plotDataPerTech(index, label_short, label, eroi_list):
    #techs = ('Wind-onshore','Wind-offshore','PT-oil','PT-salt-TES','ST-salt-TES','Poly-Si-PV','Mono-Si-PV')
    techs = 'Wind' # Wind-onshore','Wind-offshore','ST-salt-TES','Mono-Si-PV')
    #for t in techs:
    for e in eroi_list:
        plotData('wind_solar_grid/'+techs+"_"+e,index, output='wind_solar_images/'+techs+"_"+e+'_'+label_short, xLabel=label, show=False)
    
    
if __name__ == "__main__":
    sys.exit(main())
