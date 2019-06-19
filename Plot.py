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
     #plotData("../results_server/netESimple_total_kepersquaremeter", 2+2, "wind", xLabel="Wind", save = True, subplot=False)
     #plotData("../results_server/netESimple_total_kepersquaremeter", 5+2, "pv", xLabel="Solar PV", save = True, subplot=False)
     #plotData("../results_server/netESimple_total_kepersquaremeter", 8+2, "csp", xLabel="Solar CSP", save = True, subplot=False)
    
     #plotData("results_eco/results_simple_x_cell_total_10_05", 3, "max_net_e", xLabel="Wind",subplot_size=3,subplot_index=1)
     #plotData("results_eco/results_simple_x_cell_total_10_05", 4, "max_net_e",xLabel="PV",subplot_size=3,subplot_index=2)
     #plotData("results_eco/results_simple_x_cell_total_10_05", 5, "max_net_e",xLabel="CSP",subplot_size=3, subplot_index=3)

     plt.show()

def plotData(csvFile, index, output, xLabel="", save=True, subplot=True, subplot_size = 1, subplot_index = 1):
    
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
    # plt.xlabel(xLabel)
    map.drawcountries()
    # map.drawlsmask(land_color='coral',ocean_color='blue') 
    print "Maximum ", max(values), " - Minimum ", min(values)
    if(max(values)==0):
        return
    
    cs = map.contourf(lons, lats, values,np.linspace(min(values)+0.01, max(values)+0.01, 
                                  250, endpoint=True), 
                                  tri=True ) #, latlon=True)
                      #np.linspace(0.01, 3, 250, endpoint=True), tri=True ) #, latlon=True
    
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label(xLabel)
    
    tick= round(max(values),1)/5
    if tick==0: tick = max(values)/5
    
    #cbar.set_ticks(np.arange(0,max(values), tick))
    
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
