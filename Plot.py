import sys
import pygrib
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
import csv
import numpy as np
import scipy.interpolate
from matplotlib.mlab import griddata

from numpy import genfromtxt

def plotData(csvFile, index):
    data = genfromtxt('ressources/results/' + csvFile, delimiter='\t', dtype=None)
    #### data preparation 
    lats = data[:, 0] 
    # # lon => x 
    lons = data[:, 1] 
    # # values => z 
    values = data[:, index] 
    print "Lat Min :", min(lats), "Lat Max :", max(lats)
    print "Lon Min :", min(lons), "Lon Max :", max(lons)
    #### later in the defined map 
    map = Basemap(projection='cyl', llcrnrlat=min(lats), urcrnrlat=max(lats), llcrnrlon=min(lons), urcrnrlon=max(lons), resolution='l')
    
    map.drawcoastlines()
    map.drawstates()
    map.drawcountries()
    print "Maximum Speed ", max(values)
    #clevs = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14]
    # clevs = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
    # clevs = [0,1000,1250,1500,2000,2500,3000,8000]
    clevs = [44,45]
    cs = map.contourf(lons, lats, values,clevs, latlon=True, tri=True)
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label('m/s')
    plt.show()
    return

