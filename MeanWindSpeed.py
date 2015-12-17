import sys
import pygrib
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
from ncepgrib2 import dump

VERBOSE = 1  # verbose error reporting

# From observations for several years, we would like to compute the mean values of wind speed for each location
# and write that in a new grib file
# We have to have a grib file with only value for U and V wind component
def computeMeanWind(fileName,write,plot):
    grib_data = pygrib.open('ressources/emcwf/'+fileName +'.grib')
    uComponents = grib_data.select(name='10 metre U wind component')
    vComponents = grib_data.select(name='10 metre V wind component')
    
    # Geographical coordinates are the same for all the observations
    lats, lons = uComponents[0].latlons()
    meanWindSpeed = ndarray((len(lats), len(lats[1])), float)
    meanU = ndarray((len(lats), len(lats[1])), float)
    meanV = ndarray((len(lats), len(lats[1])), float)
    
    nObs = ndarray((len(lats), len(lats[1])), float)
    print(grib_data.messages)
    print len(uComponents),"-",len(vComponents)
    for m in range(0, len(uComponents)):
        print uComponents[m] , " ", vComponents[m]
        dataU, latsU, lonsU = uComponents[m].data(); dataV, latsV, lonsV = vComponents[m].data()
        for i in range(0, len(dataU)):
            # j corresponds to the geo point
            for j in range(0, len(dataU[i])):
                # Wind Speed = sqrt(u^2 + v^2)
                meanWindSpeed[i][j] = meanWindSpeed[i][j] + math.sqrt(dataU[i][j]*dataU[i][j]+dataV[i][j]*dataV[i][j])
                meanU[i][j] = meanU[i][j] + dataU[i][j]
                meanV[i][j] = meanV[i][j] + dataV[i][j]
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            meanWindSpeed[i][j] = meanWindSpeed[i][j] / nObs[i][j]
            meanU[i][j] = meanU[i][j] / nObs[i][j]
            meanV[i][j] = meanV[i][j] / nObs[i][j]
            # print(str(meanWindSpeed[i][j]) + " " + str(lats[i][j]) + " " +str( lons[i][j]) )
    if(write): writeData(fileName,[meanU,meanV,meanWindSpeed], lats, lons);
    if(plot): plotData(meanWindSpeed, lats, lons);
    
    return
  
# From observations for several years, we would like to compute the mean values of wind speed for each location
def writeData(fileName,data,lats,lons):
    res = open(fileName, 'w')
    nData = len(data)
    print(str(len(data))+"-"+str(len(lats))+"-"+str(len(lons)))
    print(str(len(data[0]))+"-"+str(len(lats[0]))+"-"+str(len(lons[0])))
    for i in range(0, len(lats)):
        for j in range(0, len(lats[0])):
            res.write(str(lats[i][j]));res.write("\t");
            res.write(str(lons[i][j]));res.write("\t")
            for n in range(0,nData-1):
                res.write(str(data[n][i][j]));res.write("\t");
            res.write(str(data[nData-1][i][j]));res.write("\n");
            
    res.close()
    return

def plotData(data,lats,lons):
    # create polar stereographic Basemap instance.
    map = Basemap(projection='cyl',llcrnrlat=min(lats),urcrnrlat=max(lats),\
            llcrnrlon=min(lons),urcrnrlon=max(lons),resolution='l')
    
    # map.drawmapboundary(fill_color='aqua')# draw coastlines, state and country boundaries, edge of map.
    # map.fillcontinents(color='coral',lake_color='aqua')
    map.drawcoastlines()
    map.drawstates()
    map.drawcountries()
    #clevs = [0,0.1,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11]
    cs = map.contourf(lats,lons,data)
    cbar = map.colorbar(cs,location='bottom',pad="5%")
    cbar.set_label('m/s')
    plt.show()
    return

