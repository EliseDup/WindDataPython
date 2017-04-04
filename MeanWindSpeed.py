import sys
import pygrib
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
from ncepgrib2 import dump

VERBOSE = 1  # verbose error reporting

def main():
    computeMeanData('../resources/wind/level58_2016','level58_2016', 'U component of wind','V component of wind')
    #computeMeanData('../dissipation40years','Instantaneous eastward turbulent surface stress','Instantaneous northward turbulent surface stress')
    #computeMeanWind('test', True, False)
# From observations for several years, we would like to compute the mean values of wind speed for each location
# and write that in a new grib file
# We have to have a grib file with only value for U and V wind component
def computeMeanWind(fileName,write,plot):
    files = ('2007_2008','2009_2010','2011_2012','2013_2014','2015_2016')
    nFiles = 5
    uComponents = list()
    vComponents= list()
    
    for f in files:
        grib_data = pygrib.open('../resources/wind/modelsLevel58/'+f+'.grib')
        uComponents.append(grib_data.select(name='U component of wind')) # To the EAST
        vComponents.append(grib_data.select(name='V component of wind')) # To the NORTH
        
        print 'Size',f,'-',grib_data.messages
    
    # Geographical coordinates are the same for all the observations
    lats, lons = uComponents[0][0].latlons()
       
    meanWindSpeed = ndarray((len(lats), len(lats[1])), float)
    meanU = ndarray((len(lats), len(lats[1])), float)
    meanV = ndarray((len(lats), len(lats[1])), float)
    
    nObs = ndarray((len(lats), len(lats[1])), float)
    for f in range(0,nFiles):
        for m in range(0, min(len(vComponents[f]),len(uComponents[f]))):
            if(m%100 == 0): print uComponents[f][m] , " ", vComponents[f][m]
            dataU, latsU, lonsU = uComponents[f][m].data(); dataV, latsV, lonsV = vComponents[f][m].data()
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
            # print lats[i][j], nObs[i][j]
            meanWindSpeed[i][j] = meanWindSpeed[i][j] / nObs[i][j]
            meanU[i][j] = meanU[i][j] / nObs[i][j]
            meanV[i][j] = meanV[i][j] / nObs[i][j]
            # print(str(meanWindSpeed[i][j]) + " " + str(lats[i][j]) + " " +str( lons[i][j]) )
    
    standardDeviation = ndarray((len(lats), len(lats[1])), float)
    for f in range(0,nFiles):
        for m in range(0, min(len(vComponents[f]),len(uComponents[f]))):
            if(m%100 == 0): print uComponents[f][m] , " ", vComponents[f][m]
            dataU, latsU, lonsU = uComponents[f][m].data(); dataV, latsV, lonsV = vComponents[f][m].data()
            for i in range(0, len(dataU)):
                 #j corresponds to the geo point
                 for j in range(0, len(dataU[i])):
                     error = (math.sqrt(dataU[i][j]*dataU[i][j]+dataV[i][j]*dataV[i][j]) - meanWindSpeed[i][j])
                     standardDeviation[i][j] = standardDeviation[i][j] + error*error
      # Do the actual mean
    #for i in range(0, 10):#len(lats)):
    #    for j in range(0, 10):#len(lats[i])):
    #            standardDeviation[i][j] = math.sqrt(standardDeviation[i][j] / (nObs[i][j]))
                
    if(write): writeData('results/'+fileName,[meanU,meanV,meanWindSpeed,standardDeviation,nObs], lats, lons);
    if(plot): plotData(meanWindSpeed, lats, lons);
    
    return

def computeMeanDissipation(fileName,write,plot):
    grib_data = pygrib.open(fileName+'.grib')
    
    uComponents = grib_data.select(name='10 metre U wind component')
    vComponents = grib_data.select(name='10 metre V wind component')
        
    uStressComponents = grib_data.select(name='Eastward turbulent surface stress')
    vStressComponents = grib_data.select(name='Northward turbulent surface stress')
    
    print len(uComponents), len(vComponents), len(uStressComponents), len(vStressComponents)
    
    # Geographical coordinates are the same for all the observations
    lats, lons = uComponents[0].latlons()
    meanDissipation = ndarray((len(lats), len(lats[1])), float)
    
    nObs = ndarray((len(lats), len(lats[1])), float)
    print(grib_data.messages)
    print len(uComponents),"-",len(vComponents)
    for m in range(0, min(len(vComponents),len(uComponents))):
        if(m%100 == 0):print uComponents[m] , " ", uStressComponents[m]
        windU, latsU, lonsU = uComponents[m].data(); windV, latsV, lonsV = vComponents[m].data()
        stressU, latsU1, lonsU1 = uStressComponents[m].data(); stressV, latsV1, lonsV1 = vStressComponents[m].data()
        for i in range(0, len(windU)):
            # j corresponds to the geo point
            for j in range(0, len(windU[i])):
                # Dissipation = tau . v = SQRT(tu*u + tv*v)
                # It should be 86400 (= 24 h ) and not 8760 !
                instantDissipation = math.sqrt(max(0,windU[i][j]*(stressU[i][j]/(24*60*60))+windV[i][j]*(stressV[i][j]/(24*60*60))))
                meanDissipation[i][j] = meanDissipation[i][j] + instantDissipation
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            # print lats[i][j], nObs[i][j]
            meanDissipation[i][j] = meanDissipation[i][j] / nObs[i][j]

            # print(str(meanWindSpeed[i][j]) + " " + str(lats[i][j]) + " " +str( lons[i][j]) )
    if(write): writeData(fileName,[meanDissipation], lats, lons);
    if(plot): plotData(meanDissipation, lats, lons);
    
    return

def computeMeanData(fileName,write,plot,componentName='Surface net solar radiation'):
    grib_data = pygrib.open(fileName +'.grib')
    
    components = grib_data.select(name=componentName)
    
    # Geographical coordinates are the same for all the observations
    lats, lons = components[0].latlons()
    meanWindSpeed = ndarray((len(lats), len(lats[1])), float)
    mean = ndarray((len(lats), len(lats[1])), float)
    
    nObs = ndarray((len(lats), len(lats[1])), float)
    print grib_data.messages, len(components)
    for m in range(0, len(components)):
        print components[m]
        data, lats, lons = components[m].data()
        for i in range(0, len(data)):
            # j corresponds to the geo point
            for j in range(0, len(data[i])):
                mean[i][j] = mean[i][j] + data[i][j]
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            # print lats[i][j], nObs[i][j]
            mean[i][j] = mean[i][j] / nObs[i][j]
  
    if(write): writeData(fileName,[mean], lats, lons);
    if(plot): plotData(mean, lats, lons);
    
    return

def computeMeanData(fileName, outputFileName, uComponentName, vComponentName):
    grib_data = pygrib.open(fileName+'.grib')
    
    uComponents = grib_data.select(name=uComponentName)
    vComponents = grib_data.select(name=vComponentName)
  
    print 'NVALUE',grib_data.messages
    
    # Geographical coordinates are the same for all the observations
    lats, lons = uComponents[0].latlons()
    meanData = ndarray((len(lats), len(lats[1])), float)
    meanU = ndarray((len(lats), len(lats[1])), float)
    meanV = ndarray((len(lats), len(lats[1])), float)
    
    nObs = ndarray((len(lats), len(lats[1])), float)
    print(grib_data.messages)
    print len(uComponents),"-",len(vComponents)
    for m in range(0, min(len(vComponents),len(uComponents))):
        if(m%100 == 0): print uComponents[m] , " ", vComponents[m]
        dataU, latsU, lonsU = uComponents[m].data(); dataV, latsV, lonsV = vComponents[m].data()
        for i in range(0, len(dataU)):
            # j corresponds to the geo point
            for j in range(0, len(dataU[i])):
                # Wind Value = sqrt(u^2 + v^2)
                meanData[i][j] = meanData[i][j] + math.sqrt(dataU[i][j]*dataU[i][j]+dataV[i][j]*dataV[i][j])
                meanU[i][j] = meanU[i][j] + dataU[i][j]
                meanV[i][j] = meanV[i][j] + dataV[i][j]
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            # print lats[i][j], nObs[i][j]
            meanData[i][j] = meanData[i][j] / nObs[i][j]
            meanU[i][j] = meanU[i][j] / nObs[i][j]
            meanV[i][j] = meanV[i][j] / nObs[i][j]
            # print(str(meanWindSpeed[i][j]) + " " + str(lats[i][j]) + " " +str( lons[i][j]) ) 
    
    standardDeviation = ndarray((len(lats), len(lats[1])), float)
    for m in range(0, min(len(vComponents),len(uComponents))):
        if(m%100 == 0): print uComponents[m] , " ", vComponents[m]
        dataU, latsU, lonsU = uComponents[m].data(); dataV, latsV, lonsV = vComponents[m].data()
        for i in range(0, len(dataU)):
             #j corresponds to the geo point
            for j in range(0, len(dataU[i])):
                error = (math.sqrt(dataU[i][j]*dataU[i][j]+dataV[i][j]*dataV[i][j]) - meanData[i][j])
                standardDeviation[i][j] = standardDeviation[i][j] + error*error
    
    #for i in range (0, len(dataU)):
    #    for j in range(0, len(dataU[i])):
    #        standardDeviation[i][j] = math.sqrt(standardDeviation[i][j]/nObs[i][j])
                
    writeData('results/'+outputFileName,[meanData,standardDeviation, nObs], lats, lons);
    
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
    map = Basemap(projection='cyl', llcrnrlat=min(lats), urcrnrlat=max(lats), llcrnrlon=min(lons), urcrnrlon=max(lons), resolution='l')
    map.drawcoastlines()
    map.drawstates()
    map.drawcountries()
    print "Maximum ", max(data), " - Minimum ", min(data)
    cs = map.contourf(lons, lats, data, np.linspace(min(data), max(data), 250, endpoint=True), tri=True ) #, latlon=True)
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label('')
    plt.show()
    return

if __name__ == "__main__":
    sys.exit(main())
 
