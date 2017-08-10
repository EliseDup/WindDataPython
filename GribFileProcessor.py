import sys
import pygrib
from numpy import ndarray
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import math
from ncepgrib2 import dump

VERBOSE = 1  # verbose error reporting
windU = '10 metre U wind component'
windV = '10 metre V wind component'
solarDwn = 'Surface solar radiation downwards'
solarNet = 'Surface net solar radiation'
solarNetClearSky = 'Surface net solar radiation, clear sky'
eraFolder = '../resources/era_interim_data/'

def main():
    computeMeanWind('europe_level57', 57, True, False)
    #writeMeanByMonth(eraFolder+'solar/radiationDwn_40years', 'radiationDwn_40yearsByMonth', solarDwn, 24 * 60 * 60)
    #writeMean(eraFolder+'solar/radiationDwn_40years', 'radiationDwn_40years', solarDwn, 24 * 60 * 60)

def getComponent(gribFile, componentName):
    return pygrib.open(gribFile + '.grib').select(name=componentName)

def getLatLons(gribFile):
    return pygrib.open(gribFile + '.grib').message(1).latlons()

# From a list of messages in components, returns the mean for each latitude, longitude pair
def computeMean(components, lats, lons, factor=1.0):
    print 'Compute means for ', len(components), ' values'
    mean = ndarray((len(lats), len(lats[1])), float)
    nObs = ndarray((len(lats), len(lats[1])), float)
    for m in range(0, len(components)):
        if(m % 100 == 0): print components[m]
        data, lats, lons = components[m].data()
        for i in range(0, len(data)):
            # j corresponds to the geo point
            for j in range(0, len(data[i])):
                mean[i][j] = mean[i][j] + data[i][j] / factor
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            mean[i][j] = mean[i][j] / nObs[i][j]
    return mean

def computeMeanByMonth(components, lats, lons, factor=1.0):
    print 'Compute means by month for ', len(components), ' values'
    mean = ndarray((12, len(lats), len(lats[1])), float)
    nObs = ndarray((12, len(lats), len(lats[1])), float)
    month = 0
    for m in range(0, len(components)):
        if(m % 100 == 0): print components[m]
        data, lats, lons = components[m].data()
        for i in range(0, len(data)):
            # j corresponds to the geo point
            for j in range(0, len(data[i])):
                mean[month][i][j] = mean[month][i][j] + data[i][j] / factor
                nObs[month][i][j] = nObs[month][i][j] + 1
        if month == 11: month = 0
        else: month = month + 1
    # Do the actual mean
    for month in range(0, 12):
        for i in range(0, len(lats)):
            for j in range(0, len(lats[i])):
                mean[month][i][j] = mean[month][i][j] / nObs[month][i][j]
    return mean

def computeMeanUV(uComponents, vComponents, lats, lons, factor=1.0):
    print 'Compute means (SQRT(u^2 + v^2)) for ', len(uComponents), ' values'
    mean = ndarray((len(lats), len(lats[1])), float)
    nObs = ndarray((len(lats), len(lats[1])), float)
    for m in range(0, min(len(vComponents), len(uComponents))):
        if(m % 100 == 0): print uComponents[m] , " ", vComponents[m]
        dataU, latsU, lonsU = uComponents[m].data(); dataV, latsV, lonsV = vComponents[m].data()
        for i in range(0, len(dataU)):
            # j corresponds to the geo point
            for j in range(0, len(dataU[i])):
                # Value = sqrt(u^2 + v^2)
                mean[i][j] = mean[i][j] + math.sqrt(dataU[i][j] / factor * dataU[i][j] / factor + dataV[i][j] / factor * dataV[i][j] / factor)
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            mean[i][j] = mean[i][j] / nObs[i][j]
            
    standardDeviation = ndarray((len(lats), len(lats[1])), float)
    for m in range(0, min(len(vComponents), len(uComponents))):
        if(m % 100 == 0): print 'Standard deviation for - ', uComponents[m] , " ", vComponents[m]
        dataU, latsU, lonsU = uComponents[m].data(); dataV, latsV, lonsV = vComponents[m].data()
        for i in range(0, len(dataU)):
            # j corresponds to the geo point
            for j in range(0, len(dataU[i])):
                error = (math.sqrt(dataU[i][j] / factor * dataU[i][j] / factor + dataV[i][j] / factor * dataV[i][j] / factor) - mean[i][j])
                standardDeviation[i][j] = standardDeviation[i][j] + error * error
    
    return mean, standardDeviation, nObs

def writeMean(input, output, componentName, factor=1.0):
    data = getComponent(input, componentName)
    lats, lons = getLatLons(input)
    mean = computeMean(data, lats, lons, factor)
    writeData(output, [mean], lats, lons);
    return

def writeMeanUV(input, output, uComponentName, vComponentName, factor=1.0):
    dataU = getComponent(input, uComponentName);dataV = getComponent(input, vComponentName)
    lats, lons = getLatLons(input)
    mean, standardDeviation, nObs = computeMeanUV(dataU, dataV, lats, lons, factor)
    writeData(output, [mean, standardDeviation, nObs], lats, lons);
    return

def writeMeanByMonth(input, output, componentName, factor=1.0):
    data = getComponent(input, componentName)
    lats, lons = getLatLons(input)
    # mean = computeMean(data, lats, lons, factor)
    meanByMonth = computeMeanByMonth(data, lats, lons, factor)
    writeData(output, meanByMonth, lats, lons);
    return

# From observations for several years, we would like to compute the mean values of wind speed for each location
def writeData(fileName, data, lats, lons):
    res = open(fileName, 'w')
    nData = len(data)
    print(str(len(data)) + "-" + str(len(lats)) + "-" + str(len(lons)))
    print(str(len(data[0])) + "-" + str(len(lats[0])) + "-" + str(len(lons[0])))
    for i in range(0, len(lats)):
        for j in range(0, len(lats[0])):
            res.write(str(lats[i][j]));res.write("\t");
            res.write(str(lons[i][j]));res.write("\t")
            for n in range(0, nData - 1):
                res.write(str(data[n][i][j]));res.write("\t");
            res.write(str(data[nData - 1][i][j]));res.write("\n");
            
    res.close()
    return

def plotData(data, lats, lons):
    map = Basemap(projection='cyl', llcrnrlat=min(lats), urcrnrlat=max(lats), llcrnrlon=min(lons), urcrnrlon=max(lons), resolution='l')
    map.drawcoastlines()
    map.drawstates()
    map.drawcountries()
    print "Maximum ", max(data), " - Minimum ", min(data)
    cs = map.contourf(lons, lats, data, np.linspace(min(data), max(data), 250, endpoint=True), tri=True)  # , latlon=True)
    cbar = map.colorbar(cs, location='bottom', pad="5%")
    cbar.set_label('')
    plt.show()
    return




# From observations for several years, we would like to compute the mean values of wind speed for each location
# and write that in a new grib file
# We have to have a grib file with only value for U and V wind component
def computeMeanWind(outputName, level=57, write=True, plot=False):
    # files = ('2007_2008', '2009_2010', '2011_2012', '2013_2014', '2015_2016')
    files = ('europe2012','europe2013','europe2014','europe2015','europe2016') # ('2012','2013','2014','2015','2016'),'europe2014','europe2014'
    folder = ('/Users/Elise/Desktop/wind_eu/')
    nFiles = 5
    uComponents = list()
    vComponents = list()
    
    for f in files:
        grib_data = pygrib.open(folder + f + '.grib')
        uComponents.append(grib_data.select(level=level, name='U component of wind'))  # To the EAST
        vComponents.append(grib_data.select(level=level, name='V component of wind'))  # To the NORTH
        
        print 'Size', f, '-', grib_data.messages
    
    # Geographical coordinates are the same for all the observations
    lats, lons = uComponents[0][0].latlons()
       
    meanWindSpeed = ndarray((len(lats), len(lats[1])), float)
    meanU = ndarray((len(lats), len(lats[1])), float)
    meanV = ndarray((len(lats), len(lats[1])), float)
    
    nObs = ndarray((len(lats), len(lats[1])), float)
    for f in range(0, nFiles):
        for m in range(0, min(len(vComponents[f]), len(uComponents[f]))):
            if(m % 100 == 0): print uComponents[f][m] , " ", vComponents[f][m]
            dataU, latsU, lonsU = uComponents[f][m].data(); dataV, latsV, lonsV = vComponents[f][m].data()
            for i in range(0, len(dataU)):
                # j corresponds to the geo point
                for j in range(0, len(dataU[i])):
                    # Wind Speed = sqrt(u^2 + v^2)
                    meanWindSpeed[i][j] = meanWindSpeed[i][j] + math.sqrt(dataU[i][j] * dataU[i][j] + dataV[i][j] * dataV[i][j])
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
    for f in range(0, nFiles):
        for m in range(0, min(len(vComponents[f]), len(uComponents[f]))):
            if(m % 100 == 0): print uComponents[f][m] , " ", vComponents[f][m]
            dataU, latsU, lonsU = uComponents[f][m].data(); dataV, latsV, lonsV = vComponents[f][m].data()
            for i in range(0, len(dataU)):
                 # j corresponds to the geo point
                 for j in range(0, len(dataU[i])):
                     error = (math.sqrt(dataU[i][j] * dataU[i][j] + dataV[i][j] * dataV[i][j]) - meanWindSpeed[i][j])
                     standardDeviation[i][j] = standardDeviation[i][j] + error * error
      # Do the actual mean
    # for i in range(0, 10):#len(lats)):
    #    for j in range(0, 10):#len(lats[i])):
    #            standardDeviation[i][j] = math.sqrt(standardDeviation[i][j] / (nObs[i][j]))
                
    if(write): writeData('results/' + outputName, [meanU, meanV, meanWindSpeed, standardDeviation, nObs], lats, lons);
    if(plot): plotData(meanWindSpeed, lats, lons);
    
    return

def computeMeanDissipation(fileName, write, plot):
    grib_data = pygrib.open(fileName + '.grib')
    
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
    print len(uComponents), "-", len(vComponents)
    for m in range(0, min(len(vComponents), len(uComponents))):
        if(m % 100 == 0):print uComponents[m] , " ", uStressComponents[m]
        windU, latsU, lonsU = uComponents[m].data(); windV, latsV, lonsV = vComponents[m].data()
        stressU, latsU1, lonsU1 = uStressComponents[m].data(); stressV, latsV1, lonsV1 = vStressComponents[m].data()
        for i in range(0, len(windU)):
            # j corresponds to the geo point
            for j in range(0, len(windU[i])):
                # Dissipation = tau . v = SQRT(tu*u + tv*v)
                # It should be 86400 (= 24 h ) and not 8760 !
                instantDissipation = math.sqrt(max(0, windU[i][j] * (stressU[i][j] / (24 * 60 * 60)) + windV[i][j] * (stressV[i][j] / (24 * 60 * 60))))
                meanDissipation[i][j] = meanDissipation[i][j] + instantDissipation
                nObs[i][j] = nObs[i][j] + 1
    # Do the actual mean
    for i in range(0, len(lats)):
        for j in range(0, len(lats[i])):
            # print lats[i][j], nObs[i][j]
            meanDissipation[i][j] = meanDissipation[i][j] / nObs[i][j]

            # print(str(meanWindSpeed[i][j]) + " " + str(lats[i][j]) + " " +str( lons[i][j]) )
    if(write): writeData(fileName, [meanDissipation], lats, lons);
    if(plot): plotData(meanDissipation, lats, lons);
    
    return

if __name__ == "__main__":
    sys.exit(main())
 
