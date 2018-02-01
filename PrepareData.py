##
## From a Grib File with detailed observation about U and V wind components, 
## do all the transformation
## 1) Calculate min wind speed
## 2) For each grid cell associate a land cover class
import os, sys, time, gdal, pygrib
import numpy as np

#from MeanWindSpeed import computeMeanWind, computeMeanData, computeMeanDissipation
from LandCover import computeLandCover, detailedLandCover
from Plot import plotData
from gdalconst import *
from osgeo import gdal
from osgeo import osr
from LandCover import latLonToPixel
from LandCover import value

def main():
    computeValueFromTiff('grid', 'grid_sea','../resources/elevation/seaLevel.tif', 0.75, 15)
    computeValueFromTiff('grid', 'grid_coast','../resources/coast_distance/coastDistance.tif', 0.75, 15)
   
def gribToFile(name):
    grib_data = pygrib.open(name +'.grib')
    data, lats, lons = grib_data.message(1).data()
    res = open(name, 'w')
    for i in range(0, len(lats)):
        for j in range(0, len(lats[0])):
            res.write(str(lats[i][j]));res.write("\t");
            res.write(str(lons[i][j]));res.write("\t");
            res.write(str(data[i][j]));res.write("\n");
            
    res.close()
    return
    
def computeValueFromTiff(fromFile, toFile, tiffFile, resolution, n):
    # sea = 'ressources/seaLevel.tif'
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    ds = gdal.Open(tiffFile, GA_ReadOnly)
    rows = ds.RasterYSize;
    cols = ds.RasterXSize;
    res = open(toFile, 'w')
    with  open(fromFile) as f:
        for line in f:
            values = line.split("\t")
            lat = float(values[0]); long = float(values[1]);
            # print lat, " , ", long
            latStart = lat-resolution
            lonStart = long-resolution
            arcsec = resolution / n
            mean = 0
            nObs = 0
            for i in range(0,n):
                for j in range(0,n):
                    latIt = latStart + i*arcsec
                    lonIt = lonStart + j*arcsec
                    pixel = latLonToPixel(ds, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows): 
                        mean = mean + value(ds, pixel)
                        nObs = nObs + 1
            
            res.write(line.split("\n")[0]);res.write("\t");
            res.write(str(mean));res.write("\t");res.write(str(nObs));res.write("\n");
                
            
    res.close()  

if __name__ == "__main__":
    sys.exit(main())
 

