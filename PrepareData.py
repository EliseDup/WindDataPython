##
## From a Grib File with detailed observation about U and V wind components, 
## do all the transformation
## 1) Calculate min wind speed
## 2) For each grid cell associate a land cover class
import os, sys, time, gdal, pygrib
import numpy as np

from MeanWindSpeed import computeMeanWind, computeMeanData, computeMeanDissipation
from LandCover import computeLandCover, detailedLandCover
from Plot import plotData
from gdalconst import *
from osgeo import gdal
from osgeo import osr
from LandCover import latLonToPixel
from LandCover import value


name = '../test'
    
def main():
    gribToFile('dec2016')
    #detailedLandCover('../test.txt', '../res.txt', 0.75, 15)
    #computeMeanData(name, 'U component of wind', 'V component of wind') 
    #computeLandCover(name, name+'lc', 0.5)
    #computeValueFromTiff(name+'lc', name+'sea','../resources/elevation/seaLevel.tif')
    #computeValueFromTiff(name+'sea', name+'coast','../resources/coastDistance.tif')
   
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
    
def computeValueFromTiff(fromFile, toFile, tiffFile):
    # sea = 'ressources/seaLevel.tif'
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    ds = gdal.Open(tiffFile, GA_ReadOnly)
    rows = ds.RasterYSize;
    cols = ds.RasterXSize;
    res = open('results/'+toFile, 'w')
    with  open('results/'+fromFile) as f:
        for line in f:
            values = line.split("\t")
            lat = float(values[0]); long = float(values[1]);
            # First try CLC
            pixel = latLonToPixel(ds, [lat, long], True) 
            if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows): 
                level = value(ds, pixel)
            else: level = "NA"
        
            if(level != "NA"):
                res.write(line.split("\n")[0]);res.write("\t");
                res.write(str(level));res.write("\n");
                
            
    res.close()  

if __name__ == "__main__":
    sys.exit(main())
 

