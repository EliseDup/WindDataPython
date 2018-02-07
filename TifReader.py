import os, sys, time, gdal
from gdalconst import *
from osgeo import gdal
from osgeo import osr
from numpy import zeros
import numpy as np
import LandCover as lc
import struct

def main():
   computeIrradiance(360,"solar_1deg")
       
def readFile(name,globCover,output):
    res = open(output, 'w')
    ds = gdal.Open(name, GA_ReadOnly)
    cover = gdal.Open(globCover, GA_ReadOnly)
    lc.printInfo(ds)
    rows = ds.RasterYSize; cols = ds.RasterXSize;
    arcres = 1200
    for lat in range(-55,60):
        for i in range(0,3600/arcres):
            latIt = lat+i*arcres/3600.0;
            for lon in range(-180,180):
                for j in range(0,3600/arcres):
                    lonIt = lon+i*arcres/3600.0;
                    pixel = lc.latLonToPixel(ds, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ds, pixel)
                        if(not np.isnan(value)):
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(latIt) + "\t" + str(lonIt) + "\t" + str(value) + "\t" + str(valueCover) + "\n")
                            
                        #else : res.write(str(latIt) + "\t" + str(lonIt) + "\t" + "0.0" +"\n")
def computeIrradiance(arcsec, output):
    file_ghi = '../resources/solar/GHI/GHI.tif'
    file_dni = '../resources/solar/DNI/DNI.tif'
    globCover = '../resources/landCover/globCover/GLOBCOVER_L4_200901_200912_V2.3.tif'
   
    # register all of the drivers
    gdal.AllRegister()
   
    res = open(output, 'w')
    ghi = gdal.Open(file_ghi, GA_ReadOnly); dni = gdal.Open(file_dni, GA_ReadOnly); cover = gdal.Open(globCover, GA_ReadOnly)
    rows = ghi.RasterYSize; cols = ghi.RasterXSize;

    for lat in range(-55,60):
        for i in range(0,3600/arcsec):
            latIt = lat+i*arcsec/3600.0;
            for lon in range(-180,180):
                for j in range(0,3600/arcsec):
                    lonIt = lon+i*arcsec/3600.0;
                    pixel = lc.latLonToPixel(ghi, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ghi, pixel)
                        if(not np.isnan(value)):
                            direct = lc.value(dni, pixel)
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(latIt) + "\t" + str(lonIt) + "\t" + str(value) + "\t" + str(direct) + "\t" + str(valueCover) + "\n")
                            


if __name__ == "__main__":
    sys.exit(main())