import os, sys, time, gdal
from gdalconst import *
from osgeo import gdal
from osgeo import osr
from numpy import zeros
import numpy as np
import LandCover as lc
import struct
from numpy import genfromtxt
def main():
    #computeIrradiance(3600,"solar_1deg")
    computeIrradiance(180,"test")
      
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
                    lonIt = lon+j*arcres/3600.0;
                    pixel = lc.latLonToPixel(ds, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ds, pixel)
                        if(not np.isnan(value)):
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(lonIt) + "\t" + str(latIt) + "\t" + str(value) + "\t" + str(valueCover) + "\n")
                        #else : 
                        #    pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                        #    valueCover = lc.value(cover, pixelCover)
                        #    res.write(str(lonIt) + "\t" + str(latIt) + "\t" + "0.0" + "\t" + str(valueCover) + "\n")
                            
def computeIrradiance(arcsec, output):
    file_ghi = '../resources/solar/GHI/GHI.tif'
    file_dni = '../resources/solar/DNI/DNI.tif'
    globCover = '../resources/landCover/globCover/GLOBCOVER_L4_200901_200912_V2.3.tif'
   
    # register all of the drivers
    gdal.AllRegister()
   
    res = open(output, 'w')
    ghi = gdal.Open(file_ghi, GA_ReadOnly); dni = gdal.Open(file_dni, GA_ReadOnly); cover = gdal.Open(globCover, GA_ReadOnly)
    rows = ghi.RasterYSize; cols = ghi.RasterXSize;

    countries = genfromtxt("../countries_0_05deg_3", delimiter="\t", dtype=str)
    print(len(countries))
    k = 0
    for lat in range(-55,60):
        if(lat%10 == 0): print "Progress", lat
        for i in range(0,3600/arcsec):
            latIt = lat+i*arcsec/3600.0;
            for lon in range(-180,180):
                for j in range(0,3600/arcsec):
                    lonIt = lon+j*arcsec/3600.0;
                    pixel = lc.latLonToPixel(ghi, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ghi, pixel)
                        if(not np.isnan(value)):
                            direct = lc.value(dni, pixel)
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(lonIt) + "\t" + str(latIt) + "\t" + str(value) + "\t" + str(direct) + "\t" + str(valueCover) + "\t" + str(countries[k].replace('"', '')) + "\n")
                            k = k+1

if __name__ == "__main__":
    sys.exit(main())