import os, sys, time, gdal
from gdalconst import *
from osgeo import gdal
from osgeo import osr
from numpy import zeros
import numpy as np
import LandCover as lc
import struct
from numpy import genfromtxt

# The fastest way to have the data from R (extraction from shp files)
# and from python (extraction from tif files), is to first launch this code to have the grid
# (i.e. all the couple (latitude,longitude) that are relevant (with value != 0), then
# compute the value from the shp file with the grid (ex :Countrie.R), then run again this code with the correct countries.
# (In scala or R directly the concatenation is taking forever)
file_ghi = '../resources/solar/GHI/GHI.tif'
file_dni = '../resources/solar/DNI/DNI.tif'
globCover = '../resources/landCover/globCover/GLOBCOVER_L4_200901_200912_V2.3.tif'
coast_distance  = "../resources/coast_distance/coastDistance.tif"
elevation = "../resources/slope/elev_0_0833.tif"
slope_path = "../resources/slope/cl"

protectedPath= "../resources/protected_areas/"
countriesPath= "../resources/countries/"

def main():
    print 'Hello'
    #computeGrid(1800,"test_total")
    #addRResultsToGrid("test_total", countriesPath+"countries_0_5deg_total", "test_total_countries")
    addRResultsToGrid("test_total_countries", protectedPath+"protected_0_5deg_total", "test_total_protected",true_false=True)
    addSlopeClasses("test_total_protected","test_total_slopes")
    
def computeGrid(arcres, finalOutput):
    readFile(file_ghi, "temp1", arcres)
    addToFile("temp1", file_dni,"temp")
    addToFile("temp", globCover,"temp1")
    addToFile("temp1", coast_distance,"temp")
    addToFile("temp",elevation,finalOutput)
   
def addSlopeClasses(gridFile, output): 
    addToFile(gridFile, slope_path+"1_0_0833.tif", "temp1")
    for i in range(2,8):
        j = i-1
        print str(j), "TO", str(i)
        addToFile("temp"+str(j), slope_path+str(i)+"_0_0833.tif", "temp"+str(i))
    addToFile("temp7", slope_path+"8_0_0833.tif", output)
    
# Read a tif file and print a file with 3 columns : latitude, longitude, and the corresponding value. The grid is computed with the resolution given in arcsecond !
def readFile(dataset, output, arcres = 1800, factor = 1.0, all = True):
    res = open(output, 'w')
    ds = gdal.Open(dataset, GA_ReadOnly)
    lc.printInfo(ds)
    rows = ds.RasterYSize; cols = ds.RasterXSize;
    
    for lat in range(-90, 90):
        for i in range(0, 3600 / arcres):
            latIt = lat + i * arcres / 3600.0;
            for lon in range(-180, 180):
                for j in range(0, 3600 / arcres):
                    lonIt = lon + j * arcres / 3600.0;
                    pixel = lc.latLonToPixel(ds, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ds, pixel)
                        if(not np.isnan(value) and float(value)!=-9999): res.write(str(lonIt) + "\t" + str(latIt) + "\t" +str(float(factor*value)) + "\n")
                        elif all: res.write(str(lonIt) + "\t" + str(latIt) + "\t"+ "0.0" + "\n")
    res.close()

# Grid is a txt file with lat & lon in the two first column and other information after
# Dataset is a tif file with the information to add a the end of the grid file
def addToFile(grid, dataset, output, factor=1.0):
    outputFile = open(output, 'w')
    ds = gdal.Open(dataset, GA_ReadOnly)
    lc.printInfo(ds)
    rows = ds.RasterYSize; cols = ds.RasterXSize;
    with  open(grid) as f:
        for line in f:
            values = line.split("\t")
            outputFile.write(line.split("\n")[0]);outputFile.write("\t");
            lon = float(values[0]);lat = float(values[1]); 
            pixel = lc.latLonToPixel(ds, [lat, lon], False)
            if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                value = lc.value(ds, pixel)
                if(not np.isnan(value) and float(value)!=-9999): outputFile.write(str(factor*value) + "\n")
                else: outputFile.write('0.0' + "\n")
            else: outputFile.write('0.0' + "\n")
    f.close(); outputFile.close();
    
def readFileWithCover(name, globCover, output, arcres=1200):
    res = open(output, 'w')
    ds = gdal.Open(name, GA_ReadOnly)
    cover = gdal.Open(globCover, GA_ReadOnly)
    lc.printInfo(ds)
    rows = ds.RasterYSize; cols = ds.RasterXSize;
    for lat in range(-55, 60):
        for i in range(0, 3600 / arcres):
            latIt = lat + i * arcres / 3600.0;
            for lon in range(-180, 180):
                for j in range(0, 3600 / arcres):
                    lonIt = lon + j * arcres / 3600.0;
                    pixel = lc.latLonToPixel(ds, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ds, pixel)
                        if(not np.isnan(value)):
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(lonIt) + "\t" + str(latIt) + "\t" + str(value) + "\t" + str(valueCover) + "\n")
                        # else : 
                        #    pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                        #    valueCover = lc.value(cover, pixelCover)
                        #    res.write(str(lonIt) + "\t" + str(latIt) + "\t" + "0.0" + "\t" + str(valueCover) + "\n")

def addRResultsToGrid(grid, rFile, output, with_nan = False, true_false = False):
    outputFile = open(output, 'w')
    rData = genfromtxt(rFile, delimiter="\t", dtype=str)
    print(len(rData))
    k = 0;
    with  open(grid) as f:
        for line in f:
            outputFile.write(line.split("\n")[0]);outputFile.write("\t");
            if(with_nan and np.isnan(rData[k])): res ="NA"
            elif(true_false): 
                if(rData[k]=="FALSE"):res=0
                else: res=1
            else: res = rData[k].replace('"', '')
            outputFile.write(str(res) + "\n")
            k = k +1
            
    outputFile.close(); f.close();
    if(k != len(rData)): print "Incompatible Sizes", " -- ", len(rData), " != ", k

                                  
def computeIrradiance(arcsec, output, country, countryFile, all, protected, protectedFile):
    
    file_ghi = '../resources/solar/GHI/GHI.tif'
    file_dni = '../resources/solar/DNI/DNI.tif'
    globCover = '../resources/landCover/globCover/GLOBCOVER_L4_200901_200912_V2.3.tif'
   
    # register all of the drivers
    gdal.AllRegister()
   
    res = open(output, 'w')
    ghi = gdal.Open(file_ghi, GA_ReadOnly); dni = gdal.Open(file_dni, GA_ReadOnly); cover = gdal.Open(globCover, GA_ReadOnly)
    rows = ghi.RasterYSize; cols = ghi.RasterXSize;

    if country:
        countries = genfromtxt(countryFile, delimiter="\t", dtype=str)
        print(len(countries))
    if protected:
        prot = genfromtxt(protectedFile, delimiter="\t", dtype=str)
        print(len(prot))
    kcountry = 0; kprotected = 0;   
    
    for lat in range(-55, 60):
        if(lat % 10 == 0): print "Progress", lat
        for i in range(0, 3600 / arcsec):
            latIt = lat + i * arcsec / 3600.0;
            for lon in range(-180, 180):
                for j in range(0, 3600 / arcsec):
                    lonIt = lon + j * arcsec / 3600.0;
                    pixel = lc.latLonToPixel(ghi, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                        value = lc.value(ghi, pixel)
                        
                        if(not np.isnan(value)):
                            direct = lc.value(dni, pixel)
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(lonIt) + "\t" + str(latIt) + "\t" + str(value) + "\t" + str(direct) + "\t" + str(valueCover))
                        elif all:
                            pixelCover = lc.latLonToPixel(cover, [latIt, lonIt], False)
                            valueCover = lc.value(cover, pixelCover)
                            res.write(str(lonIt) + "\t" + str(latIt) + "\t" + "0.0" + "\t" + "0.0" + "\t" + str(valueCover))
                       
                        if country:
                            res.write("\t" + str(countries[kcountry].replace('"', '')))
                            kcountry = kcountry + 1
                        if protected:
                            res.write("\t" + str(prot[kprotected]))
                            kprotected = kprotected + 1
                        
                        res.write("\n")
                  
# # Here I just added the proportion of the area with slope < 5, i.e. slope clases 1,2 and 3      
def addSlopeToGrid(grid, output):
    outputFile = open(output, 'w')

 #   cl1 = gdal.Open("../resources/slope/cl1_0_5.tif", GA_ReadOnly)
 #   cl2 = gdal.Open("../resources/slope/cl2_0_5.tif", GA_ReadOnly)
 #   cl3 = gdal.Open("../resources/slope/cl3_0_5.tif", GA_ReadOnly)
    cl4 = gdal.Open("../resources/slope/cl4_0_0833.tif", GA_ReadOnly)
    cl5 = gdal.Open("../resources/slope/cl5_0_0833.tif", GA_ReadOnly)
    cl6 = gdal.Open("../resources/slope/cl6_0_0833.tif", GA_ReadOnly)
    cl7 = gdal.Open("../resources/slope/cl7_0_0833.tif", GA_ReadOnly)
    cl8 = gdal.Open("../resources/slope/cl8_0_0833.tif", GA_ReadOnly)
    
    lc.printInfo(cl4)
    rows = cl4.RasterYSize; cols = cl4.RasterXSize;
    with  open(grid) as f:
        for line in f:
            outputFile.write(line.split("\n")[0]);outputFile.write("\t");
            values = line.split("\t")
            lat = float(values[1]); lon = float(values[0]);
            val = 0.0
            pixel = lc.latLonToPixel(cl4, [lat, lon], False)
            if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                #val1 = lc.value(cl1, pixel); val2 = lc.value(cl2, pixel); val3 = lc.value(cl3, pixel); 
                val4 = lc.value(cl4, pixel); val5 = lc.value(cl5, pixel); val6 = lc.value(cl6, pixel); val7 = lc.value(cl7, pixel);  val8 = lc.value(cl8, pixel); 
                
                #if(val1 != -9999): val = val + val1 / 1000.0
                #if(val2 != -9999): val = val + val2 / 1000.0
                #if(val3 != -9999): val = val + val3 / 1000.0
                if(val4 != -9999): val = val + val4 / 100.0;
                if(val5 != -9999): val = val + val5 / 100.0;
                if(val6 != -9999): val = val + val6 / 100.0;
                if(val7 != -9999): val = val + val7 / 100.0;
                if(val8 != -9999): val = val + val8 / 100.0;
            outputFile.write(str(val) + "\n")
    f.close();                   
    outputFile.close()
    
def checkSlopeToGrid(grid, output, factor = 1.0/100.0):
    outputFile = open(output, 'w')

    cl1 = gdal.Open("../resources/slope/cl1_0_0833.tif", GA_ReadOnly)
    cl2 = gdal.Open("../resources/slope/cl2_0_0833.tif", GA_ReadOnly)
    cl3 = gdal.Open("../resources/slope/cl3_0_0833.tif", GA_ReadOnly)
    cl4 = gdal.Open("../resources/slope/cl4_0_0833.tif", GA_ReadOnly)
    cl5 = gdal.Open("../resources/slope/cl5_0_0833.tif", GA_ReadOnly)
    cl6 = gdal.Open("../resources/slope/cl6_0_0833.tif", GA_ReadOnly)
    cl7 = gdal.Open("../resources/slope/cl7_0_0833.tif", GA_ReadOnly)
    cl8 = gdal.Open("../resources/slope/cl8_0_0833.tif", GA_ReadOnly)
    
    lc.printInfo(cl1)
    rows = cl1.RasterYSize; cols = cl1.RasterXSize;
    with  open(grid) as f:
        for line in f:
            #outputFile.write(line.split("\n")[0]);outputFile.write("\t");
            values = line.split("\t")
            lat = float(values[0]); lon = float(values[1]);
            val = 0.0
            pixel = lc.latLonToPixel(cl1, [lat, lon], False)
            if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows):
                val1 = lc.value(cl1, pixel); val2 = lc.value(cl2, pixel); val3 = lc.value(cl3, pixel); 
                val4 = lc.value(cl4, pixel); val5 = lc.value(cl5, pixel); val6 = lc.value(cl6, pixel); 
                val7 = lc.value(cl7, pixel); val8 = lc.value(cl8, pixel); 
                
                if(val1 != -9999): val = val + val1 * factor
                if(val2 != -9999): val = val + val2 * factor
                if(val3 != -9999): val = val + val3 * factor
                if(val4 != -9999): val = val + val4 * factor
                if(val5 != -9999): val = val + val5 * factor
                if(val6 != -9999): val = val + val6 * factor
                if(val7 != -9999): val = val + val7 * factor
                if(val8 != -9999): val = val + val8 * factor
            
            if(val!=0): outputFile.write(str(lat) + "\t" + str(lon) + "\t" +str(val) + "\n")
    f.close();                   
    outputFile.close()
                 
if __name__ == "__main__":
    sys.exit(main())
