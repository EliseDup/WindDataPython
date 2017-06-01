import os, sys, time, gdal
from gdalconst import *
from osgeo import gdal
from osgeo import osr
from numpy import zeros

def main():
    detailedLandCover("../WindPotential/grid_1deg", "grid1deg_25x25", 1.0, 25)

# Function to associate 3 different land cover types to a grid cell of given center (latitude,longitude)
def computeLandCover(fromFile, toFile, resolution):
    clc = '../resources/landCover/clc/g100_06.tif'
    glc = '../resources/landCover/modis/LCType.tif'
    globCover = '../resources/landCover/globCover/GLOBCOVER_L4_200901_200912_V2.3.tif'
    
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    dsCLC = gdal.Open(clc, GA_ReadOnly)
    dsGLC = gdal.Open(glc, GA_ReadOnly)
    dsGlob = gdal.Open(globCover, GA_ReadOnly)
    # get image size
    rowsCLC = dsCLC.RasterYSize; rowsGLC = dsGLC.RasterYSize; rowsGlob = dsGlob.RasterYSize;
    colsCLC = dsCLC.RasterXSize; colsGLC = dsGLC.RasterXSize; colsGlob = dsGlob.RasterXSize;

    # coordinates to get pixel values for
    res = open('results/'+toFile, 'w')

    # LAT - LON - Uwind - Vwind - meanWind - CLC - GLC
    with  open('results/'+fromFile) as f:
        for line in f:
            values = line.split("\t")
            lat = float(values[0]); long = float(values[1]);
            # First try CLC
            pixel = latLonToPixel(dsCLC, [lat, long], True) 
            if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < colsCLC and pixel[1] < rowsCLC): 
                clc = value(dsCLC, pixel)
            else: clc = "NA"
        
            # Then GLC
            pixelGlobal = latLonToPixel(dsGLC, [lat, long], False)
            if(pixelGlobal[0] >= 0 and pixelGlobal[1] >= 0 and pixelGlobal[0] < colsGLC and pixelGlobal[1] < rowsGLC): 
                glc = value(dsGLC, pixelGlobal)
            else: glc = "NA"
        
           # Then GLC
            pixelGlob = latLonToPixel(dsGlob, [lat, long], False)
            if(pixelGlob[0] >= 0 and pixelGlob[1] >= 0 and pixelGlob[0] < colsGlob and pixelGlob[1] < rowsGlob): 
                glob = value(dsGlob, pixelGlob)
            else: glob = "NA"
            
            # The resolution is 10 arc seconds ! So a square of 0.5 x 0.5 degrees = 180 x 180 10ArcSeconds !
            nUrban = 0
            
            if(glob != "NA"  and glob != 210):
              latStart = lat-resolution
              lonStart = long-resolution
              arcsec = resolution / 15
              for i in range(0,15):
                  for j in range(0,15):
                      latIt = latStart + i*arcsec
                      lonIt = lonStart + j*arcsec
                      urban = 0
                      pix = latLonToPixel(dsGlob, [latIt, lonIt], False)
                      if(pix[0] >= 0 and pix[1] >= 0 and pix[0] < colsGlob and pix[1] < rowsGlob): 
                          urban = value(dsGlob, pix)
                      if(urban==190): nUrban = nUrban+1
                            
              if(nUrban > 0):
                  print "Urban ?", lat , long, nUrban, "(", glob, ")"         
             
            if(clc != "NA" or glc != "NA" or glob != "NA"):
                res.write(line.split("\n")[0]);res.write("\t");
                res.write(str(clc));res.write("\t");
                res.write(str(glob));res.write("\t");
                res.write(str(glc));res.write("\t");
                res.write(str(nUrban));res.write("\n");
                
            
    res.close()

def detailedLandCover(fromFile, toFile, resolution, n):
    globCover = 'Globcover/GLOBCOVER_L4_200901_200912_V2.3.tif'
    # register all of the drivers
    gdal.AllRegister()
    # open the image
    ds = gdal.Open(globCover, GA_ReadOnly)
    # get image size
    rows = ds.RasterYSize; cols = ds.RasterXSize
    #Dictionary :
    indexes = {11:0, 14:1, 20:2, 30:3, 40:4, 50:5, 60:6, 70:7, 90:8, 100:9, 110:10, 120:11, 130:12, 140:13, 150:14, 160:15, 170:16, 180:17, 190:18, 200:19,210:20,220:21,230:22}
    # coordinates to get pixel values for
    res = open(toFile, 'w')
    # LAT - LON - Uwind - Vwind - meanWind - CLC - GLC
    pts = 0
    with  open(fromFile) as f:
        for line in f:
            pts = pts +1
            if pts%1000 == 0: print "Progress",pts
            values = line.split(",")
            lat = float(values[0]); long = float(values[1]);
            latStart = lat-resolution
            lonStart = long-resolution
            arcsec = resolution / n
            lcs = zeros(23, int)
            for i in range(0,n):
                for j in range(0,n):
                    latIt = latStart + i*arcsec
                    lonIt = lonStart + j*arcsec
                    pixel = latLonToPixel(ds, [latIt, lonIt], False)
                    if(pixel[0] >= 0 and pixel[1] >= 0 and pixel[0] < cols and pixel[1] < rows): 
                        clc = value(ds, pixel)
                        lcs[indexes[clc]] += 1
                        
            res.write(str(lat) + "\t" + str(long))
            for i in lcs:
                res.write("\t" + str(i))
            res.write("\n")
    
    res.close()
    
# The following method translates given latitude/longitude pairs into pixel locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#      latLonPairs - The decimal lat/lon pairings to be translated in the form [[lat1,lon1],[lat2,lon2]]
# OUTPUT: The pixel translation of the lat/lon pairings in the form [[x1,y1],[x2,y2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough 
#      image resolution for pixel size to be insignificant
def latLonToPixel(ds, latLon, lcl):
    # Load the image dataset
    # ds = gdal.Open(geotifAddr)
    # Get a geo-transform of the dataset
    gt = ds.GetGeoTransform()
    # Create a spatial reference object for the dataset
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    # Set up the coordinate transformation object
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srsLatLong, srs)
    # Change the point locations into the GeoTransform space
    # It seems there is an error with that method in case the DS is already defined in degrees we do not need this !!
    if(lcl): (latLon[1], latLon[0], holder) = ct.TransformPoint(latLon[1], latLon[0])
    # Translate the x and y coordinates into pixel values
    x = (latLon[1] - gt[0]) / gt[1]
    y = (latLon[0] - gt[3]) / gt[5]
    return [int(x), int(y)]

# The following method translates given pixel locations into latitude/longitude locations on a given GEOTIF
# INPUTS: geotifAddr - The file location of the GEOTIF
#      pixelPairs - The pixel pairings to be translated in the form [[x1,y1],[x2,y2]]
# OUTPUT: The lat/lon translation of the pixel pairings in the form [[lat1,lon1],[lat2,lon2]]
# NOTE:   This method does not take into account pixel size and assumes a high enough 
#      image resolution for pixel size to be insignificant
def pixelToLatLon(ds, pixel):
    # Get a geo-transform of the dataset
    gt = ds.GetGeoTransform()
    # Create a spatial reference object for the dataset
    srs = osr.SpatialReference()
    srs.ImportFromWkt(ds.GetProjection())
    # Set up the coordinate transformation object
    srsLatLong = srs.CloneGeogCS()
    ct = osr.CoordinateTransformation(srs, srsLatLong)
    # Translate the pixel pairs into untranslated points
    ulon = pixel[0] * gt[1] + gt[0]
    ulat = pixel[1] * gt[5] + gt[3]
    # Transform the points to the space
    (lon, lat, holder) = ct.TransformPoint(ulon, ulat)
    return [lat, lon]

def value(ds, pixel):
    band = ds.GetRasterBand(1)  # 1-based index
    # read data and add the value to the string
    data = band.ReadAsArray(pixel[0], pixel[1], 1, 1)
    return data[0, 0]

def printInfo(ds):
    # get georeference info
    transform = ds.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]
    print 'Transform =', transform
    print 'Origin = (', transform[0], ',', transform[3], ')'
    print 'Pixel Size = (', transform[1], ',', transform[5], ')'
    
if __name__ == "__main__":
    sys.exit(main())