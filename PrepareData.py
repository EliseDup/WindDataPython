##
## From a Grib File with detailed observation about U and V wind components, 
## do all the transformation
## 1) Calculate min wind speed
## 2) For each grid cell associate a land cover class
from MeanWindSpeed import computeMeanWind
from LandCover import computeLandCover
from Plot import plotData
import sys

name = 'worldDaily2002'

def main():
    prepareData(name)
    # plotData(name+"final",6)

def prepareData(grib):
    computeMeanWind(grib,True,False)
    computeLandCover(grib)
    
if __name__ == "__main__":
    sys.exit(main())