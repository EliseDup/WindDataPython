import numpy as np
import sys
import linecache

name = "../resources/slope/GloSlopesCl3_5min.asc"
res = 0.5

def main():
    ncols = int(linecache.getline(name,1).split()[1])
    nrows = int(linecache.getline(name,2).split()[1])
    xllcorner = int(linecache.getline(name,3).split()[1])
    yllcorner = int(linecache.getline(name,4).split()[1])
    cellsize = float(linecache.getline(name,5).split()[1])
    nodata = linecache.getline(name,6)
    output = open('cl3_0_5deg','w')
    
    step = int(res/cellsize)
    print step
    print 180/res
    print 360/res
    
    for y in range(0, int(180/res)):
        row = linecache.getline(name,y*step+7).split()
        for x in range(0, int(360/res)):
            lat = yllcorner+y*res
            lon = xllcorner+x*res
            
            output.write(str(lat) + "\t" + str(lon) + "\t" + str(float(row[x*step])/1000) + "\n")
        
    output.close()

if __name__ == "__main__":
    sys.exit(main())