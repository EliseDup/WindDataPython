#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
 
server = ECMWFDataServer()
server.retrieve({
    "class": "ei",
    "dataset": "interim",
    "date": "2015-08-01/to/2018-07-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "area": "73.5/-27/33/45",
    "levelist": "57/58",
    "levtype": "ml",
    "param": "131.128/132.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
    "target": "WindEU-2015-2018.grib",
})