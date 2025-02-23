#!/usr/bin/python3

import locale
import os
import unittest

import time
import sys
import configobj
sys.path.append('/usr/share/weewx')
sys.path.append('bin/user')
import weewx
import weewx.almanac
import weewx.units
import skyfieldalmanac

LATITUDE = 50.0
LONGITUDE = 15.0
ALTITUDE = 169
TIME_TS = 1761951600
START_OF_YEAR = 1735729200
TEMPERATURE_C = 15.0
PRESSURE_MBAR = 1013.25
CONFIG = configobj.ConfigObj({
  'DatabaseTypes': {
    'SQLite': {
      'SQLITE_ROOT':'tests'
    }
  },
  'Almanac': {
    'Skyfield': {
      'update_interval':0
    }
  }
})

default_formatter = weewx.units.get_default_formatter()

srv = skyfieldalmanac.SkyfieldService(None,CONFIG)
del weewx.almanac.almanacs[-1]
print(weewx.almanac.almanacs)

while skyfieldalmanac.eph is None:
    time.sleep(1)

print("day")
print(" h     ha       sun     diff  sundial")
for i in range(24):
    alm = weewx.almanac.Almanac(TIME_TS+i*3600,LATITUDE,LONGITUDE,altitude=169,temperature=TEMPERATURE_C,pressure=PRESSURE_MBAR,horizon=0.0,formatter=default_formatter)
    ha_vh = alm.sun.hour_angle
    ha = (ha_vh.raw-180.0)%360.0
    print("%2d  %5s  %8.4f  %.4f  %s" % (i,ha_vh,ha/15,(ha/15-i)*60,alm.solar_datetime.format("%d.%m.%Y %H:%M:%S")))
print("")


print("year")
print(time.strftime("%Y %H:%M:%S",time.localtime(START_OF_YEAR)))
print("day     ha       sun      diff   sundial              altitude    azimuth")
for i in range(365):
    alm = weewx.almanac.Almanac(START_OF_YEAR+i*86400,LATITUDE,LONGITUDE,altitude=169,temperature=TEMPERATURE_C,pressure=PRESSURE_MBAR,horizon=0.0,formatter=default_formatter)
    ha_vh = alm.sun.hour_angle
    ha = (ha_vh.raw-180.0)%360.0
    print("%3d  %5s  %8.4f  %8.4f   %17s   %s   %s" % (i+1,ha_vh,ha/15,(ha/15-12)*60,alm.solar_datetime.format("%d.%m.%Y %H:%M:%S"),alm.sun.altitude.format("%6.3f",add_label=False),alm.sun.azimuth.format("%7.3f",add_label=False)))

srv.shutDown()
