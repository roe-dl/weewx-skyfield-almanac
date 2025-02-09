#!/usr/bin/python3

import locale
import os
import unittest

locale.setlocale(locale.LC_ALL, 'C')
os.environ['LANG'] = 'C'

import time
import sys
import configobj
sys.path.append('/usr/share/weewx')
sys.path.append('bin/user')
import weewx
import weewx.almanac
import weewx.units
import skyfieldalmanac

LATITUDE = 51.123
LONGITUDE = 13.040
ALTITUDE = 169
TIME_TS = 1739098800
TEMPERATURE_C = 15.0
PRESSURE_MBAR = 1013.25
CONFIG = configobj.ConfigObj({
  'DatabaseTypes': {
    'SQLite': {
      'SQLITE_ROOT':'tests'
    }
  },
  'Almanac': {
    'update_interval':0
  }
})

default_formatter = weewx.units.get_default_formatter()

srv = skyfieldalmanac.SkyfieldService(None,CONFIG)
del weewx.almanac.almanacs[-1]
print(weewx.almanac.almanacs)

while skyfieldalmanac.eph is None:
    time.sleep(1)


class AlmanacTest(unittest.TestCase):

    def setUp(self):
        os.environ['TZ'] = 'Europe/Berlin'
        time.tzset()
        self.alm = weewx.almanac.Almanac(TIME_TS,LATITUDE,LONGITUDE,altitude=169,temperature=TEMPERATURE_C,pressure=PRESSURE_MBAR,horizon=0.0,formatter=default_formatter)
    
    def test_sidereal_time(self):
        self.assertAlmostEqual(self.alm.sidereal_time,317.8321726967708,5)
        self.assertEqual(str(self.alm.sidereal_angle),'318°')
    
    def test_sun(self):
        # Test backwards compatibility
        self.assertEqual(str(self.alm.sunrise),'07:33:20')
        self.assertEqual(str(self.alm.sunset),'17:11:23')
        
        # Use Skyfield
        self.assertEqual(str(self.alm.sun.rise),'07:33:20')
        self.assertEqual(str(self.alm.sun.transit),'12:22:00')
        self.assertEqual(str(self.alm.sun.set),'17:11:23')
        
        # Equinox / solstice
        self.assertEqual(str(self.alm.next_vernal_equinox),'03/20/25 10:01:28')
        self.assertEqual(str(self.alm.next_autumnal_equinox),'09/22/25 20:19:20')
        self.assertEqual(str(self.alm.next_summer_solstice),'06/21/25 04:42:15')
        self.assertEqual(str(self.alm.next_winter_solstice),'12/21/25 16:03:05')
        self.assertEqual(str(self.alm.previous_winter_solstice),'12/21/24 10:20:34')
        
        self.assertEqual(str(self.alm.sun.altitude),'24°')
        self.assertEqual(str(self.alm.sun.azimuth),'174°')
        self.assertAlmostEqual(self.alm.sun.alt,24.220897454265796,2)
        self.assertAlmostEqual(self.alm.sun.az,174.1585381089835,2)
    
    def test_moon(self):
        self.assertEqual(alm.moon_phase,'waxing gibbous (increasing to full)')
        self.assertEqual(alm.moon_index,3)
        self.assertEqual(alm.moon_fullness,90)
        self.assertAlmostEqual(alm.moon.moon_fullness,89.82441559659536,2)
        
        self.assertEqual(str(alm.next_new_moon),'02/28/25 01:44:49')
        self.assertEqual(str(alm.next_full_moon),'02/12/25 14:53:23')
        self.assertEqual(str(alm.previous_full_moon),'01/13/25 23:26:54')
        self.assertEqual(str(alm.next_first_quarter_moon),'03/06/25 17:31:37')
        self.assertEqual(str(alm.next_last_quarter_moon),'02/20/25 18:32:32')
    
    def test_exceptions(self):
        with self.assertRaises(weewx.UnknownType):
            self.alm.sun.foo
        
alm = weewx.almanac.Almanac(TIME_TS,LATITUDE,LONGITUDE,altitude=169,temperature=TEMPERATURE_C,pressure=PRESSURE_MBAR,horizon=0.0,formatter=default_formatter)

#print(str(alm.sunrise))
#print(alm.sunset)
#print(alm.sun.altitude,alm.sun.azimuth,alm.sun.alt,alm.sun.az)
print(alm.sun.visible)
print(alm.sun.visible_change())

if __name__ == '__main__':
    unittest.main()

srv.shutDown()
