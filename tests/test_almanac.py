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

from skyfield.earthlib import refraction, refract

LATITUDE = 50.0
LONGITUDE = 13.0
ALTITUDE = 754
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
    'Skyfield': {
      'update_interval':0
    }
  }
})

default_formatter = weewx.units.get_default_formatter()

srv = skyfieldalmanac.SkyfieldService(None,CONFIG)
del weewx.almanac.almanacs[-1]
print(weewx.almanac.almanacs)

while skyfieldalmanac.sun_and_planets is None:
    time.sleep(1)


class AlmanacTest(unittest.TestCase):

    def setUp(self):
        os.environ['TZ'] = 'Europe/Berlin'
        time.tzset()
        self.alm = weewx.almanac.Almanac(TIME_TS,LATITUDE,LONGITUDE,altitude=169,temperature=TEMPERATURE_C,pressure=PRESSURE_MBAR,horizon=0.0,formatter=default_formatter)
    
    def test_sidereal_time(self):
        self.assertAlmostEqual(self.alm.sidereal_time,317.79217269677076,5)
        self.assertEqual(str(self.alm.sidereal_angle),'318°')
    
    def test_refraction(self):
        observer, horizon, body = skyfieldalmanac._get_observer(alm,'sun',False)
        #self.assertAlmostEqual(refr,0.5660773571020249,3)
        #self.assertAlmostEqual(refract(-refr,TEMPERATURE_C,PRESSURE_MBAR),0,3)
        self.assertAlmostEqual(horizon,-0.8928867722171352)
        
    def test_sun(self):
        # Test backwards compatibility
        self.assertEqual(str(self.alm.sunrise),'07:28:22')
        self.assertEqual(str(self.alm.sunset),'17:16:39')
        
        # Use Skyfield
        self.assertEqual(str(self.alm.sun.rise),'07:28:22')
        self.assertEqual(str(self.alm.sun.transit),'12:22:10')
        self.assertEqual(str(self.alm.sun.set),'17:16:39')
        self.assertEqual(str(self.alm.sun.antitransit),'00:22:09')
        self.assertEqual(str(self.alm.sun.previous_antitransit),'00:22:09')
        self.assertEqual(str(self.alm.sun.next_antitransit),'00:22:10')
        
        # Equinox / solstice
        self.assertEqual(str(self.alm.next_vernal_equinox),'03/20/25 10:01:28')
        self.assertEqual(str(self.alm.next_autumnal_equinox),'09/22/25 20:19:20')
        self.assertEqual(str(self.alm.next_summer_solstice),'06/21/25 04:42:15')
        self.assertEqual(str(self.alm.next_winter_solstice),'12/21/25 16:03:05')
        self.assertEqual(str(self.alm.previous_winter_solstice),'12/21/24 10:20:34')
        
        # Altitude / azimuth
        self.assertEqual(str(self.alm.sun.altitude),'25°')
        self.assertEqual(str(self.alm.sun.azimuth),'174°')
        self.assertEqual(str(self.alm.sun.alt_distance),'91714331.2 miles')
        self.assertAlmostEqual(self.alm.sun.alt,25.333616717531466,2)
        self.assertAlmostEqual(self.alm.sun.az,174.0628847017934,2)
        self.assertAlmostEqual(self.alm.sun.alt_dist,147599908.7,0)
        
        # Topocentric right ascension / declination / distance
        self.assertAlmostEqual(self.alm.sun.ra,323.33557380612825,3)
        self.assertAlmostEqual(self.alm.sun.dec,-14.516113980648846,3)
        self.assertAlmostEqual(self.alm.sun.dist,147599908.7,0)
        self.assertEqual(str(self.alm.sun.topo_ra),'323°')
        self.assertEqual(str(self.alm.sun.topo_dec),'-15°')
        self.assertEqual(str(self.alm.sun.topo_dist),'91714331.2 miles')
        
        # Hour angle / declination / distance
        self.assertAlmostEqual(self.alm.sun.ha,-5.5434011060793855,3)
        self.assertAlmostEqual(self.alm.sun.ha_dec,-14.516113980648846,3)
        self.assertAlmostEqual(self.alm.sun.ha_dist,147599908.7,0)
        self.assertEqual(str(self.alm.sun.hour_angle),'-06°')
        self.assertEqual(str(self.alm.sun.ha_declination),'-15°')
        self.assertEqual(str(self.alm.sun.ha_distance),'91714331.2 miles')
        
        # Astrometric right ascension / declination / distance
        self.assertAlmostEqual(self.alm.sun.a_ra,323.34100033196574,3)
        self.assertAlmostEqual(self.alm.sun.a_dec,-14.512050426594605,3)
        self.assertAlmostEqual(self.alm.sun.a_dist,147602648.27668718,0)
        self.assertEqual(str(self.alm.sun.astro_ra),'323°')
        self.assertEqual(str(self.alm.sun.astro_dec),'-15°')
        self.assertEqual(str(self.alm.sun.astro_dist),'91716033.5 miles')
        
        # Apparent astrometric right ascension / declination / distance
        self.assertAlmostEqual(self.alm.sun.g_ra,323.33535575318166,3)
        self.assertAlmostEqual(self.alm.sun.g_dec,-14.513890168137353,3)
        self.assertAlmostEqual(self.alm.sun.g_dist,147602648.27668718,0)
        self.assertEqual(str(self.alm.sun.geo_ra),'323°')
        self.assertEqual(str(self.alm.sun.geo_dec),'-15°')
        self.assertEqual(str(self.alm.sun.geo_dist),'91716033.5 miles')
    
    def test_moon(self):
        # rising and setting
        self.assertEqual(str(alm.moon.rise),'13:18:47')
        self.assertEqual(str(alm.moon.set),'06:13:17')
        
        # phase
        self.assertEqual(alm.moon_phase,'waxing gibbous (increasing to full)')
        self.assertEqual(alm.moon_index,3)
        self.assertEqual(alm.moon_fullness,90)
        self.assertAlmostEqual(alm.moon.moon_fullness,89.82441559659536,2)
        
        # events
        self.assertEqual(str(alm.next_new_moon),'02/28/25 01:44:49')
        self.assertEqual(str(alm.next_full_moon),'02/12/25 14:53:23')
        self.assertEqual(str(alm.previous_full_moon),'01/13/25 23:26:54')
        self.assertEqual(str(alm.next_first_quarter_moon),'03/06/25 17:31:37')
        self.assertEqual(str(alm.next_last_quarter_moon),'02/20/25 18:32:32')
    
    def test_exceptions(self):
        with self.assertRaises(AttributeError):
            self.alm.sun.foo
        with self.assertRaises(AttributeError):
            self.alm.foo.rise
        
alm = weewx.almanac.Almanac(TIME_TS,LATITUDE,LONGITUDE,altitude=169,temperature=TEMPERATURE_C,pressure=PRESSURE_MBAR,horizon=0.0,formatter=default_formatter)

#print(str(alm.sunrise))
#print(alm.sunset)
#print(alm.sun.altitude,alm.sun.azimuth,alm.sun.alt,alm.sun.az)
#print(alm.sun.visible)
#print(alm.sun.visible_change())
print(skyfieldalmanac.sun_and_planets)
print(alm.jupiter.rise)

if __name__ == '__main__':
    unittest.main()

srv.shutDown()
