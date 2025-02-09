#!/usr/bin/python3
# Almanac extension to WeeWX using Skyfield
# Copyright (C) 2025 Johanna Roedenbeck

"""

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

"""
    documentation about ephemeris files:
    https://rhodesmill.org/skyfield/planets.html
    
    dates, times, and timescales:
    https://rhodesmill.org/skyfield/time.html
    
    Skyfield comes with builtin timescale data. They are updated with
    every new Skyfield version. But you can also download them.
    
    Timescale data are provided by IERS. The file is called `finals2000A.all`.
    
    
    Configuration:
    
    [Almanac]
        # which ephemeris file to use
        ephemeris = de440s.bsp  # or de440.bsp or de441.bsp
        # use builtin timescale data or download it from IERS
        use_builtin_timescale = true
        # update interval 1 year (set to 0 for no update)
        update_interval = 31557600
    
"""

VERSION = "0.1"

import time
import datetime
import configobj
import threading

import weewx
from weewx.almanac import AlmanacType, almanacs, timestamp_to_djd
from weewx.engine import StdService
from weewx.units import ValueTuple, ValueHelper
import weeutil

from skyfield import almanac
from skyfield.api import N, S, E, W, Loader, wgs84
from skyfield.earthlib import refraction

ts = None
eph = None


def timestamp_to_skyfield_time(timestamp):
    return ts.utc(1970,1,1+timestamp/86400.0)

def skyfield_time_to_djd(ti):
    return ti.ut1-ti.dut1/86400.0-2415020.0

def _get_observer(almanac_obj, time_ts):
    # Build an ephem Observer object
    observer = eph['Earth'] + wgs84.latlon(almanac_obj.lat,almanac_obj.lon,elevation_m=almanac_obj.altitude)
    refr = refraction(
        almanac_obj.horizon,
        temperature_C=almanac_obj.temperature,
        pressure_mbar=almanac_obj.pressure
    )
    return observer, refr


class SkyfieldAlmanacType(AlmanacType):
    """ Almanac extension to provide the date of the Easter sunday """
    
    EVENTS = {
        'vernal_equinox':(0,),
        'summer_solstice':(1,),
        'autumnal_equinox':(2,),
        'winter_solstice':(3,),
        'equinox':(0,2),
        'solstice':(1,3),
        'new_moon':(0,),
        'first_quarter_moon':(1,),
        'full_moon':(2,),
        'last_quarter_moon':(3,)
    }

    @property
    def hasExtras(self):
        """ PyEphem provides extras. """
        return True

    def get_almanac_data(self, almanac_obj, attr):
        """ calculate attribute """
        if ts is None or eph is None:
            raise weewx.UnknownType(attr)
        time_ti = timestamp_to_skyfield_time(almanac_obj.time_ts)
        if attr=='sunrise':
            return almanac_obj.sun.rise
        elif attr=='sunset':
            return almanac_obj.sun.set
        elif attr=='moon_fullness':
            return int(almanac_obj.moon.moon_fullness + 0.5)
        elif attr in ('moon_phase','moon_index'):
            position = almanac.moon_phase(eph, time_ti).degrees/360.0
            moon_index = int((position * 8) + 0.5) & 7
            if attr=='moon_index': return moon_index
            return almanac_obj.moon_phases[moon_index]
        elif attr in {'previous_equinox', 'next_equinox',
                      'previous_solstice', 'next_solstice',
                      'previous_autumnal_equinox', 'next_autumnal_equinox',
                      'previous_vernal_equinox', 'next_vernal_equinox',
                      'previous_winter_solstice', 'next_winter_solstice',
                      'previous_summer_solstice', 'next_summer_solstice',
                      'previous_new_moon', 'next_new_moon',
                      'previous_first_quarter_moon', 'next_first_quarter_moon',
                      'previous_full_moon', 'next_full_moon',
                      'previous_last_quarter_moon', 'next_last_quarter_moon'}:
            # This is how you call a function on an instance when all you have
            # is the function's name as a string
            previous = attr.startswith('previous')
            if attr.endswith('_moon'):
                # moon phases
                interval = 2592000
                func = almanac.moon_phases
            else:
                # seasons
                interval = 31557600
                func = almanac.seasons
            if previous:
                # looking for events before the given timestamp ("previous")
                t0 = -interval
                t1 = 0
                x = attr[9:]
            else:
                # looking for events after the given timestamp ("next")
                t0 = 0
                t1 = interval
                x = attr[5:]
            # get the Skyfield event codes
            event = SkyfieldAlmanacType.EVENTS.get(x)
            if event is None: raise weewx.UnknownType(attr)
            # time interval to look for events
            t0 = timestamp_to_skyfield_time(almanac_obj.time_ts+t0)
            t1 = timestamp_to_skyfield_time(almanac_obj.time_ts+t1)
            # find the events
            t, y  = almanac.find_discrete(t0,t1,func(eph))
            # in case of previous events search from the last event on
            if previous:
                t = reversed(t)
                y = reversed(y)
            # look for the event
            for ti, yi in zip(t,y):
                if yi in event:
                    djd = skyfield_time_to_djd(ti)
                    break
            return weewx.units.ValueHelper(ValueTuple(djd, "dublin_jd", "group_time"),
                                           context="ephem_year",
                                           formatter=almanac_obj.formatter,
                                           converter=almanac_obj.converter)
        # Check to see if the attribute is a sidereal angle
        elif attr == 'sidereal_time' or attr == 'sidereal_angle':
            # Local Apparent Sidereal Time (LAST)
            # sidereal time is obtained from an ephem Observer object...
            observer = wgs84.latlon(almanac_obj.lat,almanac_obj.lon,elevation_m=almanac_obj.altitude)
            # ... then get the angle in degrees ...
            val = observer.lst_hours_at(time_ti)*15.0
            # ... finally, depending on the attribute name, pick the proper return type:
            if attr == 'sidereal_time':
                return val
            else:
                vt = ValueTuple(val, 'degree_compass', 'group_direction')
                return weewx.units.ValueHelper(vt,
                                               context = 'ephem_day',
                                               formatter=almanac_obj.formatter,
                                               converter=almanac_obj.converter)
        else:
            # The attribute must be a heavenly body (such as 'sun', or 'jupiter').
            # Bind the almanac and the heavenly body together and return as an
            # AlmanacBinder
            return SkyfieldAlmanacBinder(almanac_obj, attr)
        # `attr` is not provided by this extension. So raise an exception.
        raise weewx.UnknownType(attr)


class SkyfieldAlmanacBinder:
    """This class binds the observer properties held in Almanac, with the heavenly
    body to be observed."""

    def __init__(self, almanac, heavenly_body):
        self.almanac = almanac

        # Calculate and store the start-of-day in Dublin Julian Days. 
        y, m, d = time.localtime(self.almanac.time_ts)[0:3]
        #self.sod_djd = timestamp_to_djd(time.mktime((y, m, d, 0, 0, 0, 0, 0, -1))

        self.heavenly_body = heavenly_body
        self.use_center = False

    def __call__(self, use_center=False):
        self.use_center = use_center
        return self
    
    @property
    def visible(self):
        """Calculate how long the body has been visible today"""
        observer, refr = _get_observer(self.almanac,self.almanac.time_ts)
        body = eph[self.heavenly_body]
        timespan = weeutil.weeutil.archiveDaySpan(self.almanac.time_ts)
        t0 = timestamp_to_skyfield_time(timespan[0])
        t1 = timestamp_to_skyfield_time(timespan[1])
        tr, yr = almanac.find_risings(observer, body, t0, t1, horizon_degrees=-refr)
        ts, ys = almanac.find_settings(observer, body, t0, t1, horizon_degrees=-refr)
        if len(tr)<1 or len(ts)<1:
            visible = None
        elif yr[-1] and ys[-1]:
            visible = (ts.ut1-tr.ut1) * weewx.units.SECS_PER_DAY
        else:
            #TODO always up and always down
            visible = 0
        return weewx.units.ValueHelper(ValueTuple(visible, "second", "group_deltatime"),
                                       context="day",
                                       formatter=self.almanac.formatter,
                                       converter=self.almanac.converter)
        
    def visible_change(self, days_ago=1):
        """Change in visibility of the heavenly body compared to 'days_ago'.
           Copyright (C) Tom Keffer
        """
        # Visibility for today:
        today_visible = self.visible
        # The time to compare to
        then_time = self.almanac.time_ts - days_ago * 86400
        # Get a new almanac, set up for the time back then
        then_almanac = self.almanac(almanac_time=then_time)
        # Find the visibility back then
        then_visible = getattr(then_almanac, self.heavenly_body).visible
        # Take the difference
        diff = today_visible.raw - then_visible.raw
        return weewx.units.ValueHelper(ValueTuple(diff, "second", "group_deltatime"),
                                       context="hour",
                                       formatter=self.almanac.formatter,
                                       converter=self.almanac.converter)

    def __getattr__(self, attr):
        """Get the requested observation, such as when the body will rise."""

        # Don't try any attributes that start with a double underscore, or any of these
        # special names: they are used by the Python language:
        if attr.startswith('__') or attr in ['mro', 'im_func', 'func_code']:
            raise AttributeError(attr)
        
        observer, refr = _get_observer(self.almanac,self.almanac.time_ts)
        body = eph[self.heavenly_body]
        
        previous = attr.startswith('previous_')
        next = attr.startswith('next_')
        
        if previous:
            # get the last event before the given timestamp
            t0 = timestamp_to_skyfield_time(self.almanac.time_ts-86400)
            t1 = timestamp_to_skyfield_time(self.almanac.time_ts)
            evt = attr[9:]
        elif next:
            # get the next event after the given timestamp
            t0 = timestamp_to_skyfield_time(self.almanac.time_ts)
            t1 = timestamp_to_skyfield_time(self.almanac.time_ts+86400)
            evt = attr[5:]
        elif attr in ('rise','set','transit','antitransit'):
            # get the event within the day the timestamp is in
            timespan = weeutil.weeutil.archiveDaySpan(self.almanac.time_ts)
            t0 = timestamp_to_skyfield_time(timespan[0])
            t1 = timestamp_to_skyfield_time(timespan[1])
            evt = attr
        else:
            # convert given timestamp
            ti = timestamp_to_skyfield_time(self.almanac.time_ts)
            position = observer.at(ti).observe(body).apparent()
            if attr=='moon_fullness':
                return position.fraction_illuminated(eph['sun'])*100.0
            if attr in ('az','alt','azimuth','altitude'):
                alt, az, distance = position.altaz(temperature_C=self.almanac.temperature,pressure_mbar=self.almanac.pressure)
                if attr=='az':
                    return az.degrees
                elif attr=='alt':
                    return alt.degrees
                else:
                    if attr=='azimuth':
                        vt = ValueTuple(az.degrees,'degree_compass','group_direction')
                    elif attr=='altitude':
                        vt = ValueTuple(alt.radians,'radian','group_angle')
                    return ValueHelper(vt,
                                               context="ephem_day",
                                               formatter=self.almanac.formatter,
                                               converter=self.almanac.converter)

            # `attr` is not provided by this extension. So raise an exception.
            raise weewx.UnknownType("%s.%s" % (self.heavenly_body,attr))

        # Note: In case of polar day or night, y is False and t is the time
        #       of transit or antitransit.
        
        t = None
        if evt in ('rise','rising'):
            t, y = almanac.find_risings(observer, body, t0, t1, horizon_degrees=-refr)
        elif evt in ('set','setting'):
            t, y = almanac.find_settings(observer, body, t0, t1, horizon_degrees=-refr)
        elif evt=='transit':
            t = almanac.find_transits(observer, body, t0, t1)
            y = True
        elif evt=='antitransit':
            raise weewx.UnknownType("%s.%s" % (self.heavenly_body,attr))
        else:
            # `attr` is not provided by this extension. So raise an exception.
            raise weewx.UnknownType(attr)
        if t is not None:
            time_djd = skyfield_time_to_djd(t[-1]) if len(t)>=1 and y else None
            return weewx.units.ValueHelper(ValueTuple(time_djd, "dublin_jd", "group_time"),
                                           context="ephem_day",
                                           formatter=self.almanac.formatter,
                                           converter=self.almanac.converter)


class SkyfieldAlmanacThread(threading.Thread):
    """ Thread to download and update ephemeris """
    
    def __init__(self, alm_conf_dict, path):
        """ init thread
        
            Args:
                alm_conf_dict(configobj): almanac configuration
                path(str): directory where to save downloaded files
        """
        super(SkyfieldAlmanacThread,self).__init__(name='SkyfieldThread')
        self.path = path
        self.eph_file = alm_conf_dict.get('ephemeris','de440s.bsp')
        self.builtin = weeutil.weeutil.to_bool(alm_conf_dict.get('use_builtin_timescale',True))
        self.update_interval = weeutil.weeutil.to_int(alm_conf_dict.get('update_interval',31557600))
        self.evt = threading.Event()
        self.running = True
    
    def shutDown(self):
        """ shut down thread """
        self.running = False
        self.evt.set()
    
    def run(self):
        """ """
        try:
            while self.running:
                self.init_skyfield()
                if not self.update_interval: break
                self.evt.wait(self.update_interval)
        except Exception as e:
            print(e)
            pass

    def init_skyfield(self):
        """ download ephemeris data or read them from file """
        global ts, eph
        # instanciate the loader
        load = Loader(self.path,verbose=False)
        # load timescale
        ts = load.timescale(builtin=self.builtin)
        # load ephemeris
        eph = load(self.eph_file)


class SkyfieldService(StdService):
    """ Service to initialize the Skyfield almanac extension """

    def __init__(self, engine, config_dict):
        global almanacs
        self.path = config_dict.get('DatabaseTypes',configobj.ConfigObj()).get('SQLite',configobj.ConfigObj()).get('SQLITE_ROOT','.')
        alm_conf_dict = config_dict.get('Almanac',configobj.ConfigObj())
        # thread to initialize Skyfield
        self.skyfield_thread = SkyfieldAlmanacThread(alm_conf_dict,self.path)
        self.skyfield_thread.start()
        # instantiate the Skyfield almanac
        self.skyfield_almanac = SkyfieldAlmanacType()
        # add to the list of almanacs
        almanacs.insert(0,self.skyfield_almanac)
    
    def shutDown(self):
        global almanacs
        # find the Skyfield almanac in the list of almanac
        idx = almanacs.index(self.skyfield_almanac)
        # remove it from the list
        del almanacs[idx]
        # stop thread
        if self.skyfield_thread.is_alive():
            self.skyfield_thread.shutDown()
