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
        # URL(s) of the timescale file (optional)
        timescale_url = '...'
        # whether to log FTP responses (optional)
        log_ftp = false
        # update interval 1 year (set to 0 for no update)
        update_interval = 31557600
    
    Skyfield downloads the timescale file finals2000A.all from a server that
    is temporarily down. If you set `use_builtin_timescale` to `false` and
    get permanent download errors, try another origin by setting
    `timescale_url` to an appropriate URL.
    
    Attributes of `almanac_obj`:
    
    lat(float):         latitude in degrees
    lon(float):         longitude in degrees
    altitude(float):    altitude of the location in **meters**
    temperature(float): temperature in **degrees Celsius**
    pressure(float):    pressure in **mbar**
    horizon(float):     horizon in degrees
    
"""

VERSION = "0.1"

# IERS timescale file as hardcoded in Skyfield
TIMESCALE_FILE = 'finals2000A.all'

import time
import datetime
import configobj
import threading
import os
import os.path
import requests
from ftplib import FTP, FTP_TLS

import weewx
from weewx.almanac import AlmanacType, almanacs, timestamp_to_djd
from weewx.engine import StdService
from weewx.units import ValueTuple, ValueHelper, std_groups, Formatter
import weeutil

# Import Skyfield modules
import numpy
from skyfield import almanac
from skyfield.api import N, S, E, W, Loader, wgs84
from skyfield.earthlib import refraction
from skyfield.searchlib import find_discrete, find_maxima
from skyfield.data import iers
from skyfield.constants import DAY_S

# Global variables
ts = None
eph = None

# Unit group and unit used for true solar time and local mean time
for _, unitgroup in weewx.units.std_groups.items():
    unitgroup['group_localtime'] = 'local_djd'

# Logging
import weeutil.logger
import logging
log = logging.getLogger("user.skyfieldalmanac")

def logdbg(msg):
    log.debug(msg)

def loginf(msg):
    log.info(msg)

def logerr(msg):
    log.error(msg)


def timestamp_to_skyfield_time(timestamp):
    """ convert Unix timestamp to Skyfield Time
    
        See https://github.com/skyfielders/python-skyfield/discussions/1027
        for why the timestamp is divided and added to the day of month
        instead of using the `second` parameter of `ts.utc`.
        
        Args:
            timestamp(int, float): Unix timestamp
        
        Returns:
            skyfield.units.Time: the same timestamp converted
    """
    return ts.utc(1970,1,1+timestamp/DAY_S)

def skyfield_time_to_djd(ti):
    """ convert Skyfield timestamp to Dublin Julian Date
    
        Args:
            ti(skyfield.units.Time): timestamp to convert
        
        Returns:
            float: the same timestamp as Dublin Julian Date
    """
    return ti.ut1-ti.dut1/DAY_S-2415020.0

def _get_observer(almanac_obj, target, use_center):
    """ get observer object and refraction angle """
    # a location on earth surface
    observer = eph['Earth'] + wgs84.latlon(almanac_obj.lat,almanac_obj.lon,elevation_m=almanac_obj.altitude)
    # calculate refraction angle
    if almanac_obj.pressure or almanac_obj.horizon:
        horizon = almanac_obj.horizon
        # Using `refraction()` switches off Skyfield's own target radius
        # calculation. So you have to take into account the body's radius 
        # at your own.
        # https://github.com/skyfielders/python-skyfield/discussions/1042
        if not use_center and horizon>-1.0:
            if target.lower()=='sun': horizon -= 16.0/60.0
        # `refraction()` returns the influence of refraction only. It does
        # not include the horizon into the result.
        refr = refraction(
            horizon,
            temperature_C=almanac_obj.temperature,
            pressure_mbar=almanac_obj.pressure
        )
        horizon -= refr
    else:
        horizon = None
    return observer, horizon, eph[target.replace('_',' ')]


class SkyfieldAlmanacType(AlmanacType):
    """ Almanac extension to use the Skyfield module for almanac computation """
    
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
        """ Skyfield provides extras. 
        
            Depending on the ephemeris file chosen, Skyfield takes some time
            to initialize after the start of WeeWX. Initialization is 
            finished when `eph` is not `None` any more.
        
        """
        return eph is not None

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
        elif attr=='solar_time' or attr=='solar_angle':
            val = almanac_obj.sun.ha+180.0
            if attr == 'solar_time':
                return val
            else:
                vt = ValueTuple(val, 'degree_compass', 'group_direction')
                return weewx.units.ValueHelper(vt,
                                               context = 'ephem_day',
                                               formatter=almanac_obj.formatter,
                                               converter=almanac_obj.converter)
        elif attr=='solar_datetime':
            # True Apparent Solar Time (sundial time)
            # We use NumPy because the hour angle is already
            # calculated using NumPy and a NumPy data type.
            # convert hour angle to part of day
            ha = almanac_obj.sun.ha/360.0
            # local Julian Date
            ti = time_ti.ut1+almanac_obj.lon/360.0
            # find solar Julian Day from mean Julian Date
            if ha>0.5:
                ti = numpy.floor(ti)
            elif ha<-0.5:
                ti = numpy.ceil(ti)
            else:
                ti = numpy.round(ti,0)
            # Julian Day (not Julian Date) + hour angle, converted
            # to local solar Dublin Julian Date
            vt = ValueTuple(ti+ha-2415020.0,'local_djd','group_localtime')
            formatter = SkyfieldFormatter(
                unit_label_dict=almanac_obj.formatter.unit_label_dict,
                time_format_dict=almanac_obj.formatter.time_format_dict,
                ordinate_names=almanac_obj.formatter.ordinate_names,
                deltatime_format_dict=almanac_obj.formatter.deltatime_format_dict
            )
            return ValueHelper(vt,
                                context="day",
                                formatter=formatter,
                                converter=almanac_obj.converter)
        elif attr in eph:
            # The attribute is a heavenly body (such as 'sun', or 'jupiter').
            # Bind the almanac and the heavenly body together and return as an
            # SkyfieldAlmanacBinder
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
        observer, horizon, body = _get_observer(self.almanac,self.heavenly_body,self.use_center)
        timespan = weeutil.weeutil.archiveDaySpan(self.almanac.time_ts)
        t0 = timestamp_to_skyfield_time(timespan[0])
        t1 = timestamp_to_skyfield_time(timespan[1])
        tr, yr = almanac.find_risings(observer, body, t0, t1, horizon_degrees=horizon)
        ts, ys = almanac.find_settings(observer, body, t0, t1, horizon_degrees=horizon)
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
        
        if attr in ('astro_ra','astro_dec','astro_dist','a_ra','a_dec','a_dist',
                    'geo_ra','geo_dec','geo_dist','g_ra','g_dec','g_dist'):
            t = timestamp_to_skyfield_time(self.almanac.time_ts)
            body = eph[self.heavenly_body.replace('_',' ')]
            astrometric = eph['Earth'].at(t).observe(body)
            if attr in ('geo_ra','geo_dec','geo_dist','g_ra','g_dec','g_dist'):
                astrometric = astrometric.apparent()
            ra, dec, distance = astrometric.radec(epoch='date')
            if attr in ('a_ra','g_ra'):
                return ra._degrees
            elif attr in ('a_dec','g_dec'):
                return dec.degrees
            elif attr in ('a_dist','g_dist'):
                return distance.km
            if attr in ('astro_ra','geo_ra'):
                vt = ValueTuple(ra._degrees,'degree_compass','group_direction')
            elif attr in ('astro_dec','geo_dec'):
                vt = ValueTuple(dec.radians,'radian','group_angle')
            else:
                vt = ValueTuple(distance.km,'km','group_distance')
            return ValueHelper(vt,
                                               context="ephem_day",
                                               formatter=self.almanac.formatter,
                                               converter=self.almanac.converter)
        
        observer, horizon, body = _get_observer(self.almanac,self.heavenly_body,self.use_center)

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
        elif attr in ('rise','set','transit','antitransit','max_alt','max_alt_time','max_altitude'):
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
            if attr in ('az','alt','alt_dist','azimuth','altitude','alt_distance'):
                alt, az, distance = position.altaz(temperature_C=self.almanac.temperature,pressure_mbar=self.almanac.pressure)
                if attr=='az':
                    return az.degrees
                elif attr=='alt':
                    return alt.degrees
                elif attr=='alt_dist':
                    return distance.km
                else:
                    if attr=='azimuth':
                        vt = ValueTuple(az.degrees,'degree_compass','group_direction')
                    elif attr=='altitude':
                        vt = ValueTuple(alt.radians,'radian','group_angle')
                    elif attr=='alt_distance':
                        vt = ValueTuple(distance.km,'km','group_distance')
                    return ValueHelper(vt,
                                               context="ephem_day",
                                               formatter=self.almanac.formatter,
                                               converter=self.almanac.converter)
            if attr in ('ra','dec','dist','topo_ra','topo_dec','topo_dist'):
                ra, dec, distance = position.radec('date')
                if attr=='ra':
                    return ra._degrees
                elif attr=='dec':
                    return dec.degrees
                elif attr=='dist':
                    return distance.km
                else:
                    if attr=='topo_ra':
                        vt = ValueTuple(ra._degrees,'degree_compass','group_direction')
                    elif attr=='topo_dec':
                        vt = ValueTuple(dec.radians,'radian','group_angle')
                    elif attr=='topo_dist':
                        vt = ValueTuple(distance.km,'km','group_distance')
                    return ValueHelper(vt,
                                               context="ephem_day",
                                               formatter=self.almanac.formatter,
                                               converter=self.almanac.converter)
            if attr in ('ha','ha_dec','ha_dist',
                        'hour_angle','ha_declination','ha_distance'):
                # measured from the plane of the Earth's physical geographic
                # equator. The coordinates are not adjusted for atmospheric
                # refraction near the horizon.
                # https://rhodesmill.org/skyfield/api-position.html#skyfield.positionlib.ICRF.hadec
                ha, dec, distance = position.hadec()
                if attr=='ha':
                    return ha._degrees
                elif attr=='ha_dec':
                    return dec.degrees
                elif attr=='ha_dist':
                    return distance.km
                else:
                    if attr=='hour_angle':
                        vt = ValueTuple(ha._degrees,'degree_compass','group_direction')
                    elif attr=='ha_declination':
                        vt = ValueTuple(dec.radians,'radian','group_angle')
                    elif attr=='ha_distance':
                        vt = ValueTuple(distance.km,'km','group_distance')
                    else:
                        vt = None
                    return ValueHelper(vt,
                                       context="ephem_day",
                                       formatter=self.almanac.formatter,
                                       converter=self.almanac.converter)
            # `attr` is not provided by this extension. So raise an exception.
            raise AttributeError("%s.%s" % (self.heavenly_body,attr))

        # Note: In case of polar day or night, y is False and t is the time
        #       of transit or antitransit.
        
        t = None
        if evt in ('rise','rising'):
            # rising
            t, y = almanac.find_risings(observer, body, t0, t1, horizon_degrees=horizon)
        elif evt in ('set','setting'):
            # setting
            t, y = almanac.find_settings(observer, body, t0, t1, horizon_degrees=horizon)
        elif evt=='transit':
            # meridian transit
            t = almanac.find_transits(observer, body, t0, t1)
            y = True
        elif evt=='antitransit':
            # antitransit
            def az_degrees(t):
                position = observer.at(t).observe(body).apparent()
                ha, _, _ = position.hadec()
                return ha.radians>=0
            az_degrees.step_days = 0.5
            t, val = find_discrete(t0, t1, az_degrees)
            if len(t)>1:
                if val[0]==1:
                    t = t[1:]
                    val = val[1:]
                else:
                    t = t[0:1]
                    val = val[0:1]
            y = len(t)>=1
        elif evt in ('max_alt','max_alt_time','max_altitude'):
            def alt_degrees(t):
                position = observer.at(t).observe(body).apparent()
                alt, _, _ = position.altaz()
                return alt.radians
            alt_degrees.step_days = 0.5
            t, val = find_maxima(t0, t1, alt_degrees)
            if evt=='max_alt' or evt=='max_altitude':
                val = val[-1] if len(val)>=1 else None
                if evt=='max_alt': return val
                return ValueHelper(ValueTuple(val,"radian","group_angle"),
                                   context="day",
                                   formatter=self.almanac.formatter,
                                   converter=self.almanac.converter)
            y = len(t)>=1
        else:
            # `attr` is not provided by this extension. So raise an exception.
            raise AttributeError("%s.%s" % (self.heavenly_body,attr))
        if t is not None:
            time_djd = skyfield_time_to_djd(t[-1]) if len(t)>=1 and y else None
            return weewx.units.ValueHelper(ValueTuple(time_djd, "dublin_jd", "group_time"),
                                           context="ephem_day",
                                           formatter=self.almanac.formatter,
                                           converter=self.almanac.converter)


class SkyfieldFormatter(weewx.units.Formatter):
    """ special formatter including solar time """

    def _to_string(self, val_t, context='current', addLabel=True,
                   useThisFormat=None, None_string=None,
                   localize=True):
        if val_t is not None and val_t[0] is not None:
            if val_t[1]=='local_djd':
                ti = time.gmtime((val_t[0]-25567.5)*DAY_S)
                if useThisFormat is None:
                    val_str = time.strftime(self.time_format_dict.get(context, "%d-%b-%Y %H:%M"),
                                            ti)
                else:
                    val_str = time.strftime(useThisFormat, ti)
                return val_str
        return super(SkyfieldFormatter,self)._to_string(val_t,
                                                context=context,
                                                addLabel=addLabel,
                                                useThisFormat=useThisFormat,
                                                None_string=None_string,
                                                localize=localize)


class SkyfieldMaintenanceThread(threading.Thread):
    """ Thread to download and update ephemeris and timescales 
    
        If initializing Skyfield requires downloading files, this can take
        too long to do it during WeeWX initialization. So it is put into
        a separate thread.
        
        Ephemeris and timescale files are updated from time to time. This
        thread re-downloads them at a given interval if configured to do so.
    """
    
    def __init__(self, alm_conf_dict, path):
        """ init thread
        
            Args:
                alm_conf_dict(configobj): almanac configuration
                path(str): directory where to save downloaded files
        """
        super(SkyfieldMaintenanceThread,self).__init__(name='SkyfieldMaintenanceThread')
        self.path = path
        logdbg("path to save Skyfield files: '%s'" % self.path)
        self.eph_file = alm_conf_dict.get('ephemeris','de440s.bsp')
        self.builtin = weeutil.weeutil.to_bool(alm_conf_dict.get('use_builtin_timescale',True))
        self.update_interval = weeutil.weeutil.to_int(alm_conf_dict.get('update_interval',31557600))
        if self.update_interval:
            self.update_interval = max(self.update_interval,86400)
        self.log_ftp = weeutil.weeutil.to_bool(alm_conf_dict.get('log_ftp',False))
        self.ts_urls = alm_conf_dict.get('timescale_url',None)
        if self.ts_urls and not isinstance(self.ts_urls,list):
            self.ts_urls = [self.ts_urls]
        loginf("ephemeris file: '%s', timescale: %s, update interval: %.2f days" % (self.eph_file,'builtin' if self.builtin else 'IERS file',self.update_interval/DAY_S))
        self.evt = threading.Event()
        self.running = True
        self.last_ts_update = 0
        self.last_eph_update = 0
        logdbg("thread '%s': initialized" % self.name)
    
    def shutDown(self):
        """ shut down thread """
        self.running = False
        self.evt.set()
        loginf("thread '%s': shutdown requested" % self.name)
    
    def run(self):
        """ Skyfield database maintenance """
        loginf("thread '%s': starting" % (self.name,))
        try:
            while self.running:
                # initialize Skyfield or update its database
                success = self.init_skyfield()
                logdbg("thread '%s': Initialization/update was%s successful." % (self.name,'' if success else ' not'))
                # If no updating is required, the thread can be closed now.
                if not self.update_interval: break
                # Wait for the update interval to pass.
                self.evt.wait(self.update_interval if success else 300)
        except Exception as e:
            logerr("thread '%s': %s - %s" % (self.name,e.__class__.__name__,e))
        finally:
            loginf("thread '%s': stopped" % self.name)

    def init_skyfield(self):
        """ download ephemeris data or read them from file """
        global ts, eph
        # instantiate the loader
        load = Loader(self.path,verbose=False)
        # get current time
        now = time.time()-86400
        # load timescale
        if self.last_ts_update<=now or ts is None:
            try:
                if not self.builtin and self.ts_urls:
                    # download timescale from a different location
                    file = self.download(self.ts_urls,'timescale.tmp')
                    if file:
                        os.rename(file,os.path.join(self.path,TIMESCALE_FILE))
                _ts = load.timescale(builtin=self.builtin)
                if _ts: 
                    ts = _ts
                    self.last_ts_update = time.time()
                    loginf("thread '%s': timescale initialized or updated" % self.name)
                if not self.builtin:
                    url = load.build_url(TIMESCALE_FILE)
                    with load.open(url) as f:
                        finals_data = iers.parse_x_y_dut1_from_finals_all(f)
                    iers.install_polar_motion_table(ts, finals_data)
                    loginf("thread '%s': installed polar motion table" % self.name)
            except OSError as e:
                logerr("thread '%s': error downloading timescale %s - %s" % (self.name,e.__class__.__name__,e))
        # load ephemeris
        if self.last_eph_update<=now or eph is None:
            try:
                _eph = load(self.eph_file)
                if _eph: 
                    eph = _eph
                    self.last_eph_update = time.time()
                    loginf("thread '%s': ephemeris initialized or updated" % self.name)
            except OSError as e:
                logerr("thread '%s': error downloading ephemeris %s - %s" % (self.name,e.__class__.__name__,e))
        # `eph` and `ts` are up to date if they were updated less than 24 
        # hours ago.
        return self.last_ts_update>now and self.last_eph_update>now
    
    def download(self, urls, filename=None):
        # try URLs in order
        for url in urls:
            x = url.split(':')[0].lower()
            if x in ('ftp','ftps'):
                fn = self.ftp_download(url, filename)
            elif x in ('http','https'):
                fn = self.http_download(url, filename)
            else:
                fn = None
            if fn is not None: return fn
        return None
    
    def http_download(self, url, filename=None):
        try:
            # target path
            if filename is None:
                x = url.split('/')
                if len(x)>=4: filename = x[-1]
            filename = os.path.join(self.path,filename)
            # download
            headers = {'User-Agent':'weewx-skyfieldalmanac'}
            reply = requests.get(url, headers=headers, timeout=5)
            if reply.status_code==200:
                with open(filename,'wb') as f:
                    f.write(reply.content)
                if self.log_ftp:
                    loginf("thread '%s': successfully downloaded '%s' to '%s'" % (self.name,url,filename))
                return filename
            logerr("thread '%s': HTTP download failed with error code %s" % (self.name,reply.status_code))
        except (OSError,requests.exceptions.Timeout) as e:
            logerr("thread '%s': HTTP download %s - %s" % (self.name,e.__class__.__name__,e))
        return None
    
    def ftp_download(self, url, filename=None):
        """ Download from FTP server """
        try:
            # split URL
            x = url.split(':')
            protocol = x[0]
            x = x[1].split('/')
            if protocol.lower() not in ('ftp','ftps'): raise OSError('unvalid protocol')
            if len(x)<4 or x[0]!='' or x[1]!='': raise OSError('invalid URL')
            server = x[2]
            path = '/'.join(x[3:])
            logdbg("thread '%s': download protocol='%s', server='%s', path='%s'" % (self.name,protocol,server,path))
            # target path
            if filename is None: filename = x[-1]
            filename = os.path.join(self.path, filename)
            # open connection
            if protocol=='ftps':
                ftp = FTP_TLS(server)
            else:
                ftp = FTP(server)
            # login
            x = ftp.login()
            if self.log_ftp: loginf(x)
            # switch to private connection in case of FTPS
            if protocol=='ftps':
                x = ftp.prot_p()
                if self.log_ftp: loginf(x)
            # download
            with open(filename,'wb') as f:
                ftp.retrbinary('RETR %s' % path,f.write)
            # close connection
            x = ftp.quit()
            if self.log_ftp: loginf(x)
            if self.log_ftp:
                loginf("thread '%s': successfully downloaded '%s' to '%s'" % (self.name,url,filename))
            return filename
        except OSError as e:
            logerr("thread '%s': FTP download %s - %s" % (self.name,e.__class__.__name__,e))
            return None
        

class SkyfieldService(StdService):
    """ Service to initialize the Skyfield almanac extension """

    def __init__(self, engine, config_dict):
        """ init this extension """
        global almanacs
        super(SkyfieldService,self).__init__(engine, config_dict)
        # directory to save ephemeris and IERS files
        self.path = config_dict.get('DatabaseTypes',configobj.ConfigObj()).get('SQLite',configobj.ConfigObj()).get('SQLITE_ROOT','.')
        # configuration
        alm_conf_dict = config_dict.get('Almanac',configobj.ConfigObj())
        # thread to initialize Skyfield
        self.skyfield_thread = SkyfieldMaintenanceThread(alm_conf_dict,self.path)
        self.skyfield_thread.start()
        # instantiate the Skyfield almanac
        self.skyfield_almanac = SkyfieldAlmanacType()
        # add to the list of almanacs
        almanacs.insert(0,self.skyfield_almanac)
        logdbg("%s started" % self.__class__.__name__)
    
    def shutDown(self):
        """ remove this extension from the almanacs list and shut down the
            maintenance thread
        """
        global almanacs
        # find the Skyfield almanac in the list of almanacs
        idx = almanacs.index(self.skyfield_almanac)
        # remove it from the list
        del almanacs[idx]
        # stop thread
        if self.skyfield_thread.is_alive():
            self.skyfield_thread.shutDown()
        logdbg("%s stopped" % self.__class__.__name__)


# log version info at startup
loginf("%s version %s" % ("Skyfield almanac extension",VERSION))
