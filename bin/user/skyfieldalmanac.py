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
        [[Skyfield]]
            # Use this almanac
            enable = true
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
    moon_phases(list):  moon phases in local language
    venus_phases(list): venus phases in local language
    mercury_phases(list): mercury phases in local language
    planet_names(list): names of the planets in local language
    formatter (Formatter): 
    converter (Converter):
    
"""

VERSION = "0.4"

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
import json

import weewx
from weewx.almanac import AlmanacType, almanacs, timestamp_to_djd, PyEphemAlmanacType
from weewx.engine import StdService
from weewx.units import ValueTuple, ValueHelper, std_groups, Formatter
import weeutil.weeutil
import weewx.defaults

# Import Skyfield modules
import numpy
from skyfield import VERSION as SKYFIELD_VERSION
from skyfield import almanac
try:
    from skyfield.api import N, S, E, W, Loader, wgs84, EarthSatellite, Star, Angle
except (ImportError,PermissionError):
    N = E = +1.0
    S = W = -1.0
    from skyfield.iokit import Loader
    from skyfield.toposlib import wgs84
    from skyfield.sgp4lib import EarthSatellite
    from skyfield.starlib import Star
    from skyfield.units import Angle
from skyfield.earthlib import refraction
from skyfield.searchlib import find_discrete, find_maxima, find_minima
from skyfield.data import iers
from skyfield.constants import DAY_S, DEG2RAD, RAD2DEG, pi, tau
from skyfield.iokit import parse_tle_file
from skyfield.magnitudelib import planetary_magnitude
from skyfield.named_stars import named_star_dict
from skyfield.framelib import ecliptic_frame
from skyfield.nutationlib import iau2000b_radians
from skyfield.positionlib import position_of_radec
from skyfield.trigonometry import position_angle_of

try:
    import pandas
    has_pandas = True
except ImportError:
    has_pandas = False

# Constants
# Note: Despite Pluto misses the current definition of a planet, it is
#       included in JPL's ephemeris files of planets. Therefore it is
#       included here, too.
PLANETS = ('mercury','venus','earth','mars','jupiter','saturn','uranus','neptune','pluto')
PLANETS_IDX = {j:i for i,j in enumerate(PLANETS)}
PLANETS_IDX.update({'%s_barycenter' % i:j for i,j in PLANETS_IDX.items() if j>2})
SUN = 'sun'
EARTH = PLANETS[2]
EARTHMOON = 'moon'

# Size of heavenly bodies
# Note: Unfortunately the SPICE kernel files do not contain information 
#       about the size of the heavenly bodies, but PyEphem does.
SIZES = {
    #           equator.r  pole.r    mean.r [all km]
    # https://iopscience.iop.org/article/10.1088/0004-637X/750/2/135
    # Note: The source does not differentiate between equatorial, polar, and
    #       mean radius of the sun. Therefore the same value is provided
    #       for all three radii here.
    'sun':     (696342,696342,696342),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
    'mercury': ( 2440.5,   2438.3,   2439.7  ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
    'venus':   ( 6051.8,   6051.8,   6051.8  ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
    'earth':   ( 6378.137, 6356.752, 6371.000),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
    'moon':    ( 1738.1,   1736.0,   1737.4  ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
    'mars':    ( 3396.2,   3376.2,   3389.5  ),
    # https://www.nature.com/articles/nature03938
    'ceres':   (  487.3,    454.7,    476,2  ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
    'jupiter': (71492,    66854,    69911    ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
    'saturn':  (60268,    54364,    58232    ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranusfact.html
    'uranus':  (25559,    24973,    25362    ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html
    'neptune': (24764,    24341,    24622    ),
    # https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
    'pluto':   ( 1188,     1188,     1188    ),
}
MEAN_MOON_RADIUS_KM = SIZES[EARTHMOON][2]
SUN_RADIUS_KM = SIZES[SUN][2]

DEFAULT_PHASES = weewx.defaults.defaults['Almanac']['moon_phases']

# Global variables
ts = None
ephemerides = None
sun_and_planets = None
stars = None
starnames = dict() # HIP number to name
starids = dict()   # name to HIP number
constellation_at = None # constellation function
constellation_names = None # dictionary of abbrevations to names
planets_list = []  # list of planets with available ephemeris
satcatalogues = set()

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


def _get_config(config_dict):
    """ get almanac configuration """
    conf_dict = config_dict.get('Almanac',configobj.ConfigObj()).get('Skyfield',configobj.ConfigObj())
    alm_conf_dict = weeutil.config.accumulateLeaves(conf_dict)
    alm_conf_dict['enable'] = weeutil.weeutil.to_bool(conf_dict.get('enable',True))
    alm_conf_dict['log_success'] = weeutil.weeutil.to_bool(alm_conf_dict.get('log_success',True))
    alm_conf_dict['log_failure'] = weeutil.weeutil.to_bool(alm_conf_dict.get('log_failure',True))
    alm_conf_dict['EarthSatellites'] = conf_dict.get('EarthSatellites',configobj.ConfigObj())
    return alm_conf_dict

def timestamp_to_skyfield_time(timestamp, offset=0):
    """ convert Unix timestamp to Skyfield Time
    
        See https://github.com/skyfielders/python-skyfield/discussions/1027
        for why the timestamp is divided and added to the day of month
        instead of using the `second` parameter of `ts.utc`.
        
        `None` is not a valid timestamp, but it is not easy to find the 
        reason if the exception is raised here. So return a valid `Time`
        object with an invalid value instead of raising TypeError.
        
        Args:
            timestamp(int, float): Unix timestamp
            offset(int, float): offset in seconds to add
        
        Returns:
            skyfield.units.Time: the same timestamp converted
    """
    if timestamp is None: return ts.utc(None)
    return ts.utc(1970,1,1+(timestamp+offset)/DAY_S)

def skyfield_time_to_djd(ti):
    """ convert Skyfield timestamp to Dublin Julian Date
    
        Args:
            ti(skyfield.units.Time): timestamp to convert
        
        Returns:
            float: the same timestamp as Dublin Julian Date
    """
    return ti.ut1-ti.dut1/DAY_S-2415020.0

def skyfield_time_to_timestamp(ti):
    return (ti.ut1-ti.dut1/DAY_S-2440587.5)*DAY_S

def hip_to_starname(hip, default=None):
    """ convert Hipparcos catalogue number to the name of the star """
    global starnames
    return starnames.get(hip,default)

def _get_body(body):
    """ get the ephemeris of the body """
    # Note: sun_and_planets['jupiter barycenter'] and sun_and_planets['jupiter_barycenter'] both work.
    global ephemerides, stars
    try:
        if body.startswith('HIP'):
            x = weeutil.weeutil.to_int(body[3:])
            return Star.from_dataframe(stars.loc[x])
        return ephemerides[body]
    except KeyError:
        raise AttributeError("unknown heavenly body '%s'" % body)

def _get_observer(almanac_obj, target, use_center, with_refraction=True):
    """ get observer object and refraction angle """
    global ephemerides
    # a location on earth surface
    observer = ephemerides[EARTH] + wgs84.latlon(almanac_obj.lat,almanac_obj.lon,elevation_m=almanac_obj.altitude)
    # calculate refraction angle
    if almanac_obj.pressure or almanac_obj.horizon:
        horizon = almanac_obj.horizon
        # Using `refraction()` switches off Skyfield's own target radius
        # calculation. So you have to take into account the body's radius 
        # at your own.
        # https://github.com/skyfielders/python-skyfield/discussions/1042
        if not use_center and horizon>-1.0:
            if target.lower()==SUN: horizon -= 16.0/60.0
        # `refraction()` returns the influence of refraction only. It does
        # not include the horizon into the result.
        # Zero refraction is returned both for objects very near the zenith,
        # as well as for objects more than one degree below the horizon.
        if with_refraction:
            refr = refraction(
                horizon,
                temperature_C=almanac_obj.temperature,
                pressure_mbar=almanac_obj.pressure
            )
            horizon -= refr
    else:
        horizon = None
    return observer, horizon, _get_body(target)

def _get_phases(almanac_obj, body):
    """ list of phases """
    key = '%s_phases' % body
    return almanac_obj.__dict__.get(key,almanac_obj.__dict__.get('texts',dict()).get(key,DEFAULT_PHASES))

def _pyephem_elongation(position):
    """ calculate elongation
    
        PyEphem states that the elongation is set to negative
        when the body is at the morning side of the sky.
    
        There are 2 different definitions of elongation. One is the
        the apparent angle between the body and the Sun. The other 
        one is the angle projected to the ecliptic. I do not know
        for sure which one of them PyEphem uses, but here the real 
        apparent angle is used.
    """
    # Get elongation
    earth = sun_and_planets[position.center].at(position.t)
    sun = earth.observe(ephemerides[SUN]).apparent()
    elong = position.separation_from(sun)
    # Check morning or evening side
    # Note: The extrema of the elongation are not exactly at the difference
    #       of the ecliptic longitude of 0.0° or 180.0°. Thus, if we are near
    #       the extremum, we have to check whether we are before or after it.
    _, lon_earth, _ = earth.frame_latlon(ecliptic_frame)
    _, lon_body, _ = sun_and_planets[position.target].at(position.t).frame_latlon(ecliptic_frame)
    dir = lon_body.degrees-lon_earth.degrees
    if dir<0.0: dir += 360.0
    logdbg('elong %s dir orig %s' % (elong,dir))
    if dir<2.0 or dir>358.0:
        # near lower conjunction or opposition
        def func(t):
            p = sun_and_planets[position.center].at(t)
            return p.observe(sun_and_planets[position.target]).apparent().separation_from(p.observe(ephemerides[SUN]).apparent()).degrees
        func.step_days = 30
        if elong.degrees<90.0:
            # lower conjunction 
            t, v = find_minima(position.t-20, position.t+20, func)
        else:
            # opposition
            t, v = find_maxima(position.t-20, position.t+20, func)
        logdbg('elong position.t %s extremum.t %s' % (position.t,t))
        if len(t)==1:
            # check whether the actual timestamp is before or after the
            # extremum
            dir = 359.0 if position.t<t[0] else 1.0
            logdbg('elong dir %s' % dir)
    elif 178.0<dir<182.0:
        # near upper conjunction
        def func(t):
            p = sun_and_planets[position.center].at(t)
            return p.observe(sun_and_planets[position.target]).apparent().separation_from(p.observe(ephemerides[SUN]).apparent()).degrees
        func.step_days = 30
        t, v = find_minima(position.t-20, position.t+20, func)
        loginf('elong position.t %s t %s' % (position.t,t))
        if len(t)==1:
            # check whether the actual timestamp is before or after the
            # extremum
            dir = 179.0 if position.t<t[0] else 181.0
            loginf('elong dir %s' % dir)
    if dir<180.0:
        # morning side
        return Angle(radians=-elong.radians)
    else:
        # evening side
        return elong

def planet_phase(planet, t):
    """ phase of a planet at time t 
    
        Returns:
            phase_angle: phase angle, 0=full 180°=new
            dir: difference of ecliptic longitude
            idx: phase index 0=new ... 4=full ...
    """
    # position of Earth at time t
    e0 = ephemerides[EARTH].at(t)
    # time the light needs to get from the planets to Earth
    diff_t = e0.observe(planet).apparent().light_time
    # timestamp when the light started at the planet to arrive at Earth at time t
    try:
        len(t)
        tv = ts.tt_jd([x.tt for x in (t-diff_t)])
    except TypeError:
        tv = t - diff_t
    # position of the planet at that time
    v = planet.at(tv)
    # direction of the Sun seen from the planet at that time
    s = v.observe(ephemerides[SUN]).apparent()
    # direction of the light starting at the planet to arrive at the position
    # the Earth will be later on at time t
    e = e0-v
    # get the angle between the direction to the Sun and the direction
    # to where the Earth will be
    phase_angle = s.separation_from(e)
    # difference of ecliptic longitude
    _, e0_lon, e0_dist = e0.frame_latlon(ecliptic_frame)
    _, v_lon, v_dist = v.frame_latlon(ecliptic_frame)
    dir = v_lon.degrees-e0_lon.degrees
    dir = numpy.where(dir>=0,dir,dir+360.0)
    # phase index
    if numpy.any(v_dist.au<e0_dist.au):
        # inner planets
        idx = numpy.where(dir<180.0,180.0-phase_angle.degrees,180.0+phase_angle.degrees)
    else:
        # outer planets
        idx = numpy.fill_like(dir,180.0)
    idx = (idx*8.0/360.0+0.5).astype(numpy.uint8)&7
    # return the angle between the direction to the Sun and the direction
    # to where the Earth will be
    return phase_angle, dir, idx

def moon_tilt(moon_alt_radians, sun_alt_radians, az_radians):
    """ Calculate crescent moon tilt
    
        formula according to Karlheinz Schott 2007
        https://falsche-mondneigung.jimdofree.com/b-geometrische-darstellung-und-berechnung/
        and Andrea K. Myers-Beaghton et al. 30.06.2014
        https://www.seas.upenn.edu/~amyers/MoonPaperOnline.pdf
        
        Note that there is a singularity in this formula at new and full moon.
    
        Args:
            moon_alt_radians(float): current moon altitude
            sun_alt_radians(float): current sun altitude
            az_radians(float): azimuth difference moon - sun 
        
        Returns:
            float: crescent moon tilt in radians, 0 is ☽, π is ☾
    """
    try:
        if az_radians<-numpy.pi: az_radians += numpy.pi*2.0
        if az_radians>numpy.pi: az_radians -= numpy.pi*2.0
        a = numpy.arctan((numpy.cos(az_radians)*numpy.sin(moon_alt_radians)-numpy.tan(sun_alt_radians)*numpy.cos(moon_alt_radians))/numpy.sin(az_radians))
        if az_radians>=0.0: a += numpy.pi
        return a
    except ArithmeticError:
        return 0.5*numpy.pi*(1 if az_radians>0 else -1)

def _database_refraction(archive, ti, alt_degrees):
    """ lookup temperature and pressure for the timestamps ti """
    temperature = []
    pressure = []
    for event in skyfield_time_to_timestamp(ti):
        rec = archive.getRecord(event, max_delta=3600)
        if rec is not None:
            if 'outTemp' in rec:
                temp = weewx.units.convert(weewx.units.as_value_tuple(rec, 'outTemp'), "degree_C")[0]
                if temp is None: temp = 15.0
            else:
                temp = 15.0
            if 'barometer' in rec:
                press = weewx.units.convert(weewx.units.as_value_tuple(rec, 'barometer'), "mbar")[0]
                if press is None: press = 1013.25
            else:
                press = 1013.25
        else:
            temp = 15.0
            press = 1013.25
        temperature.append(temp)
        pressure.append(press)
    temperature = numpy.array(temperature)
    pressure = numpy.array(pressure)
    return refraction(alt_degrees,temperature,pressure)

def _fastify(t):
    """ Copyright Brandon Rhodes """
    t._nutation_angles_radians = iau2000b_radians(t)

def _setting_hour_angle(latitude, declination, altitude_radians):
    """Return the hour angle, in radians, when a body reaches the horizon.
       Copyright Brandon Rhodes
    """
    lat = latitude.radians
    dec = declination.radians
    numerator = numpy.sin(altitude_radians) - numpy.sin(lat) * numpy.sin(dec)
    denominator = numpy.cos(lat) * numpy.cos(dec)
    ha = numpy.arccos(numpy.clip(numerator / denominator, -1.0, 1.0))
    return ha

def _rising_hour_angle(latitude, declination, altitude_radians):
    """Return the hour angle, in radians, when a body reaches the horizon.
       Copyright Brandon Rhodes
    """
    return - _setting_hour_angle(latitude, declination, altitude_radians)

def _transit_ha(latitude, declination, altitude_radians):
    """Copyright Brandon Rhodes """
    return 0.0

def _q(a, b, c, sign):
    """Copyright Brandon Rhodes """
    discriminant = numpy.maximum(b*b - 4*a*c, 0.0)  # avoid tiny negative results
    return - 2*c / (b + sign * numpy.sqrt(discriminant))

def _intersection(y0, y1, v0, v1):
    """ Copyright Brandon Rhodes """
    # Return x at which a curve reaches y=0, given its position and
    # velocity y0,v0 at x=0 and y1,v1 at x=1.  For details, see
    # `design/intersect_function.py` in the Skyfield repository.
    sign = 1 - 2 * (y0 > y1)
    return _q(y1 - y0 - v0, v0, y0, sign)

_MICROSECOND = 1 / 24.0 / 3600.0 / 1e6
_clip_lower = -1.0
_clip_upper = +2.0

def _adjust_to_refraction(observer, body, t, y, horizon_degrees, func):
    """ adjust to refraction 
        Copyright Brandon Rhodes
        
        This is copied because we already have an array of times and a
        corresponding array of horizons which is not supported by the
        original function.
    """
    horizon_radians = horizon_degrees / 360.0 * tau
    geo = observer.vector_functions[-1]  # should we check observer.center?
    latitude = geo.latitude
    _fastify(t)
    ha, dec, distance = observer.at(t).observe(body).apparent().hadec()
    old_ha_radians = ha.radians
    old_t = t
    def normalize_zero_to_tau(radians):
        return radians % tau
    def normalize_plus_or_minus_pi(radians):
        return (radians + pi) % tau - pi
    normalize = normalize_zero_to_tau
    if func is _setting_hour_angle:
        t += 0.000694444444444
    else:
        t -= 0.000694444444444
    for i in 0, 1, 2:
        _fastify(t)
        apparent = observer.at(t).observe(body).apparent()
        ha, dec, distance = apparent.hadec()
        desired_ha = func(latitude, dec, horizon_radians)
        ha_adjustment = desired_ha - ha.radians
        ha_adjustment = (ha_adjustment + pi) % tau - pi
        if i < 2:
            ha_diff = normalize(ha.radians - old_ha_radians)
            t_diff = t - old_t
            ha_per_day = ha_diff / t_diff
        old_ha_radians = ha.radians
        old_t = t
        timebump = ha_adjustment / ha_per_day
        timebump[timebump == 0.0] = _MICROSECOND   # avoid divide-by-zero
        previous_t = t
        t = ts.tt_jd(t.whole, t.tt_fraction + timebump)
        normalize = normalize_plus_or_minus_pi
    if func is _transit_ha:
        return t
    v = observer.vector_functions[-1]
    altitude0, _, distance0, rate0, _, _ = (
        apparent.frame_latlon_and_rates(v))
    t.M = previous_t.M
    t._nutation_angles_radians = previous_t._nutation_angles_radians
    apparent = observer.at(t).observe(body).apparent()
    altitude1, _, distance1, rate1, _, _ = (
        apparent.frame_latlon_and_rates(v))
    tdiff = t - previous_t
    t_scaled_offset = _intersection(
        altitude0.radians - horizon_radians,
        altitude1.radians - horizon_radians,
        rate0.radians.per_day * tdiff,
        rate1.radians.per_day * tdiff,
    )
    t_scaled_offset = numpy.clip(t_scaled_offset, _clip_lower, _clip_upper)
    w, f = numpy.divmod(t_scaled_offset * tdiff, 1.0)
    t = previous_t.ts.tt_jd(previous_t.whole + w, previous_t.tt_fraction + f)
    is_above_horizon =  (
        (desired_ha % pi != 0.0)
        | ((t_scaled_offset > _clip_lower) & (t_scaled_offset < _clip_upper))
    )
    return t, is_above_horizon


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
        'last_quarter_moon':(3,),
        'first_quarter_venus':(False,),
        'last_quarter_venus':(True,)
    }

    @property
    def hasExtras(self):
        """ Skyfield provides extras. 
        
            Depending on the ephemeris file chosen, Skyfield takes some time
            to initialize after the start of WeeWX. Initialization is 
            finished when `ephemerides` is not `None` any more.
        
        """
        global ephemerides
        return ephemerides is not None

    def get_almanac_data(self, almanac_obj, attr):
        """ calculate attribute """
        global ephemerides, planets_list, satcatalogues
        if ts is None or ephemerides is None:
            raise weewx.UnknownType(attr)
        if attr=='sunrise':
            return almanac_obj.sun.rise
        elif attr=='sunset':
            return almanac_obj.sun.set
        elif attr=='moon_fullness':
            return int(almanac_obj.moon.moon_fullness + 0.5)
        elif attr=='planets':
            return planets_list
        elif attr in {'previous_equinox', 'next_equinox',
                      'previous_solstice', 'next_solstice',
                      'previous_autumnal_equinox', 'next_autumnal_equinox',
                      'previous_vernal_equinox', 'next_vernal_equinox',
                      'previous_winter_solstice', 'next_winter_solstice',
                      'previous_summer_solstice', 'next_summer_solstice',
                      'previous_new_moon', 'next_new_moon',
                      'previous_first_quarter_moon', 'next_first_quarter_moon',
                      'previous_full_moon', 'next_full_moon',
                      'previous_last_quarter_moon', 'next_last_quarter_moon',
                      'previous_first_quarter_venus','next_first_quarter_venus',
                      'previous_last_quarter_venus','next_last_quarter_venus'}:
            previous = attr.startswith('previous')
            if attr.endswith('_moon'):
                # moon phases
                interval = 2592000
                func = almanac.moon_phases
            elif attr.endswith('_venus'):
                # venus phases
                interval = 50457600
                def func(eph):
                    venus = eph['venus']
                    def phase_at(t):
                        return planet_phase(venus,t)[0].degrees>=90.0
                    phase_at.step_days = 30
                    return phase_at
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
            t0 = timestamp_to_skyfield_time(almanac_obj.time_ts,t0)
            t1 = timestamp_to_skyfield_time(almanac_obj.time_ts,t1)
            # find the events
            t, y  = almanac.find_discrete(t0,t1,func(ephemerides))
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
        elif attr in {'previous_aphelion','next_aphelion',
                      'previous_perihelion','next_perihelion',
                      'previous_apogee_moon','next_apogee_moon',
                      'previous_perigee_moon','next_perigee_moon',
                      'previous_new_venus','next_new_venus',
                      'previous_full_venus','next_full_venus'}:
            earth = ephemerides[EARTH]
            if attr.endswith('helion'):
                # aphelion, perihelion
                interval = 31557600
                sun = ephemerides[SUN]
                # Note: In method `skyfield.searchlib._identify_maxima()` there
                #       is a comment reading "Also choose the midpoint between
                #       the edges of a plateau, if both edges are in view."
                #       Therefore we assume we can round the result to avoid
                #       comparing tenth of millimetres which are for sure no
                #       real differences.
                def func(t):
                    return numpy.round(earth.at(t).observe(sun).apparent().distance().km,1)
                func.step_days = 90.0
            elif attr.endswith('gee_moon'):
                # apogee, perigee of the Moon
                interval = 2592000
                moon = ephemerides[EARTHMOON]
                def func(t):
                    return numpy.round(earth.at(t).observe(moon).apparent().distance().km,1)
                func.step_days = 7.0
            elif attr.endswith('_venus'):
                # phases of the Venus
                interval = 50457600
                venus = ephemerides['venus']
                def func(t):
                    return planet_phase(venus,t)[0].degrees
                func.step_days = 30
            if attr.startswith('previous_'):
                # last previous event
                t0 = -interval
                t1 = 0
                idx = -1
            else:
                # next event
                t0 = 0
                t1 = interval
                idx = 0
            # time interval to look for events
            t0 = timestamp_to_skyfield_time(almanac_obj.time_ts,t0)
            t1 = timestamp_to_skyfield_time(almanac_obj.time_ts,t1)
            # find event
            if attr.endswith('aphelion') or 'apogee' in attr or 'new' in attr:
                t, v = find_maxima(t0, t1, func)
            else:
                t, v = find_minima(t0, t1, func)
            """
            try:
                if len(t)>1:
                    for tt,vv in zip(t.ut1,v):
                        loginf("%s %s %s" % (attr,tt,vv))
            except Exception as e:
                logerr("%s %s" % (e.__class__.__name__,e))
            """
            djd = skyfield_time_to_djd(t)
            return weewx.units.ValueHelper(ValueTuple(djd, "dublin_jd", "group_time"),
                                           context="ephem_year",
                                           formatter=almanac_obj.formatter,
                                           converter=almanac_obj.converter)
        time_ti = timestamp_to_skyfield_time(almanac_obj.time_ts)
        if attr in {'moon_phase','moon_index'}:
            position = almanac.moon_phase(ephemerides, time_ti).degrees/360.0
            moon_index = int((position * 8) + 0.5) & 7
            if attr=='moon_index': return moon_index
            return _get_phases(almanac_obj,'moon')[moon_index]
        elif attr in {'venus_phase','venus_index'}:
            _, _, idx = planet_phase(ephemerides['venus'],time_ti)
            if attr=='venus_index': return idx
            return _get_phases(almanac_obj,'venus')[idx]
        elif attr in {'mercury_phase','mercury_index'}:
            _, _, idx = planet_phase(ephemerides['mercury'],time_ti)
            if attr=='mercury_index': return idx
            return _get_phases(almanac_obj,'mercury')[idx]
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
        elif attr.lower() in ephemerides:
            # The attribute is a heavenly body (such as 'sun', or 'venus').
            # Bind the almanac and the heavenly body together and return as a
            # SkyfieldAlmanacBinder
            return SkyfieldAlmanacBinder(almanac_obj, attr.lower())
        elif (attr.lower()+'_barycenter' in ephemerides and 
             attr in {'mars','jupiter','saturn','uranus','neptune','pluto'}):
            # The attribute is a heavenly body (such as 'jupiter'), but its
            # barycentre is available only. So map the name.
            # Bind the almanac and the heavenly body together and return as a
            # SkyfieldAlmanacBinder
            return SkyfieldAlmanacBinder(almanac_obj,
                                                    attr.lower()+'_barycenter')
        elif '_' in attr and attr.lower().split('_')[0] in satcatalogues:
            # The attribute points to an Earth satellite that is not included
            # in `ephemerides`. Bind it to get the right exception.
            return SkyfieldAlmanacBinder(almanac_obj, attr.lower())
        elif attr.startswith('HIP') and attr[3:].isdigit():
            # The attribute is a star. Bind the almanac and the star together
            # and return as a SkyfieldAlmanacBinder.
            return SkyfieldAlmanacBinder(almanac_obj, attr)
        elif attr.capitalize() in starids:
            # PyEphem maintains a list of 115 wellknown stars by name.
            # To replace PyEphem, name resolution for star names is required
            # here, too. Unfortunately the Hipparcos catalog does not contain
            # names at all.
            return SkyfieldAlmanacBinder(almanac_obj,
                                          'HIP%s' % starids[attr.capitalize()])
        # `attr` is not provided by this extension. So raise an exception.
        raise weewx.UnknownType(attr)


class SkyfieldAlmanacBinder:
    """This class binds the observer properties held in Almanac, with the 
    heavenly body to be observed."""

    def __init__(self, almanac, heavenly_body):
        self.almanac = almanac
        self.heavenly_body = heavenly_body
        self.use_center = False

    def __call__(self, use_center=False):
        self.use_center = use_center
        return self

    def __str__(self):
        """ SkyfieldAlmanacBinder cannot be printed itself. It always needs an 
            attribute.
        """
        raise AttributeError(self.heavenly_body)

    @property
    def visible(self):
        """ Calculate how long the body has been visible today
        
            This property returns the total time the body is up during the 
            archive day the almanac time is in. The body may rise in the 
            morning and set in the evening with a period of visibility in 
            between. But the body can also be up at midnight, then set, and 
            later rise again. In that case there are two separate periods 
            of visibility which will be summarized. Finally the body can
            be up or down the whole day in which case 24 hours or 0 hours
            is returned, respectively. In case of leap seconds or daylight
            savings time switches the result can be more than 86400 seconds.
        """
        observer, horizon, body = _get_observer(self.almanac,self.heavenly_body,self.use_center)
        timespan = weeutil.weeutil.archiveDaySpan(self.almanac.time_ts)
        t0 = timestamp_to_skyfield_time(timespan[0])
        t1 = timestamp_to_skyfield_time(timespan[1])
        if SKYFIELD_VERSION<(1,47):
            # outdated function
            tx, yx = almanac.find_discrete(t0, t1, almanac.risings_and_settings(ephemerides, body, wgs84.latlon(self.almanac.lat,self.almanac.lon,elevation_m=self.almanac.altitude)))
            change = [(t,y==1) for t,y in zip(tx.tai,yx) if y in (0,1)]
        else:
            # actual function
            tr, yr = almanac.find_risings(observer, body, t0, t1, horizon_degrees=horizon)
            up = [(t,True) for t,y in zip(tr.tai,yr) if y]
            tg, yg = almanac.find_settings(observer, body, t0, t1, horizon_degrees=horizon)
            down = [(t,False) for t,y in zip(tg.tai,yg) if y]
            change = up+down
        if change:
            # There is at least one change in visibility.
            change.sort(key=lambda x:x[0])
            change = [(t0.tai,not change[0][1])] + change + [(t1.tai,not change[-1][1])]
            visible = 0
            for idx in range(len(change)-1):
                if change[idx][1] and not change[idx+1][1]:
                    visible += (change[idx+1][0]-change[idx][0]) * weewx.units.SECS_PER_DAY
        else:
            # Always up or always down
            alt, _, _ = observer.at(timestamp_to_skyfield_time(0.5*(timespan[0]+timespan[1]))).observe(body).altaz()
            if alt>0:
                visible = (t1.tai-t0.tai) * weewx.units.SECS_PER_DAY
            else:
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
        if today_visible.raw is None or then_visible.raw is None:
            diff = None
        else:
            diff = today_visible.raw - then_visible.raw
        return weewx.units.ValueHelper(ValueTuple(diff, "second", "group_deltatime"),
                                       context="hour",
                                       formatter=self.almanac.formatter,
                                       converter=self.almanac.converter)

    @property
    def name(self):
        """ name of the heavenly object
        
            The `name` attribute is supported by the PyEphem almanac, too, but
            it is not documented in the customization guide of WeeWX.
        """
        planet_names = self.almanac.__dict__.get('planet_names',self.almanac.__dict__.get('texts',dict()).get('planet_names'))
        if planet_names:
            # If there is a list of the planets' names in local language
            # in the language file, try it.
            idx = PLANETS_IDX.get(self.heavenly_body)
            if idx is not None and idx<len(planet_names):
                return planet_names[idx]
        body = _get_body(self.heavenly_body)
        if isinstance(body,Star):
            # If `_get_body()` returned an instance of class `Star`,
            # we know for sure that `self.heavenly_body` starts
            # with `HIP`.
            hip = weeutil.weeutil.to_int(self.heavenly_body[3:])
            return hip_to_starname(hip, self.heavenly_body)
        if isinstance(body,EarthSatellite):
            # There is always a name attached to an earth satellite,
            # and it is language independent.
            return body.name
        # If no other source of the name is available, use the name
        # of the attribute that specifies the heavenly body which
        # is mostly the English name of it.
        return self.heavenly_body.split('_')[0].capitalize()

    @property
    def hip_number(self):
        """ Hipparcos catalogue number in case of stars """
        if (self.heavenly_body.startswith('HIP') and
                                             self.heavenly_body[3:].isdigit()):
            return weeutil.weeutil.to_int(self.heavenly_body[3:])
        else:
            return None

    def __getattr__(self, attr):
        """Get the requested observation, such as when the body will rise."""
        global ephemerides, constellation_names
        # Don't try any attributes that start with a double underscore, or any of these
        # special names: they are used by the Python language:
        if attr.startswith('__') or attr in ['mro', 'im_func', 'func_code']:
            raise AttributeError(attr)
        
        # The `sun_distance` attribute is supported by the PyEphem almanac, 
        # too, but it is not documented in the customization guide of WeeWX.
        if attr=='sun_distance':
            t = timestamp_to_skyfield_time(self.almanac.time_ts)
            body = _get_body(self.heavenly_body)
            astrometric = ephemerides[SUN].at(t).observe(body)
            return astrometric.distance().au
        
        if attr in {'hlat','hlon','hlatitude','hlongitude'}:
            # astrometric heliocentric latitude and longitude
            # https://github.com/skyfielders/python-skyfield/discussions/1029
            t = timestamp_to_skyfield_time(self.almanac.time_ts)
            body = _get_body(self.heavenly_body)
            #hlat, hlon, _ = body.at(t).frame_latlon(ecliptic_frame)
            hlat, hlon, _ = ephemerides[SUN].at(t).observe(body).frame_latlon(ecliptic_frame)
            if attr=='hlat':
                return hlat.degrees
            elif attr=='hlon':
                return hlon.degrees
            elif attr=='hlatitude':
                vt = ValueTuple(hlat.radians,'radian','group_angle')
            else:
                vt = ValueTuple(hlon.radians,'radian','group_angle')
            return ValueHelper(vt,
                               context="ephem_year",
                               formatter=self.almanac.formatter,
                               converter=self.almanac.converter)
        
        if attr in {'astro_ra','astro_dec','astro_dist','a_ra','a_dec','a_dist',
                    'geo_ra','geo_dec','geo_dist','g_ra','g_dec','g_dist',
                    'earth_distance','elong','elongation','mag',
                    'sublat','sublong','sublatitude','sublongitude',
                    'elevation','size','radius','radius_size',
                    'constellation','constellation_abbr'}:
            t = timestamp_to_skyfield_time(self.almanac.time_ts)
            body = _get_body(self.heavenly_body)
            if isinstance(body,EarthSatellite):
                astrometric = body.at(t)
            else:
                astrometric = ephemerides[EARTH].at(t).observe(body)
                if attr in {'geo_ra','geo_dec','geo_dist','g_ra','g_dec','g_dist','elong','elongation'}:
                    astrometric = astrometric.apparent()
            if attr in {'elong','elongation'}:
                try:
                    elong = _pyephem_elongation(astrometric)
                except Exception as e:
                    logerr('elongation %s %s' % (e.__class__.__name__,e))
                    raise
                if attr=='elong':
                    return elong.degrees
                vt = ValueTuple(elong.radians,'radian','group_angle')
            elif attr=='mag':
                if isinstance(body,Star):
                    return stars.loc[int(self.heavenly_body[3:])]['magnitude']
                try:
                    return float(planetary_magnitude(astrometric))
                except (ValueError,TypeError):
                    return None
            elif attr in {'constellation','constellation_abbr'}:
                # Which constellation the position belongs to?
                if constellation_at:
                    abbr = constellation_at(astrometric)
                    if attr=='constellation_abbr': return abbr
                    if abbr in self.almanac.__dict__.get('texts',dict()).get('Constellations',dict()):
                        return self.almanac.__dict__['texts']['Constellations'][abbr]
                    if abbr and constellation_names:
                        return constellation_names[abbr]
                return 'N/A'
            elif attr in {'sublat','sublong','sublatitude','sublongitude','elevation'}:
                # https://rhodesmill.org/skyfield/coordinates.html#geographic-itrs-latitude-and-longitude
                point = wgs84.geographic_position_of(astrometric)
                if attr=='sublat':
                    return point.latitude.degrees
                elif attr=='sublong':
                    return point.longitude.degrees
                elif attr=='elevation':
                    return point.elevation.m
                elif attr=='sublatitude':
                    vt = ValueTuple(point.latitude.radians,'radian','group_angle')
                elif attr=='sublongitude':
                    vt = ValueTuple(point.longitude.degrees,'degree_compass','group_direction')
            else:
                ra, dec, distance = astrometric.radec(epoch='date')
                if attr in {'a_ra','g_ra'}:
                    return ra._degrees
                elif attr in {'a_dec','g_dec'}:
                    return dec.degrees
                elif attr in {'a_dist','g_dist'}:
                    return distance.km
                elif attr=='earth_distance':
                    return distance.au
                elif attr in {'size','radius','radius_size'}:
                    radius = SIZES.get(self.heavenly_body.split('_')[0])
                    if radius is None or radius[0] is None or not distance.km:
                        radius = None
                    else:
                        radius = radius[0]/distance.km
                    if attr=='radius': 
                        return radius
                    if attr=='size': 
                        return radius*2.0*RAD2DEG*3600.0 if radius is not None else None
                    # `radius_size`
                    vt = ValueTuple(radius,'radian','group_angle')
                if attr in {'astro_ra','geo_ra'}:
                    vt = ValueTuple(ra._degrees,'degree_compass','group_direction')
                elif attr in {'astro_dec','geo_dec'}:
                    vt = ValueTuple(dec.radians,'radian','group_angle')
                else:
                    vt = ValueTuple(distance.km,'km','group_distance')
            return ValueHelper(vt,
                               context="ephem_day",
                               formatter=self.almanac.formatter,
                               converter=self.almanac.converter)
        
        observer, horizon, body = _get_observer(
                               self.almanac,self.heavenly_body,self.use_center)
        
        previous = attr.startswith('previous_')
        next = attr.startswith('next_')
        
        if previous:
            # get the last event before the given timestamp
            t0 = timestamp_to_skyfield_time(self.almanac.time_ts,-86400)
            t1 = timestamp_to_skyfield_time(self.almanac.time_ts)
            evt = attr[9:]
            idx = -1
        elif next:
            # get the next event after the given timestamp
            t0 = timestamp_to_skyfield_time(self.almanac.time_ts)
            t1 = timestamp_to_skyfield_time(self.almanac.time_ts,86400)
            evt = attr[5:]
            idx = 0
        elif attr in {'rise','set','transit','antitransit','day_max_alt','day_max_alt_time','day_max_altitude'}:
            # get the event within the day the timestamp is in
            timespan = weeutil.weeutil.archiveDaySpan(self.almanac.time_ts)
            t0 = timestamp_to_skyfield_time(timespan[0])
            t1 = timestamp_to_skyfield_time(timespan[1])
            evt = attr
            idx = 0
        else:
            # convert given timestamp
            ti = timestamp_to_skyfield_time(self.almanac.time_ts)
            if isinstance(body,EarthSatellite):
                station = wgs84.latlon(self.almanac.lat,self.almanac.lon,elevation_m=self.almanac.altitude)
                position = (body-station).at(ti)
            else:
                position = observer.at(ti).observe(body).apparent()
            if attr=='parallactic_angle':
                # https://github.com/skyfielders/python-skyfield/issues/819#issuecomment-2001972677
                station = wgs84.latlon(self.almanac.lat,self.almanac.lon,elevation_m=self.almanac.altitude)
                north_pole = position_of_radec(
                    ra_hours=0.0, dec_degrees=90.0, epoch=ti,
                    center=station, t=ti,
                )
                pa = position_angle_of(
                    position.altaz(), 
                    north_pole.altaz()
                )
                vt = ValueTuple(pa.radians,'radian','group_angle')
                return ValueHelper(vt,
                               context="ephem_day",
                               formatter=self.almanac.formatter,
                               converter=self.almanac.converter)
            if attr==('%s_fullness' % self.heavenly_body.lower()):
                # `moon_fullness`, `venus_fullness`, `mercury_fullness`
                return position.fraction_illuminated(ephemerides[SUN])*100.0
            if attr=='moon_tilt':
                # tilt of the moon crescent
                # 0 = illuminated side to the right on northern hemisphere
                # π = illuminated side to the left on northern hemisphere
                alt_moon, az_moon, _ = position.altaz()
                alt_sun, az_sun, _ = observer.at(ti).observe(ephemerides[SUN]).apparent().altaz()
                a = moon_tilt(alt_moon.radians,alt_sun.radians,az_moon.radians-az_sun.radians)
                vt = ValueTuple(a,'radian','group_angle')
                return ValueHelper(vt,
                                               context="ephem_day",
                                               formatter=self.almanac.formatter,
                                               converter=self.almanac.converter)
            if attr in {'az','alt','alt_dist','azimuth','altitude','alt_distance'}:
                alt, az, distance = position.altaz(temperature_C=self.almanac.temperature,pressure_mbar=self.almanac.pressure if self.almanac.pressure else 'standard')
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
            if attr in {'ra','dec','dist','topo_ra','topo_dec','topo_dist'}:
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
            if attr in {'ha','ha_dec','ha_dist',
                        'hour_angle','ha_declination','ha_distance'}:
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
        
        if SKYFIELD_VERSION<(1,47):
            station = wgs84.latlon(self.almanac.lat,self.almanac.lon,elevation_m=self.almanac.altitude)
        t = None
        if evt in {'rise','rising'}:
            # rising
            try:
                if SKYFIELD_VERSION<(1,47):
                    if self.heavenly_body==SUN and horizon<(-0.8333):
                        f = almanac.dark_twilight_day(ephemerides, station)
                        what = int(4+horizon/6.0)
                    else:
                        f = almanac.risings_and_settings(ephemerides, body, station)
                        what = 1
                    t, y = almanac.find_discrete(t0, t1, f)
                    t = [i for i,j in zip(t,y) if j==what]
                    y = True
                else:
                    t, y = almanac.find_risings(observer, body, t0, t1, horizon_degrees=horizon)
                    if (t is not None and len(t)>=1 and 
                            self.almanac.horizon==0 and 
                            self.heavenly_body.lower()==EARTHMOON and 
                            not self.use_center):
                        _, _, distance = observer.at(t).observe(body).apparent().hadec()
                        #horizon = self.almanac.horizon-almanac._moon_radius_m/distance.m*RAD2DEG
                        horizon = self.almanac.horizon-MEAN_MOON_RADIUS_KM/distance.km*RAD2DEG
                        horizon -= refraction(horizon,self.almanac.temperature,self.almanac.pressure)
                        t, y = almanac.find_risings(observer, body, t0, t1, horizon_degrees=horizon)
                    y = len(y)>=1 and y[idx]
            except ValueError as e:
                logerr("%s.%s: %s - %s" % (self.heavenly_body,attr,e.__class__.__name__,e))
                t = None
        elif evt in {'set','setting'}:
            # setting
            try:
                if SKYFIELD_VERSION<(1,47):
                    f = almanac.risings_and_settings(ephemerides, body, station)
                    t, y = almanac.find_discrete(t0, t1, f)
                    t = [i for i,j in zip(t,y) if j==0]
                    y = True
                else:
                    t, y = almanac.find_settings(observer, body, t0, t1, horizon_degrees=horizon)
                    if (t is not None and len(t)>=1 and 
                            self.almanac.horizon==0 and 
                            self.heavenly_body.lower()==EARTHMOON and 
                            not self.use_center):
                        _, _, distance = observer.at(t).observe(body).apparent().hadec()
                        #horizon = self.almanac.horizon-almanac._moon_radius_m/distance.m*RAD2DEG
                        horizon = self.almanac.horizon-MEAN_MOON_RADIUS_KM/distance.km*RAD2DEG
                        horizon -= refraction(horizon,self.almanac.temperature,self.almanac.pressure)
                        t, y = almanac.find_settings(observer, body, t0, t1, horizon_degrees=horizon)
                    y = len(y)>=1 and y[idx]
            except ValueError as e:
                logerr("%s.%s: %s - %s" % (self.heavenly_body,attr,e.__class__.__name__,e))
                t = None
        elif evt=='transit':
            # meridian transit
            if SKYFIELD_VERSION<(1,47):
                f = almanac.meridian_transits(ephemerides, body, station)
                t, y = almanac.find_discrete(t0, t1, f)
                t = [i for i,j in zip(t,y) if j==1]
                y = len(t)>=1
            else:
                # If no transit is found, the result is a value of type `Time`
                # that conatains an empty list.
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
        elif evt in {'day_max_alt','day_max_alt_time','day_max_altitude'}:
            def alt_degrees(t):
                position = observer.at(t).observe(body).apparent()
                alt, _, _ = position.altaz()
                return alt.radians
            alt_degrees.step_days = 0.5
            t, val = find_maxima(t0, t1, alt_degrees)
            if evt=='day_max_alt' or evt=='day_max_altitude':
                val = val[-1] if len(val)>=1 else None
                if evt=='day_max_alt': return val
                return ValueHelper(ValueTuple(val,"radian","group_angle"),
                                   context="day",
                                   formatter=self.almanac.formatter,
                                   converter=self.almanac.converter)
            y = len(t)>=1
        else:
            # `attr` is not provided by this extension. So raise an exception.
            raise AttributeError("%s.%s" % (self.heavenly_body,attr))
        time_djd = skyfield_time_to_djd(t[idx]) if t is not None and len(t)>=1 and y else None
        return weewx.units.ValueHelper(ValueTuple(time_djd, "dublin_jd", "group_time"),
                                           context="ephem_day",
                                           formatter=self.almanac.formatter,
                                           converter=self.almanac.converter)

    def genVisibleTimespans(self, context=None, timespan=None, archive=None):
        """ generator function returning uptimes of body
        
            Args:
                context(str): name of the timespan to search, optional
                timespan(TimeSpan): timespan to search, optional, if provided
                    `almanac_time` and `context` are ignored
                archive(Manager): database manager to look up temperature and
                    pressure
            
            Returns:
                TimeSpan: timespan the body is up
        """
        # observer
        observer, horizon, body = _get_observer(self.almanac,self.heavenly_body,self.use_center,with_refraction=False)
        # determine the timespan to search
        if timespan is None:
            if context=='month':
                timespan = weeutil.weeutil.archiveMonthSpan(self.almanac.time_ts)
            elif context=='week':
                timespan = weeutil.weeutil.archiveWeekSpan(self.almanac.time_ts)
            else:
                timespan = weeutil.weeutil.archiveYearSpan(self.almanac.time_ts)
        t0 = timestamp_to_skyfield_time(timespan[0])
        t1 = timestamp_to_skyfield_time(timespan[1])
        # get risings and settings during the timespan
        tri, yri = almanac.find_risings(observer, body, t0, t1, horizon_degrees=horizon)
        tse, yse = almanac.find_settings(observer, body, t0, t1, horizon_degrees=horizon)
        # if a database manager is available look up temperature and pressure
        # and adjust for refraction
        if archive and self.almanac.horizon>-1.0:
            if self.use_center:
                body_radius_degrees = 0.0
            else:
                # TODO: other bodies than the Sun
                body_radius_degrees = 16/60
            hori_rise = _database_refraction(archive,tri,-body_radius_degrees)+body_radius_degrees-self.almanac.horizon
            hori_set = _database_refraction(archive,tse,-body_radius_degrees)+body_radius_degrees-self.almanac.horizon
            tri, yri = _adjust_to_refraction(observer, body, tri, yri, -hori_rise, _rising_hour_angle)
            tse, yse = _adjust_to_refraction(observer, body, tse, yse, -hori_set, _setting_hour_angle)
        # convert to `unix_epoch` timestamps
        tri = skyfield_time_to_timestamp(tri)
        tse = skyfield_time_to_timestamp(tse)
        # If the list starts with a timestamp of setting, remove it
        if len(tri)>0 and len(tse)>0 and tse[0]<tri[0]: 
            tse.pop(0)
            yse.pop(0)
        # Return timespans as WeeWX type `TimeSpan`
        for rising, yrising, setting, ysetting in zip(tri,yri,tse,yse):
            if yrising and ysetting:
                yield weeutil.weeutil.TimeSpan(rising,setting)


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
    """ Thread to download and update ephemerides and timescales 
    
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
        self.log_success = alm_conf_dict['log_success']
        self.log_failure = alm_conf_dict['log_failure']
        self.path = path
        logdbg("path to save Skyfield files: '%s'" % self.path)
        ephem_files = alm_conf_dict.get('ephemeris','de440s.bsp')
        if isinstance(ephem_files,list):
            # A list of files is provided. Save it.
            self.eph_files = ephem_files
            self.spk = [None]*len(ephem_files)
        else:
            # One single file is provided. Convert it into a list.
            self.eph_files = [ephem_files]
            self.spk = [None]
        # Use Skyfield's built-in timescale or retrieve it from IERS?
        self.builtin = weeutil.weeutil.to_bool(alm_conf_dict.get('use_builtin_timescale',True))
        # Update interval for timescale (if not built-in) and ephemerides
        self.update_interval = weeutil.weeutil.to_int(alm_conf_dict.get('update_interval',31557600))
        if self.update_interval:
            self.update_interval = max(self.update_interval,86400)
        # For debegging the FTP communication can be logged.
        self.log_ftp = weeutil.weeutil.to_bool(alm_conf_dict.get('log_ftp',False))
        # URLs to try to get timescale files from
        self.ts_urls = alm_conf_dict.get('timescale_url',None)
        if self.ts_urls and not isinstance(self.ts_urls,list):
            self.ts_urls = [self.ts_urls]
        # TLE files
        self.earthsatellites = alm_conf_dict['EarthSatellites']
        loginf("timescale: %s, update interval: %.2f days" % ('builtin' if self.builtin else 'IERS file',self.update_interval/DAY_S))
        loginf("ephemeris file(s): %s" % self.eph_files)
        # Constellations
        if weeutil.weeutil.to_bool(alm_conf_dict.get('load_constellation_map',True)):
            self.init_constellation_map()
        # Used to inform the thread about shutdown
        self.evt = threading.Event()
        self.running = True
        # Variables to remember updates
        self.last_ts_update = 0
        self.last_eph_update = [0]*len(self.eph_files)
        self.last_sat_update = dict()
        logdbg("thread '%s': initialized" % self.name)
    
    def shutDown(self):
        """ shut down thread """
        self.running = False
        self.evt.set()
        loginf("thread '%s': shutdown requested" % self.name)
    
    def run(self):
        """ Skyfield database maintenance """
        stars_ok = False
        starnames_ok = False
        next_update = 0
        loginf("thread '%s': starting" % (self.name,))
        try:
            while self.running:
                # initialize list of star names (once at startup)
                if not starnames_ok:
                    starnames_ok = self.init_starnames()
                # initialize Skyfield or update its database
                if next_update<time.time():
                    success = self.init_skyfield()
                    logdbg("thread '%s': Initialization/update was%s successful." % (self.name,'' if success else ' not'))
                    if success: next_update = time.time()+self.update_interval
                # initialize/update earth satellite data (every day)
                sats_success = self.init_earth_satellites()
                # initialize star data (once at startup)
                if not stars_ok:
                    stars_ok = self.init_stars()
                # If no updating is required, the thread can be closed now.
                if not self.update_interval: break
                # Wait for the update interval to pass.
                self.evt.wait(86400 if success and sats_success and starnames_ok and stars_ok else 300)
        except Exception as e:
            logerr("thread '%s': %s - %s" % (self.name,e.__class__.__name__,e))
        finally:
            #ephemerides = None
            #for _eph in self.spk: _eph.close()
            loginf("thread '%s': stopped" % self.name)
    
    def init_constellation_map(self):
        """ load Skyfield's built-in constellation map """
        global constellation_at, constellation_names
        import skyfield.constellationlib
        constellation_at = skyfield.constellationlib.load_constellation_map()
        constellation_names = dict(skyfield.constellationlib.load_constellation_names())

    def init_starnames(self):
        """ init star names dictionary """
        global starnames, starids
        starids = named_star_dict
        starnames = {hip:name for name, hip in starids.items() }
        return True
        """
        try:
            fn = os.path.dirname(os.path.realpath(__file__))
            fn = os.path.join(fn,'starnames.dat')
            with open(fn,'rt') as f:
                for line in f:
                    x = line.split('|')
                    hip = weeutil.weeutil.to_float(x.pop(0))
                    if x:
                        starnames[hip] = x[0].strip()
                        if x[0]: starids[x[0].strip()] = hip
            loginf("thread '%s': successfully loaded starnames.dat" % self.name)
            return True
        except (OSError,TypeError,ValueError) as e:
            logerr("thread '%s': could not load '%s': %s - %s" % (self.name,fn,e.__class__.__name__,e))
            return False
        """
    
    def init_skyfield(self):
        """ download ephemeris data or read them from file """
        global ts, ephemerides, sun_and_planets, planets_list
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
                # Create timescale either from file or from builtin data
                _ts = load.timescale(builtin=self.builtin)
                if _ts:
                    s = 'initialized' if ts is None else 'updated'
                    ts = _ts
                    self.last_ts_update = time.time()
                    loginf("thread '%s': timescale %s" % (self.name,s))
                # install polar motion data
                if not self.builtin:
                    url = load.build_url(TIMESCALE_FILE)
                    with load.open(url) as f:
                        finals_data = iers.parse_x_y_dut1_from_finals_all(f)
                    iers.install_polar_motion_table(ts, finals_data)
                    loginf("thread '%s': installed polar motion table" % self.name)
            except OSError as e:
                logerr("thread '%s': error downloading timescale %s - %s" % (self.name,e.__class__.__name__,e))
        # load ephemerides
        ct = 0
        for idx, file in enumerate(self.eph_files):
            if self.last_eph_update[idx]<=now or self.spk[idx] is None:
                try:
                    # The last part of the path or URL is the file name.
                    file_name = file.split('/')[-1]
                    # Load the file
                    _eph = load(file)
                    if _eph:
                        self.spk[idx] = _eph
                        self.last_eph_update[idx] = time.time()
                        ct += 1
                except OSError as e:
                    logerr("thread '%s': error downloading ephemeris %s - %s" % (self.name,e.__class__.__name__,e))
        # merge ephemerides
        if ct and self.spk[0] is not None:
            _eph = dict()
            _sem = None
            for _spk in self.spk:
                # get the bodies supported by this ephemeris
                for nr,nm in _spk.names().items():
                    attr = nm[0].lower()
                    if attr not in _eph and nr!=0:
                        _eph[attr] = _spk[nr]
                # get the ephemeris set that covers sun, earth, and moon
                # required for seasons and moon phase calculation
                if (not _sem and 
                          SUN in _spk and EARTH in _spk and EARTHMOON in _spk):
                    _sem = _spk
            # Get a list of planets available in loaded ephemerides
            _pl = []
            for planet in PLANETS:
                if planet==EARTH:
                    # Because the coordinate systems used here refer to the
                    # Earth, it is excluded from the list.
                    continue
                if planet in _eph:
                    # The available ephemerides contain the planet itself.
                    # So include it here.
                    _pl.append(planet)
                elif ('%s_barycenter' % planet) in _eph:
                    # The available ephemerides contain the planet's barycentre 
                    # only. So use it instead.
                    _pl.append('%s_barycenter' % planet)
            # Is it the first run or some consecutive run?
            s = 'initialized' if ephemerides is None else 'updated'
            # Set the global variables
            ephemerides = _eph
            sun_and_planets = _sem
            planets_list = _pl
            loginf("thread '%s': ephemerides %s" % (self.name,s))
        # `ephemerides` and `ts` are up to date if they were updated less 
        # than 24 hours ago.
        return self.last_ts_update>now and min(self.last_eph_update)>now

    def init_earth_satellites(self):
        """ initialize earth satellites data
        
            Earth satellite ephemerides become out of date very fast.
        """
        global ephemerides, satcatalogues
        if ephemerides is None: return False
        # instantiate the loader
        load = Loader(self.path,verbose=False)
        rtn = True
        # get current time
        now = time.time()-86400
        # load earth satellites
        sats = dict()
        catalogues = []
        for file_name, url in self.earthsatellites.items():
            format = file_name.split('/')[-1].split('.')[-1].lower()
            if self.last_sat_update.get(file_name,0)<=now:
                try:
                    tmpfile = self.download([url],'%s.tmp' % file_name)
                    if tmpfile:
                        x = url.split('?')
                        x = x[1] if len(x)>1 else ""
                        for i in x.split('&'):
                            j = i.split('=')
                            if j[0].upper()=='FORMAT' and len(j)>1:
                                format = j[1].lower()
                                break
                        os.rename(tmpfile,os.path.join(self.path,file_name))
                        self.last_sat_update[file_name] = time.time()
                except OSError as e:
                    logerr("thread '%s': error downloading satellites file '%s' to '%s': %s - %s" % (self.name,url,file_name,e.__class__.__name__,e))
                    rtn = False
            try:
                with load.open(file_name) as f:
                    if format=='json':
                        data = json.load(f)
                        x = [EarthSatellite.from_omm(ts, fields) for fields in data]
                    elif format=='tle':
                        x = list(parse_tle_file(f,ts))
                    else:
                        raise TypeError('unknown file format %s' % format)
                catname = file_name.split('.')[0]
                catalogues.append(catname)
                x = {'%s_%s' % (catname,sat.model.satnum): sat for sat in x}
                sats.update(x)
                loginf("thread '%s': successfully processed satellites file '%s'" % (self.name,file_name))
            except (OSError,AttributeError,TypeError) as e:
                logerr("thread '%s': error installing satellites file %s - %s" % (self.name,e.__class__.__name__,e))
                rtn = False
        satcatalogues = set(catalogues)
        ephemerides.update(sats)
        loginf("thread '%s': %d earth satellite%s installed/updated" % (self.name,len(sats),'' if len(sats)==1 else 's'))
        return rtn
    
    def init_stars(self):
        """ load and initialize the stars dictionary
        
            See also skyfield.data.tycho2
        """
        global stars
        rtn = True
        # instantiate the loader
        load = Loader(self.path,verbose=False)
        # stars
        if has_pandas:
            import skyfield.data.hipparcos
            if not os.path.isfile(os.path.join(self.path, 'hip_main.dat')):
                try:
                    file = self.download([skyfield.data.hipparcos.URL],'hip_main.tmp')
                    if file:
                        os.rename(file,os.path.join(self.path,'hip_main.dat'))
                except (OSError,AttributeError,TypeError) as e:
                    logerr("thread '%s': error downloading star catalog %s - %s" % (self.name,e.__class__.__name__,e))
                    rtn = False
            try:
                with load.open('hip_main.dat') as f:
                    df = skyfield.data.hipparcos.load_dataframe(f)
                    stars = df[df['ra_degrees'].notnull()]
                    loginf("thread '%s': successfully installed star catalog" % self.name)
            except (OSError,AttributeError,TypeError) as e:
                logerr("thread '%s': error installing star catalog %s - %s" % (self.name,e.__class__.__name__,e))
                rtn = False
        return rtn
    
    def download(self, urls, filename=None):
        """ download file by FTP(S) or HTTP(S) """
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
        """ download file by HTTP(S) """
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
        sqlite_root = config_dict.get('DatabaseTypes',configobj.ConfigObj()).get('SQLite',configobj.ConfigObj()).get('SQLITE_ROOT','.')
        path = os.path.join(
            config_dict.get('WEEWX_ROOT','.'),
            sqlite_root,
            'skyfield'
        )
        if not os.path.isdir(path): os.mkdir(path)
        self.path = path
        self.default_almanac = None
        # configuration
        alm_conf_dict = _get_config(config_dict)
        if alm_conf_dict['enable']:
            # Is that run by `weectl report run`? If so, the maintenance
            # thread is stopped after initializing the ephemerides.
            is_interactive = SkyfieldService.is_interactive(engine)
            if is_interactive: alm_conf_dict['update_interval'] = 0
            # thread to initialize Skyfield
            self.skyfield_thread = SkyfieldMaintenanceThread(
                                                      alm_conf_dict, self.path)
            self.skyfield_thread.start()
            # in case of `weectl report run` wait for ephemerides to be
            # populated
            if is_interactive:
                loginf('waiting for ephemerides to be loaded')
                self.skyfield_thread.join()
                loginf('loading ephemerides finished')
            # instantiate the Skyfield almanac
            self.skyfield_almanac = SkyfieldAlmanacType()
            # add to the list of almanacs
            almanacs.insert(0,self.skyfield_almanac)
            logdbg("%s started" % self.__class__.__name__)
            # remove the PyEphem almanac if configured so
            if alm_conf_dict.get('disable_pyephem',False):
                if isinstance(almanacs[-1],PyEphemAlmanacType):
                    self.default_almanac = almanacs[-1]
                    del almanacs[-1]
                    loginf('Default almanac disabled')
        else:
            loginf("Skyfield almanac not enabled. Skipped.")
    
    def shutDown(self):
        """ remove this extension from the almanacs list and shut down the
            maintenance thread
        """
        global almanacs
        # if the PyEphem almanac was removed, re-add it here
        if self.default_almanac:
            almanacs.append(self.default_almanac)
        # find the Skyfield almanac in the list of almanacs
        idx = almanacs.index(self.skyfield_almanac)
        # remove it from the list
        del almanacs[idx]
        # stop thread
        if self.skyfield_thread.is_alive():
            self.skyfield_thread.shutDown()
        logdbg("%s stopped" % self.__class__.__name__)
    
    @staticmethod
    def is_interactive(engine):
        """ check whether this service is run by `weectl` or `weewxd`
        """
        return isinstance(engine,weewx.engine.DummyEngine)


class LiveService(StdService):
    """ calculate solar values for live update (optional service) 
    
        Some skins do not only update the web pages once every archive 
        interval but display live data out of LOOP packets. This service
        adds the necessary observation types for fast changing almanac 
        values to the LOOP packets and the ARCHIVE records for live updates.
        It also supports mobile stations.
    """

    def __init__(self, engine, config_dict):
        """ init """
        super(LiveService,self).__init__(engine, config_dict)
        # configuration
        alm_conf_dict = _get_config(config_dict)
        self.enabled = (alm_conf_dict['enable'] and 
            weeutil.weeutil.to_bool(alm_conf_dict.get('enable_live_data',True)))
        self.log_success = alm_conf_dict['log_success']
        self.log_failure = alm_conf_dict['log_failure']
        # station altitude
        try:
            self.altitude = weewx.units.convert(engine.stn_info.altitude_vt,'meter')[0]
        except (ValueError,TypeError,IndexError):
            self.altitude = 0.0
        loginf("Altitude %s ==> %.0f m" % (engine.stn_info.altitude_vt,self.altitude))
        # station location
        # Note: The default value of `elevation_m` is `0.0` not `None`.
        self.station = wgs84.latlon(engine.stn_info.latitude_f,
                                    engine.stn_info.longitude_f,
                                    elevation_m=self.altitude)
        # observation types
        weewx.units.obs_group_dict.setdefault('solarAltitude','group_angle')
        weewx.units.obs_group_dict.setdefault('solarAzimuth','group_direction')
        weewx.units.obs_group_dict.setdefault('solarTime','group_direction')
        weewx.units.obs_group_dict.setdefault('solarPath','group_percent')
        # instance variables
        self.last_archive_outTemp = None # degree_C
        self.last_archive_pressure = None # mbar
        self.last_almanac_error = 0
        self.last_mobile_error = 0
        self.sunrise = 0
        self.sunset = 0
        self.end_of_day = 0
        # register the methods with the kernel
        if self.enabled:
            self.bind(weewx.NEW_LOOP_PACKET, self.new_loop_packet)
            self.bind(weewx.NEW_ARCHIVE_RECORD, self.new_archive_record)
        logdbg("%s started" % self.__class__.__name__)
    
    def shutDown(self):
        """ shutdown service """
        logdbg("%s stopped" % self.__class__.__name__)
    
    def new_loop_packet(self, event):
        """ augment loop packet """
        # If there are a `latitude` and `longitude` observation type within
        # the LOOP packet we assume the station to be mobile and set the
        # current location.
        try:
            altitude_vt = weewx.units.as_value_tuple(event.packet,'altitude')
        except LookupError:
            altitude_vt = None
        try:
            self.set_current_location(
                weewx.units.as_value_tuple(event.packet,'latitude'),
                weewx.units.as_value_tuple(event.packet,'longitude'),
                altitude_vt
            )
        except LookupError:
            pass
        # calculate data and augment the LOOP packet
        self.calc_almanac(event.packet,False)
    
    def new_archive_record(self, event):
        """ augment archive record """
        try:
            val = weewx.units.as_value_tuple(event.record,'outTemp')
            self.last_archive_outTemp = weewx.units.convert(val,'degree_C')[0]
            logdbg("outTemp %s °C" % self.last_archive_outTemp)
        except (LookupError,ArithmeticError,AttributeError):
            pass 
        try:
            val = weewx.units.as_value_tuple(event.record,'pressure')
            self.last_archive_pressure = weewx.units.convert(val,'mbar')[0]
            logdbg("pressure %s mbar" % self.last_archive_pressure)
        except (LookupError,ArithmeticError,AttributeError):
            pass
        # calculate data and augment the ARCHIVE record
        self.calc_almanac(event.record,True)
    
    def calc_almanac(self, packet, archive):
        """ calculate solarAzimuth, solarAltitude, solarPath """
        global ephemerides
        # Do nothing until the Skyfield almanac ist initialized.
        if ephemerides is None: return
        try:
            sun = ephemerides[SUN]
            # target unit system
            usUnits = packet['usUnits']
            # current timestamp
            ts = packet.get('dateTime',time.time())
            ti = timestamp_to_skyfield_time(ts)
            # observer's location
            observer = ephemerides[EARTH] + self.station
            # apparent position of the sun in respect to the observer's location
            position = observer.at(ti).observe(sun).apparent()
            # solar altitude and azimuth
            alt, az, _ = position.altaz(temperature_C=self.last_archive_outTemp,pressure_mbar=self.last_archive_pressure if self.last_archive_pressure else 'standard')
            packet['solarAzimuth'] = weewx.units.convertStd(ValueTuple(az.degrees,'degree_compass','group_direction'),usUnits)[0]
            packet['solarAltitude'] = weewx.units.convertStd(ValueTuple(alt.radians,'radian','group_angle'),usUnits)[0]
            # solar time (hour angle)
            ha, _, _ = position.hadec()
            packet['solarTime'] = weewx.units.convertStd(ValueTuple(ha._degrees+180.0,'degree_compass','group_direction'),usUnits)[0]
            # solar path
            if archive and ((self.sunrise-8)<ha._degrees<self.sunrise or 
                                                           ts>self.end_of_day):
                if (self.last_archive_outTemp is not None and 
                                self.last_archive_pressure is not None):
                    horizon = -16.0/60.0-refraction(
                        -16.0/60.0,
                        temperature_C=self.last_archive_outTemp,
                        pressure_mbar=self.last_archive_pressure
                    )
                else:
                    horizon = None
                start_ts, end_ts = weeutil.weeutil.archiveDaySpan(ts)
                start_ti = timestamp_to_skyfield_time(start_ts)
                end_ti = timestamp_to_skyfield_time(end_ts)
                if SKYFIELD_VERSION<(1,47):
                    f = almanac.risings_and_settings(sun_and_planets, sun, self.station)
                    t, y = almanac.find_discrete(start_ti, end_ti, f)
                    h, _, _ = observer.at(t).observe(sun).apparent().hadec()
                    for i,j in zip(h._degrees,y):
                        if j==1:
                            self.sunrise = i
                        elif j==0:
                            self.sunset = i
                else:
                    t, y = almanac.find_risings(observer,sun,start_ti,end_ti,horizon_degrees=horizon)
                    if y[0]:
                        h, _, _ = observer.at(t[0]).observe(sun).apparent().hadec()
                        self.sunrise = h._degrees
                    else:
                        self.sunrise = 0
                    t, y = almanac.find_settings(observer,sun,start_ti,end_ti,horizon_degrees=horizon)
                    if y[0]:
                        h, _, _ = observer.at(t[0]).observe(sun).apparent().hadec()
                        self.sunset = h._degrees
                    else:
                        self.sunset = 0
                self.end_of_day = end_ts
            if (self.sunrise and self.sunset and 
                          self.sunset>self.sunrise and 
                                       self.sunrise<=ha._degrees<=self.sunset):
                sp = (ha._degrees-self.sunrise)/(self.sunset-self.sunrise)*100.0
            else:
                sp = None
            packet['solarPath'] = weewx.units.convertStd(ValueTuple(sp,'percent','group_percent'),usUnits)[0]
        except (LookupError,ArithmeticError,AttributeError,TypeError,ValueError) as e:
            # report the error at most once every 5 minutes
            if self.log_failure and time.time()>=self.last_almanac_error+300:
                logerr("live almanac error: %s" % e)
                self.last_almanac_error = time.time()
    
    def set_current_location(self, latitude_vt, longitude_vt, altitude_vt):
        """ set current location in case of a mobile station 

            If there are `latitude` and `longitude` observation types
            in the LOOP packet, we assume the station to be mobile
            and set the station's location to the new values.
        """
        try:
            if latitude_vt is not None and longitude_vt is not None:
                # altitude
                try:
                    alt = weewx.units.convert(latitude_vt,'meter')[0] if altitude_vt else 0.0
                except (LookupError,ArithmeticError,AttributeError,TypeError,ValueError):
                    alt = 0.0
                # latitude
                # Note: From WeeWX 5.0 on there is a new unit `degree_angle`
                #       and a new unit group `group_angle`, which are used
                #       for altitude and declination, and thus would be
                #       appropriate for latitude, too. Nevertheless, for
                #       reasons of backward compatibility the old unit
                #       `degree_compass` has also to be supported.
                if latitude_vt[1]=='degree_compass':
                    lat = latitude_vt[0]
                else:
                    lat = weewx.units.convert(latitude_vt,'degree_angle')[0]
                # longitude
                lon = weewx.units.convert(longitude_vt,'degree_compass')[0]
                if lat is not None and lon is not None:
                    # Note: The default value of `elevation_m` is `0.0` not
                    #       `None`.
                    self.station = wgs84.latlon(lat,lon,elevation_m=alt)
        except (LookupError,ArithmeticError,AttributeError,TypeError,ValueError) as e:
            # report the error at most once every 5 minutes
            if self.log_failure and time.time()>=self.last_mobile_error+300:
                logerr("could not set new location: %s - %s" % (e.__class__.__name__,e))
                self.last_mobile_error = time.time()
        
    
# log version info at startup
loginf("%s version %s" % ("WeeWX Skyfield almanac extension",VERSION))
loginf("Skyfield version %s" % '.'.join([str(i) for i in SKYFIELD_VERSION]))
