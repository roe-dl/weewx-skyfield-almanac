# weewx-skyfield-almanac
Almanac extension to WeeWX using Skyfield module

## Why should I use this extension?

PyEphem is deprecated. Its astronomical database is outdated and won't get
updated any more. It ends in 2018. Dates after that year are calculated
by extrapolation.

Skyfield is the successor of PyEphem. It is from the same author, Brandon
Rhodes. It uses more modern and more precise formulae and actual ephemeris
provided by NASA's JPL.

There is no other requirement than installing this extension to replace
PyEphem calculated values by Skyfield calculated values in existing skins.

## Prerequisites

WeeWX from version 5.2 on is required.

Install Skyfield and NumPy

## Installation instructions

1) download

   ```shell
   wget -O weewx-skyfield-almanac.zip https://github.com/roe-dl/weewx-skyfield-almanac/archive/master.zip
   ```

2) run the installer

   WeeWX from version 5.2 on and WeeWX packet installation

   ```shell
   sudo weectl extension install weewx-skyfield-almanac.zip
   ```

   WeeWX from version 5.2 on and WeeWX pip installation into an virtual environment

   ```shell
   source ~/weewx-venv/bin/activate
   weectl extension install weewx-skyfield-almanac.zip
   ```
   
3) restart weewx

   for SysVinit systems:

   ```shell
   sudo /etc/init.d/weewx stop
   sudo /etc/init.d/weewx start
   ```

   for systemd systems:

   ```shell
   sudo systemctl stop weewx
   sudo systemctl start weewx
   ```

## Configuration instructions

There is no need to configure anything, but there are some tuning options
available if you have special requirements.

```
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
    # enable LOOP packet augmentation
    enable_live_data = true
```

* `ephemeris`: Ephemeris file to use. Different files cover different
  scopes of heavenly bodies. Some of those files are huge. See
  [Planets and their moons: JPL ephemeris files](https://rhodesmill.org/skyfield/planets.html)
  for details.
* `use_builtin_timescale`: Use builtin timescale data or download them
  from IERS. See [UT1 and downloading IERS data](https://rhodesmill.org/skyfield/time.html#ut1-and-downloading-iers-data)
  for details.
* `timescale_url`: an URL or a list of URLs to download the timescale file from (optional). 
  There is a default URL hardcoded in Skyfield. Unfortunately the server
  is temporarily down. That's why you can specify an alternative
  source here.
* `log_ftp`: whether to log FTP responses of the server (optional).
  If you specified an alternative source for the timescale file in 
  `timescale_url` and that URL is at an FTP server, you can switch
  on logging of the server responses in case of trouble.
* `update_interval`: interval for updating ephemeris and timescale data
  (set to 0 to switch off updates)
* `enable_live_data`: enable live data for fast changing almanac values
  (default: on)

## Usage

There is no other requirement than installing this extension to replace
PyEphem calculated values by Skyfield calculated values in existing skins.

Once the weewx-skyfield-almanac extension is installed and initialized after
startup, `$almanac.hasExtras` becomes `True` and extended almanac
information is available. Initialization can take several archive
intervals to be completed at first run after installation, depending on 
configuration.

## Customization of WeeWX using this extension

See the WeeWX customization guide, section "The Cheetah generator",
sub-section
"[Almanac](http://weewx.com/docs/latest/custom/cheetah-generator/#almanac)",
for a detailed description how to use the almanac in WeeWX. This section
repeats some of that information and adds, what is specific to this
extension.

The general syntax is:

```
$almanac(almanac_time=time,            ## Unix epoch time
         lat=latitude, lon=longitude,  ## degrees
         altitude=altitude,            ## meters
         pressure=pressure,            ## mbars
         horizon=horizon,              ## degrees
         temperature=temperature_C     ## degrees C
       ).heavenly_body(use_center=[01]).attribute
```

If `almanac_time` is not specified, the actual time as returned by
`$current.dateTime` is used.

If `lat` and `lon` are not specified, the location of the station is used.

### Date and time

Date and time values do not require a heavenly body. They refer to the
location on earth that `$almanac` is bound to. You can specify the
location by parameters to `$almanac`. If you do not set them, the
location of the station is used.

WeeWX datatype   | Pure float result | Meaning
-----------------|-------------------|----------------
`sidereal_angle` | `sidereal_time`   | Local Apparent Sidereal Time
`solar_angle`    | `solar_time`      | Local Apparent Solar Time
`solar_datetime` | &mdash;           | true local solar date and time

Note: The tags `$almanac.sidereal_angle`, `$almanac.sidereal_time`,
`$almanac.solar_angle`, and `$almanac.solar_time` return a value
in decimal degrees rather than the more customary value from 0 to
24 hours.

Solar time is the time a sundial would show.

If you do not install this extension together with Skyfield but use PyEphem,
which is supported by core WeeWX, `sidereal_angle` and `sidereal_time` are
available only.

### Calendar events

Calendar events do not require a heavenly body. They refer to earth.
This extension provides the events described in the WeeWX customization
guide, but calculated using Skyfield. They happen independent of your
location on earth at the same instant all over the world. The local time
differs only. Here is a list of available events:

Previous event | Next event |
---------------|------------|
`previous_equinox` | `next_equinox`
`previous_solstice` | `next_solstice`
`previous_autumnal_equinox` | `next_autumnal_equinox`
`previous_vernal_equinox` | `next_vernal_equinox`
`previous_winter_solstice` | `next_winter_solstice`
`previous_summer_solstice` | `next_summer_solstice`
`previous_new_moon` | `next_new_moon`
`previous_first_quarter_moon` | `next_first_quarter_moon`
`previous_full_moon` | `next_full_moon`
`previous_last_quarter_moon` | `next_last_quarter_moon`

Example:
```
$almanac.next_solstice
```

### Heavenly bodies

This extension provides the attributes described in the WeeWX customization 
guide, but calculated using Skyfield. Additionally it provides some
extra attributes, that are not available with PyEphem.

Depending on the ephemeris you chose you may be required to add
`_barycenter` to the name of a heavenly body to get results
(for example `jupiter_barycenter`).

These events are supported for heavenly bodies in reference to the 
location and the timestamp specified:

Day of timestamp | Previous event | Next event | Meaning
------|----------------|------------|----------
`rise` | `previous_rising` | `next_rising` | rising of the body above the horizon
`transit` | `previous_transit` | `next_transit` | when the body transits the meridian
`set` | `previous_setting` | `next_setting` | setting of the body above the horizon
`antitransit` | `previous_antitransit` | `next_antitransit` | antitransit
`visible` | &mdash; | &mdash; | how long the body will be visible
`visible_change` | &mdash; | &mdash; | change in visbility compared to previous day

Example:
```
$almanac.sun.rise
```

Here is the list of attributes provided by this extension but not by
core WeeWX using PyEphem:

WeeWX datatype | Pure float result | Meaning
---------------|-------------------|----------------
`astro_dist`   | `a_dist`          | astrometric geocentric distance
`geo_dist`     | `g_dist`          | apparent astrometric geocentric distance
`topo_dist`    | `dist`            | apparent topocentric distance 
`alt_distance` | `alt_dist`        | distance in reference to the coordinate system of altitude and azimuth
`hour angle`   | `ha`              | topocentric hour angle
`ha_declination` | `ha_dec`        | declination in reference to the coordinate system of the hour angle
`ha_distance`  | `ha_dist`         | distance in referenc to the coordinate system of the hour angle
`max_altitude` | `max_alt`         | maximum altitude
`max_alt_time` | &mdash;           | timestamp of the maximum altitude

And these attributes are supported by both core WeeWX using PyEphem and
this extension using Skyfield:

WeeWX datatype | Pure float result | Meaning
---------------|-------------------|----------------
`azimuth` | `az` | apparent azimuth of the body in the sky
`altitude` | `alt` | apparent altitude of the body in the sky
`astro_ra` | `a_ra` | astrometric geocentric right ascension
`astro_dec` | `a_dec` | astrometric geocentric declination
`geo_ra` | `g_ra` | apparent geocentric right ascension
`geo_dec` | `g_dec` | apparent geocentric declination
`topo_ra` | `ra` | apparent topocentric right ascension
`topo_dec` | `dec` | apparent topocentric declination

### Coordinate systems

There are different coordinate systems used to express locations.

Within base plane | Rectangular to it | Base plane  | Origin | Direction
-----------|--------------|-------------|--------|-------------
longitude (-180°...+180°) | latitude (-90°...+90°) | earth's equator | earth's center | meridian of Greenwich
azimuth (0°...360°) | altitude (-90°...+90°) | horizon of the observer | observer | geographic north
right ascension (0°...360°) | declination (-90°...+90°) | earth's equator | earths's center | sun at spring equinox
hour angle (-180°...+180°) | declination (-90°...+90°) | earth's equator | earth's center | observer's meridian

Please note, that the earth's equator also moves in different ways, and
the attributes differ in which of the movements they take care of.

And those changes are also the reason you have to provide a date, called
*epoch*, along with the coordinates to fully specify a location in space.
For calculating apparent geocentric and topocentric right ascension and
declination this extension uses `almanac_time` for the epoch.

The hour angle and its declination are calculated using polar motion data
if available. `ha_declination` (`ha_dec`) is thus slightly different from
`topo_dec` (`dec`).

## PyEphem and Skyfield

If you install both PyEphem and Skyfield, Skyfield is preferred. If the
given heavenly body is available with Skyfield, the attribute is calculated
using Skyfield. Otherwise PyEphem is tried. If neither Skyfield nor 
PyEphem know about the body, an exception is raised.

## Live data

Some skins do not only update the web pages once every archive interval
but present live data out of LOOP packets. This extension can add fast 
changing almanac values to the LOOP packet for live updates. To use them
include them in MQTT output. To activate live data, set the
`enable_live_data` configuration option to `true`, which is the default.

* `solarAltitude`: current altitude of the sun
* `solarAzimuth`: current azimuth of the sun
* `solarTime`: current hour angle of the sun, counted from midnight
  (that is sundial time)
* `solarPath`: current percentage of the way of the sun from sunrise
  to sunset

If you want to use those values for further calculations - for example
for sunshine duration calculation - make sure to put 
`user.skyfieldalmanac.LiveService` in the list of data services before 
the module that uses the values.

For mobile stations the location is updated from `latitude` and `longitude`
observation types if they are present in the LOOP packet.

## FAQ

Q: I installed the Skyfield module, but no extended almanac values
   are displayed.

A: Depending on the ephemeris file you chose initialization takes more or
   less time. While it is not finished, the skin displays some core
   values only. Wait a moment, and the extended values will show up.
   If the ephemeris file is huge and it is the first run after
   installation, it can take several archive intervals to be completed.


Q: Which ephemeris file should I use?

A: The ephemeris files differ in size and coverage of heavenly bodies and
   time. The default, if you do not choose anything, is `de440s.bsp`. It
   covers sun, earth, moon, and some data of the long known planets. If
   you need more, try `de440.bsp` or even `de441.bsp`.
   From time to time the NASA releases new sets of ephemeris files. Then you
   can try the new edition.
   See [Planets and their moons: JPL ephemeris files](https://rhodesmill.org/skyfield/planets.html)
   for more details.


Q: Are there disadvantages of Skyfield?

A: Yes. Always are. Skyfield depends on NumPy, while PyEphem does not.


Q: How do I have to adapt my skin to use this extension?

A: There is nothing to do. Installing this extension is enough. But you
   could add the additional attributes supported by this extension
   to your skin to display them there.


## Links

* [WeeWX](https://weewx.com)
* [Skyfield](https://rhodesmill.org/skyfield/)
* [International Earth Rotation and Reference Service IERS](https://iers.org)
  (provides the timescale file `finals2000A.all`)
* [Jet Propulsion Laboratory JPL](https://www.jpl.nasa.gov)
  (provides the ephemeris files)
* [Issue #981: PyEphem is deprecated](https://github.com/weewx/weewx/issues/981)
* [International Terrestrial Reference System ITRS](https://en.wikipedia.org/wiki/International_Terrestrial_Reference_System_and_Frame)
