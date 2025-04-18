# What is the difference between PyEphem based and Skyfield based WeeWX almanac?

See [README](README.md) for detailled information and manual. This is just a
summary.

## History

Both PyEphem and Skyfield are written by the same author, Brandon Rhodes.
Skyfield is the successor of PyEphem.

## Formulae

PyEphem uses formulae as published by Jean Meeus in the 1980s.

Skyfield uses up to date formulae as provided by the NASA.

## Source of timescales and ephemerides

PyEphem includes built-in data which are not updated any more.

Skyfield is based on downloadable files which are updated regularly. A lot
of different files are available covering different kinds of heavenly objects
including earth satellites and stars and different timespans.

## Additions

As there is sidereal time, there is also solar time. Solar time is more
related to the daily course of the weather than civil or timezone time.

Like the seasons the apsides (perihelion, aphelion, perigee, apogee) are 
events during the course of the year or the month, respectively.

Like the Moon Venus and Mercury show phases. There is a current phase, and
there are distinct events.

## Almanac attributes

The following tables only include additional attributes available with the 
Skyfield based almanac but not with the PyEphem based one.

### General attributes

Attribute | Data type | Meaning
----------|-----------|--------
`venus_index` | int | venus phase index (0 to 7)
`venus_phase` | str | name of the actual venus phase
`mercury_index` | int | mercury phase index (0 to 7)
`mercury_phase` | str | name of the actual mercury phase

### Date and time

WeeWX datatype   | Pure float result | Meaning
-----------------|-------------------|----------------
`solar_angle`    | `solar_time`      | Local Apparent Solar Time
`solar_datetime` | &mdash;           | true local solar date and time

### Calendar events

Previous event | Next event | Meaning
---------------|------------|------------------
`previous_perihelion` | `next_perihelion` | perihelion of the Earth (when the Earth is nearest to the Sun)
`previous_aphelion` | `next_aphelion` | aphelion of the Earth (when the Earth is farthest from the Sun)
`previous_perigee_moon` | `next_perigee_moon` | perigee of the Moon (when the Moon is nearest to the Earth; in connection with full moon sometimes "supermoon")
`previous_apogee_moon` | `next_apogee_moon` | apogee of the Moon (when the Moon is farthest from the Earth)
`previous_new_venus` | `next_new_venus` | maximum of phase angle; Venus changes from evening to morning side
`previous_first_quarter_venus` | `next_first_quarter_venus` | waxing 90 degrees of phase angle
`previous_full_venus` | `next_full_venus` | minimum of phase angle
`previous_last_quarter_venus` | `next_last_quarter_venus` | waning 90 degrees of phase angle

### Heavenly bodies

WeeWX datatype | Pure float result | Meaning
---------------|-------------------|----------------
`astro_dist`   | `a_dist`          | astrometric geocentric distance
`geo_dist`     | `g_dist`          | apparent astrometric geocentric distance
`topo_dist`    | `dist`            | apparent topocentric distance 
`alt_distance` | `alt_dist`        | distance in reference to the coordinate system of altitude and azimuth
`hour angle`   | `ha`              | topocentric hour angle
`ha_declination` | `ha_dec`        | declination in reference to the coordinate system of the hour angle
`ha_distance`  | `ha_dist`         | distance in referenc to the coordinate system of the hour angle
`day_max_altitude` | `day_max_alt` | maximum altitude of the day
`day_max_alt_time` | &mdash;       | timestamp of the maximum altitude of the day
`moon_tilt` | &mdash; | crescent moon tilt angle (0 = illuminated side to the right on northern hemisphere, π = illuminated side to the left on northern hemisphere
&mdash; | `hip_number` | in case of stars the Hipparcos catalog number
&mdash; | `venus_fullness` | percentage of Venus that is illuminated
&mdash; | `mercury_fullness` | percentage of Mercury that is illuminated
&mdash; | `constellation_abbr` | abbrevation of the constellation the actual position of the body is in
&mdash; | `constellation` | name of the constellation the actual position of the body is in

### Maps

An extra WeeWX almanac extension provides a sky map and a moon phase picture
based on the WeeWX Skyfield almanac extension.

### Live data

Live data differ from almanac attributes in that they can be output by MQTT
for live website updates and can be saved to database for creating diagrams.
On the other side they cannot have parameters.

Observation type | Meaning
-----------------|-----------
`solarAltitude` | current altitude of the Sun
`solarAzimuth` | current azimuth of the Sun
`solarTime` | current hour angle of the Sun, counted from midnight (that is sundial time)
`solarPath` | current percentage of the way of the Sun from sunrise to sunset
