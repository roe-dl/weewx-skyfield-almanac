# What is the difference between PyEphem based and Skyfield based WeeWX almanac?

See [README](README.md) for detailled information and manual. This is a
summary here.

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
&mdash; | `hip_number` | in case of stars the Hipparcos catalog number
&mdash; | `venus_fullness` | percentage of Venus that is illuminated
&mdash; | `mercury_fullness` | percentage of Mercury that is illuminated

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
