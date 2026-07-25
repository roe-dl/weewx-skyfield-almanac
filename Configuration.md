# Configuration instructions

There is no need to configure anything, but there are some tuning options
available if you have special requirements.

> [!NOTE]
> If you are new to this extension, please, do NOT change the configuration
> at the beginning. The installation sets reasonable values to all 
> configuration keys to immediately use this extension.

```
[Almanac]
    [[Skyfield]]
        # use this almanac
        enable = true
        # which ephemeris file to use
        ephemeris = de440s.bsp  # or de440.bsp or de441.bsp
        # list planetary constants files
        planetaryconstants = ...
        # use builtin timescale data or download it from IERS
        use_builtin_timescale = true
        # URL(s) of the timescale file (optional)
        timescale_url = '...'
        # whether to load Skyfield's built-in constellation map
        load_constellation_map = true
        # whether to log FTP responses (optional)
        log_ftp = false
        # update interval 1 year (set to 0 for no update)
        update_interval = 31557600
        # enable LOOP packet augmentation
        enable_live_data = true
        # which observation types to calculate live data for
        live_data_observations = altitude, azimuth
        # optional list of heavenly bodies
        live_data_bodies = 
        # disable the built-in PyEphem almanac (optional)
        disable_pyephem = false
        [[[EarthSatellites]]]
            file_name1 = url1
            file_name2 = url2
            ...
        [[[Frames]]]
            heavenly_body_name1 = frame_name1
            heavenly_body_name2 = frame_name2
            ...
```

* `enable`: Enable this almanac extension
* `ephemeris`: Ephemeris (or SPICE kernel) file or a list of such files to use.
  Different files cover different scopes of heavenly bodies and/or different
  time spans. Some of those files are huge. See
  [Ephemerides](#ephemerides) and
  [Planets and their moons: JPL ephemeris files](https://rhodesmill.org/skyfield/planets.html)
  for more details. 
  The base kernel file is for example  `de440s.bsp`, `de440.bsp`, or `de441.bsp`.
  This covers the Sun, the Moon, and all the planets. You do not need to
  specify the full URL of the base kernel, as that is built into Skyfield.
  An additional file is for example
  `https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup365.bsp`
  which covers some well-known moons of Jupiter.
* `planetaryconstants`: List of URLs to `.tf`, `.tpc`, and `.bpc` files.
  Optional. Required for reference frames only.
* `use_builtin_timescale`: Use builtin timescale data or download them
  from IERS. See [UT1 and downloading IERS data](https://rhodesmill.org/skyfield/time.html#ut1-and-downloading-iers-data)
  for details.
* `timescale_url`: an URL or a list of URLs to download the timescale file from (optional). 
  There is a default URL hardcoded in Skyfield. Unfortunately the server
  is temporarily down. That's why you can specify an alternative
  source here.
* `load_constellation_map`: Whether to load Skyfield's built-in
  constellation map (optional, default `True`)
* `log_ftp`: whether to log FTP responses of the server (optional).
  If you specified an alternative source for the timescale file in 
  `timescale_url` and that URL is at an FTP server, you can switch
  on logging of the server responses in case of trouble.
* `update_interval`: interval for updating ephemerides and timescale data
  (set to 0 to switch off updates)
* `enable_live_data`: enable live data for fast changing almanac values
  (default: on)
* `live_data_observations`: list of observation types to calculate live
  data for. Optional. Default is altitude, azimuth, and distance. Possible 
  additional values are `declination`, `right ascension`, and `libration`.
* `live_data_bodies`: list of additional heavenly bodies to include in live 
  data (e.g. LOOP packets). Optional. The Sun is always included.
* `disable_pyephem`: disable the built-in PyEphem almanac 
  (optional, default: false).
  Sometimes enabling both the Skyfield and PyEphem almanac produces 
  special errors if the attribute or heavenly body is not available.
  In those cases setting this option may help.
* `[[[EarthSatellites]]]`: This section contains earth satellite data files
  to load. Each entry in the section contains of a file name and an URL
  The file name is used when saved to disk. You can find such files
  for example at [CelesTrak](https://celestrak.org/NORAD/elements/).
  See section [Earth satellites](#earth-satellites) for more datails.
* `[[[Frames]]]`: Planetary reference frames. Optional.
  Used to get a location on the surface of a heavenly body, for example
  the Moon.
  Actually required for the libration and coordinate axis calculation only.
  See
  [Skyfield: Planetary Reference Frames](https://rhodesmill.org/skyfield/planetary.html).
  for more details.

