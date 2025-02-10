# weewx-skyfield-almanac
Almanac extension to WeeWX using Skyfield module

## Prerequisites

Install Skyfield

## Installation instructions

1) download

   ```shell
   wget -O weewx-skyfield-almanac.zip https://github.com/roe-dl/weewx-skyfield-almanac/archive/master.zip
   ```

2) run the installer

   WeeWX from version 5.1 on and WeeWX packet installation

   ```shell
   sudo weectl extension install weewx-skyfield-almanac.zip
   ```

   WeeWX from version 5.1 on and WeeWX pip installation into an virtual environment

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
    # update interval 1 year (set to 0 for no update)
    update_interval = 31557600
```

* `ephemeris`: Ephemeris file to use. Different files cover different
  scopes of heavenly bodies. Some of those files are huge. See
  [Planets and their moons: JPL ephemeris files](https://rhodesmill.org/skyfield/planets.html)
  for details.
* `use_builtin_timescale`: Use builtin timescale data or download them
  from IERS. See [UT1 and downloading IERS data](https://rhodesmill.org/skyfield/time.html#ut1-and-downloading-iers-data)
  for details.
* `update_interval`: interval for updating ephemeris and timescale data
  (set to 0 to switch off updates)

## Links

* [WeeWX](https://weewx.com)
* [Skyfield](https://rhodesmill.org/skyfield/)
* [Issue #981: PyEphem is deprecated](https://github.com/weewx/weewx/issues/981)
