# installer almanac extension
# Copyright 2025 Johanna Roedenbeck
# Distributed under the terms of the GNU Public License (GPLv3)

from weecfg.extension import ExtensionInstaller

def loader():
    return SkyfieldInstaller()

class SkyfieldInstaller(ExtensionInstaller):
    def __init__(self):
        super(SkyfieldInstaller, self).__init__(
            version="0.1",
            name='Skyfield almanac',
            description='almanac extension using Skyfield mdule',
            author="Johanna Roedenbeck",
            author_email="",
            prep_services='user.skyfieldalmanac.SkyfieldService',
            config={
                'Almanac': {
                    'ephemeris':'de440.bsp',
                    'use_builtin_timescale':'true',
                    'timescale_url': [
                        'https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all',
                        'ftps://gdc.cddis.eosdis.nasa.gov/products/iers/finals.all'
                    ],
                    '#log_ftp':'false',
                    'update_interval':'31557600'
                }
            },
            files=[('bin/user', ['bin/user/skyfieldalmanac.py'])]
        )
