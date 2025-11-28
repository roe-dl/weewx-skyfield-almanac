# installer almanac extension
# Copyright 2025 Johanna Roedenbeck
# Distributed under the terms of the GNU Public License (GPLv3)

import os.path
import configobj
from weecfg.extension import ExtensionInstaller
from weeutil.config import merge_config

def loader():
    return SkyfieldInstaller()

class SkyfieldInstaller(ExtensionInstaller):
    def __init__(self):
        super(SkyfieldInstaller, self).__init__(
            version="0.3",
            name='Skyfield almanac',
            description='almanac extension using Skyfield mdule',
            author="Johanna Roedenbeck",
            author_email="",
            prep_services='user.skyfieldalmanac.SkyfieldService',
            data_services='user.skyfieldalmanac.LiveService',
            config={
                'Almanac': {
                    'Skyfield': {
                        'enable':'true',
                        'ephemeris':'de440.bsp',
                        'use_builtin_timescale':'true',
                        'timescale_url': [
                            'https://datacenter.iers.org/products/eop/rapid/standard/finals2000A.all',
                            'ftps://gdc.cddis.eosdis.nasa.gov/products/iers/finals.all'
                        ],
                        '#log_ftp':'false',
                        'update_interval':'31557600',
                        'enable_live_data':'true'
                    }
                }
            },
            files=[('bin/user', ['bin/user/skyfieldalmanac.py'])]
        )

    def configure(self, engine):
        """ special configuration """
        # WEEWX_ROOT, USER_DIR, EXT_DIR, BIN_DIR, probably SKIN_DIR
        skin_dir = engine.root_dict.get('SKIN_DIR')
        user_dir = engine.root_dict.get('USER_DIR')
        extension_dir = os.path.dirname(__file__)
        extension_lang_dir = os.path.join(extension_dir,'lang')
        if not os.path.isdir(extension_lang_dir):
            engine.printer.out('directory %s not present. Skip injecting translations' % extension_lang_dir)
        elif skin_dir and user_dir and 'StdReport' in engine.config_dict:
            # if present update the files if possible
            for skin in engine.config_dict['StdReport'].sections:
                # relative path to the skin
                skin_pth = engine.config_dict['StdReport'][skin].get('skin')
                # absolute path to the language files of the skin
                lang_dir = os.path.join(skin_dir,skin_pth,'lang') if skin_pth else None
                if lang_dir and os.path.isdir(lang_dir):
                    engine.printer.out('processing skin %s' % skin)
                    for fn in os.listdir(extension_lang_dir):
                        if fn.endswith('.conf') and fn!='lang.conf':
                            src_fn = os.path.join(extension_lang_dir,fn)
                            dest_fn = os.path.join(lang_dir,fn)
                            # check if the source file is a regular file
                            if not os.path.isfile(src_fn):
                                engine.printer.out('file %s not found' % src_fn)
                                continue
                            # The standardized language code for czech is 
                            # 'cs', but in WeeWX it is 'cz'. Both are tried.
                            if fn=='cs.conf':
                                if os.path.isfile(dest_fn):
                                    self._update_lang_file(engine, src_fn, dest_fn)
                                dest_fn = os.path.join(lang_dir,'cz.conf')
                            # check if the target file exists and update it
                            if os.path.isfile(dest_fn):
                                self._update_lang_file(engine, src_fn, dest_fn)
        # return whether changes to the configuration file were done
        return False

    def _update_lang_file(self, engine, src_fn, dest_fn):
        """ Update language definition files """
        # get the original language file
        try:
            config = configobj.ConfigObj(dest_fn,encoding='utf-8')
        except configobj.ConfigObjError as e:
            engine.printer.out('cannot merge to %s: %s %s' % (dest_fn,e.__class__.__name__,e))
            return
        # get the Skyfield additions
        to_be_merged = configobj.ConfigObj(src_fn,encoding='utf-8')
        # merge the additions to the localization config
        if 'Almanac' in config:
            engine.printer.out('merging %s to %s' % (os.path.basename(src_fn),dest_fn))
            merge_config(config['Almanac'],to_be_merged['Almanac'])
            if ('Astronomical' in config.get('Texts',dict()) and 
                'Astronomical' in to_be_merged.get('Texts',dict())):
                merge_config(config['Texts']['Astronomical'],to_be_merged['Texts']['Astronomical'])
            # save
            if engine.dry_run:
                engine.printer.out(config)
                engine.printer.out('-'*72)
            else:
                config.write()
