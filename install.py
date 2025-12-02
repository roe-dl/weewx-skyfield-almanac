# installer almanac extension
# Copyright 2025 Johanna Roedenbeck
# Distributed under the terms of the GNU Public License (GPLv3)

import os.path
import configobj
from weecfg.extension import ExtensionInstaller
from weeutil.config import conditional_merge

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
        self.is_download_ephemeris = False
        self.is_skip_localization = False
    
    def process_args(self, args):
        """ process args passed to the installer 
        
            Note: Recognized from WeeWX 5.3 on
        """
        for arg in args:
            if arg=='--download-ephemeris':
                self.is_download_ephemeris = True
            if arg=='--skip-localization':
                self.is_skip_localization = True

    def configure(self, engine):
        """ special configuration """
        # WEEWX_ROOT, USER_DIR, EXT_DIR, BIN_DIR, probably SKIN_DIR
        skin_dir = engine.root_dict.get('SKIN_DIR')
        user_dir = engine.root_dict.get('USER_DIR')
        extension_dir = os.path.dirname(__file__)
        extension_lang_dir = os.path.join(extension_dir,'lang')
        if self.is_skip_localization:
            engine.printer.out('Skip localization due to user request')
        elif not os.path.isdir(extension_lang_dir):
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
        # Download Ephemerides
        if self.is_download_ephemeris:
            data_dir = engine.config_dict.get('DatabaseTypes',
                                dict()).get('SQLite',dict()).get('SQLITE_ROOT')
            if data_dir:
                data_dir = os.path.join(
                    engine.root_dict.get('WEEWX_ROOT','.'),
                    data_dir,
                    'skyfield'
                )
                self.download_ephemerides(engine, data_dir)
            else:
                engine.printer.out('could not determine database directory')
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
            conditional_merge(config['Almanac'],to_be_merged['Almanac'])
            if ('Astronomical' in config.get('Texts',dict()) and 
                'Astronomical' in to_be_merged.get('Texts',dict())):
                conditional_merge(config['Texts']['Astronomical'],to_be_merged['Texts']['Astronomical'])
            # save
            if engine.dry_run:
                engine.printer.out(config)
                engine.printer.out('-'*72)
            else:
                config.write()
    
    def download_ephemerides(self, engine, data_dir):
        """ download ephemerides according to configuration
        
            Note: This is only required if your WeeWX installation does
                  not have permanent Internet access.
        """
        try:
            from skyfield.iokit import Loader
        except ImportError:
            engine.printer.out('cannot download ephemerides: Skyfield module not available')
            return
        alm_dict = engine.config_dict['Almanac']['Skyfield']
        engine.printer.out('download ephemerides to %s' % data_dir)
        load = Loader(data_dir)
        eph_files = alm_dict.get('ephemeris','de440s.bsp')
        if not isinstance(eph_files,list): eph_files = [eph_files]
        for eph_file in eph_files:
            eph = load(eph_file)
            engine.printer.out(eph)
