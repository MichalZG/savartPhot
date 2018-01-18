from configparser import ConfigParser
import os
from astropy.io import fits


class Star:
	def __init__(self, name):
		self.name = name
		self.files = []


class Image:
    def __init__(self, image_dir, datetime_key, jd_key, filter_key, object_key):
        self.image_dir = image_dir
        self.hdr = self._get_header()
        self.name = os.path.basename(image_dir)
        self._find_hdr_fields(datetime_key, jd_key, filter_key, object_key)

    def _find_filter(self, filter_name):

        if not any(f in filter_name for f in ['P1', 'P2', 'P3', 'P4']):
            raise ValueError('Savart name in filter key has not been found')

        return filter_name.split('-')

    def _find_hdr_fields(self, datetime_key, jd_key, filter_key, object_key):

        try:
            self.datetime = self.hdr[datetime_key]
            self.jd = self.hdr[jd_key]
            self.savart, self.filter_color = self._find_filter(self.hdr[filter_key])
            self.object = self.hdr[object_key]
        except KeyError:
            print('Some key has not been found in header')

    def _get_header(self):
        try:
            hdr = fits.getheader(self.image_dir)
        except OSError:
            print('Header has not been found.')
        return hdr

    @property
    def data(self):
        return fits.getdata(self.image_dir)


class Pair:
	def __init__(self, p1, p2, p3, p4):
		self.p1 = 0
		self.p2 = 0
		self.p3 = 0
		self.p4 = 0


class ConfigurationSection:
    def __init__(self, name):
        self.name = name
        self.keys = []

    def set(self, key, value):
        setattr(self, key, value)
        self.keys.append(key)

    def get(self, key):
        if key not in self.keys:
            return None
        return getattr(self, key)

    def has_key(self, key):
        return key in self.keys

    def __repr__(self):
        return '{0}: {1}'.format(self.name,
                                 [(x, self.get(x)) for x in self.keys])


class Configuration:
    def __init__(self, path, fields):
        if not os.path.exists(path):
            raise ValueError('Configuration file has not been found.')

        self.parser = ConfigParser()
        self.parser.read(path)
        self.sections = []
        self.fields = fields

        for section in self.parser.sections():
            tmp = ConfigurationSection(section)

            for key, target_type in self.fields:
                tmp.set(key, self._get_configuration_value(section, key,
                                                           target_type))

            self.sections.append(tmp)

    def _get_configuration_value(self, section, key, data_type):
        if not self.parser.has_section(section):
            return None

        if not self.parser.has_option(section, key):
            return None

        if data_type is bool:
            return self.parser.getboolean(section, key)

        if data_type is int:
            return self.parser.getint(section, key)

        if data_type is float:
            return self.parser.getfloat(section, key)

        return self.parser.get(section, key)

    def get_section(self, name):
        for section in self.sections:
            if section.name == name:
                return section

        return None

    def has_param(self, section_name, key):
        for section in self.sections:
            if section.name == section_name:
                return key in section.has_key(key)
        return False

    def add(self, section_name, key, value):
        for section in self.sections:
            if section.name == section_name:
                section.set(key, value)

    def get(self, section_name, key):
        for section in self.sections:
            if section.name == section_name:
                section.get(key)

    def __repr__(self):
        return 'Configuration: {0}'.format([x for x in self.sections])


