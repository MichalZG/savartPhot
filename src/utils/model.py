from configparser import ConfigParser
import os
from math import sqrt, atan, pi
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

	@property
	def shape(self):
		return self.data.shape


class Coordinate:
	def __init__(self, name, x, y):
		self.name = name
		self.x = x
		self.y = y


class SavartCounts:
	def __init__(self, name, jd, counts_tab, counts_error_tab):
		self.name = name
		self.jd = jd
		self.counts1 = counts_tab[0]
		self.error1 = counts_error_tab[0]
		self.counts2 = counts_tab[1]
		self.error2 = counts_error_tab[1]

	@property
	def stokes(self):
		return (self.counts2 - self.counts1) / (self.counts2 + self.counts1)

	@property
	def stokes_error(self):
		return sqrt(
			4*(self.counts2**2*self.error1**2 + self.counts1**2*self.error2**2)/
			(self.counts2+self.counts1)**4)


class SavartsPair:
	def __init__(self, savart1, savart2):
		self.savart1 = savart1
		self.savart2 = savart2

	@property
	def pd(self):
		return sqrt(self.savart1.stokes**2 + self.savart2.stokes**2) * 100

	@property
	def pd_error(self):
		return sqrt(
			(self.savart1.stokes**2 * self.savart1.stokes_error**2 +
			 self.savart2.stokes**2 * self.savart2.stokes_error**2) /
			(self.savart1.stokes**2 + self.savart2.stokes**2)) 

	@property
	def pa(self):
		return 0.5*atan(self.savart2.stokes/self.savart1.stokes) * 180./pi

	@property
	def pa_error(self):
		return 0.5*sqrt(
			(self.savart1.stokes**2 * self.savart2.stokes_error**2 +
			 self.savart2.stokes**2 * self.savart1.stokes_error**2) /
			(self.savart1.stokes**2 + self.savart2.stokes**2) ** 2) 

	@property
	def jd(self):
		return self.savart1.jd + (self.savart2.jd - self.savart1.jd) / 2.


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


