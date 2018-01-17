
from glob import glob
from src.utils.model import Configuration, Star
from astropy.io import fits
import os


class List_objects:
	def __init__(self, config):
		self.config = Configuration(config, [
			('pattern', str),
			('name_key', str)])
		self.config_section = self.config.get_section('list_object')
		self.stars_list = {}

		if not self.config_section:
			raise ValueError('Configuration file is not correct.')


	def create_files_list(self, work_dir):
		return sorted(glob(os.path.join(work_dir,
		 self.config_section.get('pattern'))))


	def create_stars_list(self, work_dir):
		files_list = self.create_files_list(work_dir)

		for file in files_list:
			try:
				hdr = fits.getheader(file)
			except OSError:
				if "unknown" in self.stars_list:
					self.stars_list["unknown"].files.append(file)
				else:
					self.stars_list["unknown"] = Star("unknown")
					self.stars_list["unknown"].files.append(file)

			try:
				name = hdr[self.config_section.get('name_key')]

				if name in self.stars_list:
					self.stars_list[name].files.append(file)
				else:
					self.stars_list[name] = Star(name)
					self.stars_list[name].files.append(file)

			except KeyError:
				if "unknown" in self.stars_list:
					self.stars_list["unknown"].files.append(file)
				else:
					self.stars_list["unknown"] = Star("unknown")
					self.stars_list["unknown"].files.append(file)

		return self.stars_list


	def show_objects(self, work_dir):

		self.create_stars_list(work_dir)

		for _, star in self.stars_list.items():
			print("{0}: {1}".format(star.name, len(star.files)))


