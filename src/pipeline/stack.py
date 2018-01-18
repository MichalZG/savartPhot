from src.pipeline.base import PipelineBase
from src.utils.model import Configuration, Image

import ccdproc
from ccdproc import CCDData, Combiner
import os
from glob import glob
from astropy.io import fits
from astropy.time import Time


class Stack(PipelineBase):
	def __init__(self, config, work_path, stack_size, output_directory, log_file_name='stack'):
		super(Stack, self).__init__(log_file_name, None)
		self.config = Configuration(config, [
			('savarts_to_stack', str),
			('pattern', str),
			('datetime_key', str),
			('jd_key', str),
			('filter_key', str),
			('object_key', str)])
		self.config_section = self.config.get_section('stack')

		if not self.config_section:
			raise ValueError('Configuration file is not correct.')

		self.work_path = work_path
		self.stack_size = stack_size
		self.output_directory = output_directory
		self.images_list = self._create_images_list(self.work_path)
		self._create_directory(self.output_directory)


	def _create_directory(self, directory):
		try:
			os.makedirs(directory)
			self.info('Directory "{}" has been created.'.format(directory))
		except OSError:
			self.warning(
				'Directory "{}" could not be created.'.format(directory))


	def _create_images_list(self, work_dir):
		if not os.path.exists(work_dir):
			self.error('Directory {} has not been found'.format(work_dir))
			raise ValueError('Directory has not been found')
		
		images_dirs_list = sorted(
			glob(
				os.path.join(work_dir, self.config_section.get('pattern'))))
		
		if not images_dirs_list:
			raise ValueError('Empty images dirs list')

		images_list = [] 

		for image_dir in images_dirs_list:
			image = Image(image_dir, self.config_section.get('datetime_key'),
			 						 self.config_section.get('jd_key'),
			  						 self.config_section.get('filter_key'),
			  						 self.config_section.get('object_key'))
			images_list.append(image)

		self.info('Images list has been created')
		self.info('Images list length is {}'.format(len(images_list)))

		return images_list


	def _create_stack(self, images_list, stack_name):
		
		CCD_data_table = [CCDData(im.data, unit='adu') for im in images_list]
		combiner = Combiner(CCD_data_table)
		median = combiner.median_combine()
	
		master_hdr = self._create_stack_hdr(images_list,
		 self.config_section.get('datetime_key'),
		 self.config_section.get('jd_key'))

		self._save_stack(median, stack_name, master_hdr)


	def _calculate_mid_time(self, start_im, end_im):

		start_datetime = Time(start_im.datetime, format='isot')
		end_datetime = Time(end_im.datetime, format='isot')
		start_jd = start_im.jd
		end_jd = end_im.jd

		mid_datetime = start_datetime + (end_datetime - start_datetime) / 2.
		mid_jd = start_jd + (end_jd - start_jd) / 2.

		return mid_datetime.value, mid_jd

		
	def _create_stack_hdr(self, images_list, datetime_key, jd_key):
		
		stack_hdr = fits.Header()

		start_im = images_list[0]
		end_im = images_list[-1]

		stack_hdr['DATESTAR'] =  start_im.datetime
		stack_hdr['DATEEND'] = end_im.datetime
		stack_hdr['JDSTART'] = start_im.jd
		stack_hdr['JDEND'] = end_im.jd
		mid_datetime, mid_jd = self._calculate_mid_time(start_im, end_im)

		stack_hdr['DATEMID'] = mid_datetime
		stack_hdr['JDMID'] = mid_jd

		return stack_hdr


	def _create_stack_lists(self, names_of_savarts):

		stack_list = dict((name, []) for name in names_of_savarts)

		for image in self.images_list:
			if image.savart in stack_list:
				stack_list[image.savart].append(image)

		return stack_list


	def _save_stack(self, stack_arr, stack_name, master_hdr):
		CCDData.write(stack_arr, os.path.join(self.output_directory, stack_name),
			hdu_mask=None, hdu_uncertainty=None, clobber=True)
		f = fits.open(os.path.join(self.output_directory, stack_name), mode='update')
		f[0].header = master_hdr
		f.flush()


	def process(self):

		if not self.config_section.get('savarts_to_stack'):
			self.error('No savarts to stack in configuration file')
			raise ValueError('No savarts to stack')

		stack_lists = self._create_stack_lists(
			self.config_section.get('savarts_to_stack').split(','))

		for savart_name, stack_list in stack_lists.items():
			stack_name = savart_name + '.fits'
			self._create_stack(stack_list, stack_name)
