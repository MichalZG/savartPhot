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

		self.info('Processing stack {} finished'.format(stack_name))
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

		stack_hdr['DATETIME'] = mid_datetime
		stack_hdr['JD'] = mid_jd

		stack_hdr['IMNUM'] = len(images_list)
		stack_hdr['FILTER'] = '-'.join([start_im.savart, start_im.filter_color])

		return stack_hdr


	def _check_stack_list(self, stack_lists):

		if len(stack_lists) == 0:
			self.error('No files to stack in savart stack lists')
			raise ValueError('No files to stack in savart stack lists')
		elif len(stack_lists) == 1:
			self.warning('Only one savart list to stack')
		elif len(stack_lists) > 1:
			it = iter(stack_lists.values())
			ref_len = len(next(it))
			if not all(len(l) == ref_len for l in it):
				self.warning('Savarts stacks lists are not equal')

		for savart_name, savart_list in stack_lists.items():
			if len(savart_list) % self.stack_size != 0:
				self.warning(
					'Lenght ({}) of {} savart stack list is not a multiple of {}'.format(
					len(savart_list), savart_name, self.stack_size))
				self.warning('The remaining files will be skipped')


	def _create_stack_lists(self, names_of_savarts):

		stack_lists = dict((name, []) for name in names_of_savarts)

		for image in self.images_list:
			if image.savart in stack_lists:
				stack_lists[image.savart].append(image)

		self._check_stack_list(stack_lists)

		for	savart_name, stack_list in stack_lists.items():

			if self.stack_size > 0:
				if len(stack_list) % self.stack_size == 0:
					stack_lists[savart_name] = [
					stack_list[i:i + self.stack_size] for i in range(
						0, len(stack_list), self.stack_size)]
				else:
					stack_lists[savart_name] = [
					stack_list[i:i + self.stack_size] for i in range(
						0, len(stack_list), self.stack_size)][:-1]

		self.info('Stack lists have been created')
		return stack_lists


	def _save_stack(self, stack_arr, stack_name, master_hdr):
		CCDData.write(stack_arr, os.path.join(self.output_directory, stack_name),
			hdu_mask=None, hdu_uncertainty=None, clobber=True)
		f = fits.open(os.path.join(self.output_directory, stack_name), mode='update')
		f[0].header = master_hdr
		f.flush()
		self.info('Saving stack {} finished'.format(stack_name))


	def process(self):

		if not self.config_section.get('savarts_to_stack'):
			self.error('No savarts to stack in configuration file')
			raise ValueError('No savarts to stack')

		stack_lists = self._create_stack_lists(
			self.config_section.get('savarts_to_stack').split(','))

		for savart_name, stack_list in stack_lists.items():
			for i, chunk in enumerate(stack_list):
				stack_name = '{}_{}.fits'.format(savart_name, i+1)
				self.info('Processing stack {} started'.format(stack_name))
				self._create_stack(chunk, stack_name)
				