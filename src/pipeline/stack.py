from src.pipeline.base import PipelineBase
from src.utils.model import Configuration

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
			('jd_key', str)])
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
		
		images_list = sorted(glob(os.path.join(work_dir, self.config_section.get('pattern'))))
		
		if not images_list:
			raise ValueError('Empty images list')
		self.info('Images list has been created')
		self.info('Images list length is {}'.format(len(images_list)))

		return images_list


	def _save_stack(self, stack_arr, stack_name, master_hdr):
		CCDData.write(stack_arr, os.path.join(self.output_directory, stack_name),
			hdu_mask=None, hdu_uncertainty=None, clobber=True)
		f = fits.open(os.path.join(self.output_directory, stack_name), mode='update')
		f[0].header = master_hdr
		f.flush()


	def _create_stack(self, images_list, stack_name):
		
		CCD_data_table = [CCDData.read(im, unit='adu') for im in images_list]
		combiner = Combiner(CCD_data_table, dtype='float32')
		median = combiner.median_combine()
	
		master_hdr = self._create_stack_hdr(images_list,
		 self.config_section.get('datetime_key'),
		 self.config_section.get('jd_key'))

		self._save_stack(median, stack_name, master_hdr)


	def _calculate_mid_time(self, start_hdr, end_hdr, datetime_key, jd_key):

		start_datetime = Time(start_hdr[datetime_key], format='isot')
		end_datetime = Time(end_hdr[datetime_key], format='isot')
		start_jd = start_hdr[jd_key]
		end_jd = end_hdr[jd_key]

		mid_datetime = start_datetime + (end_datetime - start_datetime) / 2.
		mid_jd = start_jd + (end_jd - start_jd) / 2.

		return mid_datetime.value, mid_jd

		
	def _create_stack_hdr(self, images_list, datetime_key, jd_key):
		
		stack_hdr = fits.Header()

		start_hdr = fits.getheader(images_list[0])
		end_hdr = fits.getheader(images_list[-1])

		stack_hdr['DATESTART'] =  start_hdr[datetime_key]
		stack_hdr['DATEEND'] = end_hdr[datetime_key]
		stack_hdr['JDSTART'] = start_hdr[jd_key]
		stack_hdr['JDEND'] = end_hdr[jd_key]
		mid_datetime, mid_jd = self._calculate_mid_time(start_hdr, end_hdr,
													    datetime_key, jd_key)

		stack_hdr['DATEMID'] = mid_datetime
		stack_hdr['JDMID'] = mid_jd

		return stack_hdr


	def process(self):
		pass

