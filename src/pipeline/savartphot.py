from src.pipeline.base import PipelineBase
from src.utils.model import Configuration, Image

import os
import math
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from collections import namedtuple, OrderedDict
from photutils import source_properties, properties_table
from photutils import CircularAperture, EllipticalAperture
from photutils import (CircularAnnulus, EllipticalAnnulus,
                       detect_sources, aperture_photometry)
from astropy.stats import gaussian_fwhm_to_sigma, sigma_clipped_stats
from astropy.table import hstack, vstack, Column
from astropy.visualization import astropy_mpl_style
from astropy.visualization import (ImageNormalize, MinMaxInterval,
                                   LinearStretch, ZScaleInterval)

plt.style.use(astropy_mpl_style)

class Savartphot(PipelineBase):
    def __init__(self, config, work_path, coordinates_file, output_directory,
     log_file_name='savartphot'):
        super(Savartphot, self).__init__(log_file_name, None)
        self.config = Configuration(config, [
            ('output_file_ext', str),
            ('pattern', str),
            ('datetime_key', str),
            ('jd_key', str),
            ('filter_key', str),
            ('object_key', str),
            ('gain_key', str),
            ('gain_default', float),
            ('sigma_clipped_sigma', float),
            ('sigma_clipped_iters', int),
            ('calculate_center', int),
            ('calculate_aperture', int),
            ('r_mask', int),
            ('r_aperture_multi', float),
            ('r_annulus_in', int),
            ('r_annulus_out', int),
            ('r_aperture', float),
            ('fwhm', float),
            ('kernel_x', float),
            ('kernel_y', float),
            ('detect_threshold', float),
            ('npixels', int),
            ('plot', int),
            ('aperture_linewidth', float),
            ('aperture_in_color', str),
            ('aperture_out_color', str),
            ('area_linewidth', float),
            ('area_color', str),
            ('area_linestyle', str),
            ('flux_type', str),
            ('flux_error_type', str)])

        self.config_section = self.config.get_section('savartphot')

        if not self.config_section:
            raise ValueError('Configuration file is not correct.')

        self.coordinates_file = coordinates_file
        self.work_path = work_path
        self.stars_coordinates = self._load_stars_coordinates()
        self.output_directory = output_directory
        self.images_list = self._create_images_list(self.work_path)
        self._create_directory(self.output_directory)
        self.measurements = {}


    def _create_directory(self, directory):

        try:
            os.makedirs(directory, exist_ok=True)
            self.info('Directory "{}" has been created.'.format(directory))
        except OSError:
            self.warning(
                'Directory "{}" could not be created.'.format(directory))


    def _load_stars_coordinates(self):

        if not os.path.exists(self.coordinates_file):
            self.error('Coordinates file {} has not been found'.format(
                self.coordinates_file))
            raise ValueError('Coordinates file has not been found')
            
        try:
            coordinates_table = np.loadtxt(
                self.coordinates_file, dtype=[('Savart', '<U5'),
                                              ('X1', '<f8'),
                                              ('Y1', '<f8'),
                                              ('X2', '<f8'),
                                              ('Y2', '<f8')])
        except TypeError:
            self.error('Problem with coordinates file')

        stars_coordinates = dict(
            (c[0], [[c[1], c[2]], [c[3], c[4]]]) for c in coordinates_table)

        self.info('Coordinates has been loaded.')

        return stars_coordinates


    def _create_images_list(self, work_dir):

        if not os.path.exists(work_dir):
            self.error('Directory {} has not been found'.format(work_dir))
            raise ValueError('Directory has not been found')
        
        images_dirs_list = sorted(
            glob(os.path.join(work_dir, self.config_section.get('pattern'))))
        
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


    def _make_apertures(self, image, image_stars_coordinates, data_shape):

        apertures = []
        masks = [self._create_mask(
            star_coo, data_shape) for star_coo in image_stars_coordinates]

        for i, mask in enumerate(masks):
            mean, median, std = sigma_clipped_stats(
                image.data, mask=mask,
                sigma=self.config_section.get('sigma_clipped_sigma'),
                iters=self.config_section.get('sigma_clipped_iters'))

            if (self.config_section.get('calculate_center') or 
            self.config_section.get('calculate_aperture')):

                star_properties = self._find_star_properties(np.copy(image.data), median,
                 mask, image_stars_coordinates[i])

            if self.config_section.get('calculate_center'):
                position = (star_properties['xcentroid'].value,
                            star_properties['ycentroid'].value)
            else:
                position = image_stars_coordinates[i]

            if self.config_section.get('calculate_aperture'):
                a = star_properties['semimajor_axis_sigma'] * self.config_section.get(
                    'r_aperture_multi')
                b = star_properties['semiminor_axis_sigma'] * self.config_section.get(
                    'r_aperture_multi')
                theta = star_properties['orientation']

                aperture = EllipticalAperture(position, a=a.value, b=b.value, theta=theta.value)
                annulus = EllipticalAnnulus(position,
                                            a_in=a.value+self.config_section.get(
                                                'r_annulus_in'),
                                            a_out=a.value+self.config_section.get(
                                                'r_annulus_out'),
                                            b_out=b.value+self.config_section.get(
                                                'r_annulus_out'),
                                            theta=theta.value)
                print(('i:{}, a:{}, b:{}, pos:{},{}').format(i+1, a, b, *position))
            else:
                aperture = CircularAperture(position, r=self.config_section.get('r_aperture'))
                annulus = CircularAnnulus(position,
                                          r_in=self.config_section.get(
                                            'r_aperture')+self.config_section.get(
                                            'r_annulus_in'),
                                          r_out=self.config_section.get(
                                            'r_aperture')+self.config_section.get(
                                            'r_annulus_out'))

            apertures.append([aperture, annulus, std])

        return apertures


    def _create_mask(self, star_coo, data_shape):
        y, x = np.ogrid[-star_coo[1]-1:data_shape[0]-star_coo[1]-1,
                        -star_coo[0]-1:data_shape[1]-star_coo[0]-1]
        mask = x * x + y * y <= self.config_section.get('r_mask') ** 2
        mask_arr = np.full(data_shape, True, dtype=bool)
        mask_arr[mask] = False

        return mask_arr


    def _find_star_properties(self, data, median, mask, star_coo):

            sigma = self.config_section.get('fwhm') * gaussian_fwhm_to_sigma 
            kernel = Gaussian2DKernel(sigma,
             x_size=self.config_section.get('kernel_x'),
             y_size=self.config_section.get('kernel_y'))
            kernel.normalize()
            data[mask] = 0
            segm = detect_sources(data,
                median*self.config_section.get('detect_threshold'),
                npixels=self.config_section.get('npixels'),
                filter_kernel=kernel)
            properties = properties_table(
                source_properties(data-np.uint64(median), segm),
                columns=['id', 'xcentroid', 'ycentroid', 'source_sum',
                         'semimajor_axis_sigma', 'semiminor_axis_sigma',
                         'orientation'])

            if len(properties) > 1:
                properties = self._find_nearest_object(star_coo, properties, data.shape)
                return properties
            else:
                return properties[0]


    def _find_nearest_object(self, star_coo, properties, data_shape):
        min_dist = max(data_shape)
        min_dist_id = 0
        x, y = star_coo

        for prop in properties:
            dist = math.sqrt((x - prop['xcentroid'].value) ** 2 +
                             (y - prop['ycentroid'].value) ** 2)

            if dist < min_dist:
                min_dist = dist
                min_dist_id = prop['id']

        return properties[min_dist_id - 1]


    def _calc_phot_error(self, hdr, aperture, phot_table, bkgflux_table):

        try:
            effective_gain = float(hdr[self.config_section.get('gain_key')])
        except KeyError:
            effective_gain = self.config_section.get('gain_default')

        err = math.sqrt(
            (phot_table['residual_aperture_sum'] / effective_gain) +
            (aperture[0].area() * aperture[2] ** 2) +
            ((aperture[0].area() ** 2 + aperture[2] ** 2) /
             (bkgflux_table['aperture_sum'] * aperture[0].area())))

        return [err]


    def _make_plot(self, data, apertures, im_name):

        norm = ImageNormalize(data, interval=ZScaleInterval(), 
                      stretch=LinearStretch())

        plt.imshow(data, cmap='Greys', origin='lower',
                   norm=norm)
        for aperture in apertures:
            aperture[0].plot(
                linewidth=self.config_section.get('aperture_linewidth'),
                color=self.config_section.get('aperture_in_color'))
            aperture[1].plot(fill=False,
                linewidth=self.config_section.get('aperture_linewidth'),
                color=self.config_section.get('aperture_out_color'))
            area = aperture[0]  # for plot mask area
            area.a, area.b = [self.config_section.get('r_mask')] * 2
            area.plot(linewidth=self.config_section.get('area_linewidth'),
                color=self.config_section.get('area_color'),
                ls=self.config_section.get('area_linestyle'))

        plt.savefig(os.path.join(self.output_directory, im_name+'.png'), dpi=900)
        plt.clf()


    def _save_table(self, out_table, output_name):

        out_table = vstack(out_table)
        out_table.add_column(
            Column(name='COUNTS', data=out_table[
                self.config_section.get('flux_type')]),
             index=0)
        out_table.add_column(
            Column(name='COUNTS_ERR', data=out_table[
                self.config_section.get('flux_error_type')]),
            index=1)

        out_table.write(os.path.join(self.output_directory, output_name),
            format='ascii', delimiter=',')


    def make_phot(self):

        for image in self.images_list:
            apertures = self._make_apertures(image, self.stars_coordinates[image.savart],
                                             image.shape)
            out_table = []

            for aperture in apertures:
                rawflux_table = aperture_photometry(image.data, aperture[0])
                bkgflux_table = aperture_photometry(image.data, aperture[1])
                phot_table = hstack([rawflux_table, bkgflux_table],
                                    table_names=['raw', 'bkg'])
                bkg_mean = phot_table['aperture_sum_bkg'] / aperture[1].area()
                bkg_sum = bkg_mean * aperture[0].area()
                final_sum = phot_table['aperture_sum_raw'] - bkg_sum
                phot_table['residual_aperture_sum'] = final_sum

                phot_table.add_column(
                    Column(name='residual_aperture_err_sum',
                           data=self._calc_phot_error(image.hdr, aperture,
                            phot_table, bkgflux_table)))

                phot_table['xcenter_raw'].shape = 1
                phot_table['ycenter_raw'].shape = 1
                phot_table['xcenter_bkg'].shape = 1
                phot_table['ycenter_bkg'].shape = 1
                out_table.append(phot_table)

            out_table = vstack(out_table)

            if self.config_section.get('plot'):
                self._make_plot(image.data, apertures, image.name)

            self._save_table(out_table, image.name + '.csv')