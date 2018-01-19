from src.pipeline.base import PipelineBase
from src.utils.model import Configuration, Image

import os
import numpy as np
from glob import glob
from astropy.io import fits
from collections import namedtuple, OrderedDict
from photutils import source_properties, properties_table
from photutils import CircularAperture, EllipticalAperture
from photutils import CircularAnnulus, EllipticalAnnulus
from photutils import detect_sources


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
            ('sigma_clipped_sigma', float),
            ('sigma_clipped_iters', int),
            ('calculate_center', int),
            ('calculate_aperture', int),
            ('r_mask', int),
            ('r_aperture_multi', int),
            ('r_annulus_in', int),
            ('r_annulus_out', int),
            ('r_aperture', int)])

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


    def _make_apertures(self, image_stars_coordinates, data_shape):

        apertures = []
        masks = [self._create_mask(
            star_coo, data_shape) for star_coo in image_stars_coordinates]

        for i, mask in enumerate(masks):
            mean, median, std = sigma_clipped_stats(
                image.data, mask=mask,
                sigma=self.config_section.get('sigma_clipped_sigma'),
                iters=self.config_section.get('sigma_clipped_iters'))

            if self.config_section.get('calculate_center') or 
            self.config_section.get('calculate_aperture'):

                props = _make_props(np.copy(image.data), median,
                 mask, image_stars_coordinates[i])

            if self.config_section.get('calculate_center'):
                position = (props['xcentroid'],
                            props['ycentroid'])
            else:
                position = image_stars_coordinates[i]

            if self.config_section.get('calculate_aperture'):
                a = props['semimajor_axis_sigma'] * self.config_section.get('r_aperture_multi')
                b = props['semiminor_axis_sigma'] * self.config_section.get('r_aperture_multi')
                theta = props['orientation']

                aperture = EllipticalAperture(position, a=a, b=b, theta=theta)
                annulus = EllipticalAnnulus(position,
                                            a_in=a+self.config_section.get('r_annulus_in'),
                                            a_out=a+self.config_section.get('r_annulus_out'),
                                            b_out=b+self.config_section.get('r_annulus_out'),
                                            theta=theta)
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
        mask = x * x + y * y <= cfg.r_mask * cfg.r_mask
        mask_arr = np.full(data_shape, True, dtype=bool)
        mask_arr[mask] = False

        return mask_arr


    def _make_props(data, median, mask, star):

            #FIXME all variables to config 
            sigma = 2.0 * gaussian_fwhm_to_sigma    # FWHM = 2.
            kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
            kernel.normalize()
            data[mask] = 0
            segm = detect_sources(data, median*3, npixels=5,
                                  filter_kernel=kernel)
            props = properties_table(
                source_properties(data-np.uint64(median), segm),
                columns=['id', 'xcentroid', 'ycentroid', 'source_sum',
                         'semimajor_axis_sigma', 'semiminor_axis_sigma',
                         'orientation'])

            if len(props) > 1:
                props = findNearestObject(star, props, data.shape)
                return props
            else:
                return props[0]

    def make_phot(self):

        for image in self.images_list:
            apertures = self._make_apertures(self.stars_coordinates[image.savart],
                                             image.shape,)
            out_table = []

        # for aperture in apertures:
        #     rawflux_table = aperture_photometry(data, aperture[0])
        #     bkgflux_table = aperture_photometry(data, aperture[1])
        #     phot_table = hstack([rawflux_table, bkgflux_table],
        #                         table_names=['raw', 'bkg'])
        #     bkg_mean = phot_table['aperture_sum_bkg'] / aperture[1].area()
        #     bkg_sum = bkg_mean * aperture[0].area()
        #     final_sum = phot_table['aperture_sum_raw'] - bkg_sum
        #     phot_table['residual_aperture_sum'] = final_sum

        #     phot_table.add_column(
        #         Column(name='residual_aperture_err_sum',
        #                data=calcPhotErr(hdr, aperture,
        #                                 phot_table, bkgflux_table)))

        #     phot_table['xcenter_raw'].shape = 1
        #     phot_table['ycenter_raw'].shape = 1
        #     phot_table['xcenter_bkg'].shape = 1
        #     phot_table['ycenter_bkg'].shape = 1
        #     out_table.append(phot_table)

        # out_table = vstack(out_table)

        # if plot:
        #     makePlot(data, apertures, hdr['FILENAME'])

        # return out_table







def makeProps(data, median, mask, star):

        #FIXME all variables to config 
        sigma = 2.0 * gaussian_fwhm_to_sigma    # FWHM = 2.
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        kernel.normalize()
        data[mask] = 0
        segm = detect_sources(data, median*3, npixels=5,
                              filter_kernel=kernel)
        props = properties_table(
            source_properties(data-np.uint64(median), segm),
            columns=['id', 'xcentroid', 'ycentroid', 'source_sum',
                     'semimajor_axis_sigma', 'semiminor_axis_sigma',
                     'orientation'])

        if len(props) > 1:
            props = findNearestObject(star, props, data.shape)
            return props
        else:
            return props[0]


def findNearestObject(star, props, data_shape):
    min_dist = max(data_shape)
    min_dist_id = 0
    x, y = star

    for prop in props:
        dist = math.sqrt((x - prop['xcentroid']) ** 2 +
                         (y - prop['ycentroid']) ** 2)

        if dist < min_dist:
            min_dist = dist
            min_dist_id = prop['id']

    return props[min_dist_id-1]


def makeApertures(data, stars):



def makePhot


def calcPhotErr(hdr, aperture, phot_table, bkgflux_table):

    try:
        effective_gain = float(hdr[cfg.gain_key])
    except KeyError:
        effective_gain = 1.0

    err = math.sqrt(
        (phot_table['residual_aperture_sum'] / effective_gain) +
        (aperture[0].area() * aperture[2] ** 2) +
        ((aperture[0].area() ** 2 + aperture[2] ** 2) /
         (bkgflux_table['aperture_sum'] * aperture[0].area())))

    return [err]


def makePlot(data, apertures, im_name):

    plt.imshow(data, cmap='Greys', origin='lower',
               norm=LogNorm())
    for aperture in apertures:
        aperture[0].plot(linewidth=0.3, color='#d62728')
        aperture[1].plot(fill=False, linewidth=0.3, color='k')
        area = aperture[0]  # for plot mask area
        area.a, area.b = cfg.r_mask, cfg.r_mask
        area.plot(linewidth=0.2, color='k', ls=':')

    plt.savefig(im_name+'.png', dpi=300)
    plt.clf()


def saveTable(out_table, output_name):

    out_table = vstack(out_table)

    out_table.add_column(
        Column(name='COUNTS', data=out_table[cfg.output_flux]), index=0)
    out_table.add_column(
        Column(name='COUNTS_ERR', data=out_table['residual_aperture_err_sum']),
        index=1)


    out_table.write(output_name+'.csv', format='ascii', delimiter=',')



def loadImages():

    images = sorted(glob.glob(
        os.path.join(os.path.curdir, cfg.images)))
    # images = [name for name in images if len(name.split('.')) < 3]

    return images


def loadStars():

    try:
        stars = np.loadtxt(
            os.path.join(os.path.curdir, cfg.stars_file), dtype='f')
        stars = stars[np.argsort(stars[:, 1])]  # sort by y-axis
    except IOError:
        print('NO REGION FILE!')
        sys.exit()

    return stars
