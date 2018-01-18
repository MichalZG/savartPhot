from src.pipeline.base import PipelineBase
from src.utils.model import Configuration

import os
import numpy as np
from glob import glob
from astropy.io import fits

class Savartphot(PipelineBase):
    def __init__(self, config, coordinates_file, log_file_name='savartphot'):
        super(Savartphot, self).__init__(log_file_name, None)
        self.config = (Configuration(config, [
            (output_file_ext, str)]))

        self.config_section = self.config.get_section('savartphot')

        if not self.config_section:
            raise ValueError('Configuration file is not correct.')

        self.coordinates_file = coordinates_file
        self.coordinates = self._load_stars_coordinates()


    def _load_stars_coordinates():
        
        if not os.path.exists(self.coordinates_file):
            self.error('Coordinates file {} has not been found'.format(
                self.coordinates_file))
            raise ValueError('Coordinates file has not been found')
        try:
            stars = np.loadtxt(
                os.path.join(os.path.curdir, cfg.stars_file), dtype='f')
            stars = stars[np.argsort(stars[:, 1])]  # sort by y-axis
        except IOError:
            print('NO REGION FILE!')
            sys.exit()

        return stars




def createMask(star, data_shape):
    y, x = np.ogrid[-star[1]-1:data_shape[0]-star[1]-1,
                    -star[0]-1:data_shape[1]-star[0]-1]
    mask = x * x + y * y <= cfg.r_mask * cfg.r_mask
    mask_arr = np.full(data_shape, True, dtype=bool)
    mask_arr[mask] = False

    return mask_arr


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
    apertures = []
    masks = [createMask(star, data.shape) for star in stars]

    for i, mask in enumerate(masks):
        mean, median, std = sigma_clipped_stats(data, mask=mask,
                                                sigma=1.0, iters=5)

        if cfg.calc_center or cfg.calc_aperture:

            props = makeProps(np.copy(data), median, mask, stars[i])

        if cfg.calc_center:

            position = (props['xcentroid'],
                        props['ycentroid'])
        else:
            position = stars[i]

        if cfg.calc_aperture:
            a = props['semimajor_axis_sigma'] * cfg.r_multi_ap
            b = props['semiminor_axis_sigma'] * cfg.r_multi_ap
            theta = props['orientation']
            aperture = EllipticalAperture(position, a=a, b=b, theta=theta)
            annulus = EllipticalAnnulus(position,
                                        a_in=a+cfg.r_ann_in,
                                        a_out=a+cfg.r_ann_out,
                                        b_out=b+cfg.r_ann_out,
                                        theta=theta)
            print(('i:{}, a:{}, b:{}, pos:{},{}').format(i+1, a, b, *position))
        else:
            aperture = CircularAperture(position, r=cfg.r_ap)
            annulus = CircularAnnulus(position,
                                      r_in=cfg.r_ap+cfg.r_ann_in,
                                      r_out=cfg.r_ap+cfg.r_ann_out)

        apertures.append([aperture, annulus, std])

    return apertures


def makePhot(data, hdr, stars, plot=False):

    apertures = makeApertures(data, stars)
    out_table = []

    for aperture in apertures:
        rawflux_table = aperture_photometry(data, aperture[0])
        bkgflux_table = aperture_photometry(data, aperture[1])
        phot_table = hstack([rawflux_table, bkgflux_table],
                            table_names=['raw', 'bkg'])
        bkg_mean = phot_table['aperture_sum_bkg'] / aperture[1].area()
        bkg_sum = bkg_mean * aperture[0].area()
        final_sum = phot_table['aperture_sum_raw'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum

        phot_table.add_column(
            Column(name='residual_aperture_err_sum',
                   data=calcPhotErr(hdr, aperture,
                                    phot_table, bkgflux_table)))

        phot_table['xcenter_raw'].shape = 1
        phot_table['ycenter_raw'].shape = 1
        phot_table['xcenter_bkg'].shape = 1
        phot_table['ycenter_bkg'].shape = 1
        out_table.append(phot_table)

    out_table = vstack(out_table)

    if plot:
        makePlot(data, apertures, hdr['FILENAME'])

    return out_table


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
