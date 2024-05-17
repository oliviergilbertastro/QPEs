# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Photutils is an Astropy affiliated package to provide tools for
detecting and performing photometry of astronomical sources.  It also
has tools for background estimation, ePSF building, PSF matching,
centroiding, and morphological measurements.
"""

import warnings

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
from ._astropy_init import *  # noqa: F401, F403

from . import aperture
from . import background
from . import detection
from . import psf
from . import segmentation

# deprecations
__depr__ = {}

__depr__[aperture] = ('BoundingBox', 'CircularMaskMixin',
                      'CircularAperture', 'CircularAnnulus',
                      'SkyCircularAperture', 'SkyCircularAnnulus', 'Aperture',
                      'SkyAperture', 'PixelAperture', 'EllipticalMaskMixin',
                      'EllipticalAperture', 'EllipticalAnnulus',
                      'SkyEllipticalAperture', 'SkyEllipticalAnnulus',
                      'ApertureMask', 'aperture_photometry',
                      'RectangularMaskMixin', 'RectangularAperture',
                      'RectangularAnnulus', 'SkyRectangularAperture',
                      'SkyRectangularAnnulus', 'ApertureStats')

__depr__[background] = ('Background2D', 'BackgroundBase', 'BackgroundRMSBase',
                        'MeanBackground', 'MedianBackground',
                        'ModeEstimatorBackground', 'MMMBackground',
                        'SExtractorBackground', 'BiweightLocationBackground',
                        'StdBackgroundRMS', 'MADStdBackgroundRMS',
                        'BiweightScaleBackgroundRMS', 'BkgZoomInterpolator',
                        'BkgIDWInterpolator')

__depr__[detection] = ('StarFinderBase', 'DAOStarFinder', 'IRAFStarFinder',
                       'find_peaks', 'StarFinder')

__depr__[psf] = ('EPSFFitter', 'EPSFBuilder', 'EPSFStar', 'EPSFStars',
                 'LinkedEPSFStar', 'extract_stars', 'DAOGroup', 'DBSCANGroup',
                 'GroupStarsBase', 'NonNormalizable', 'FittableImageModel',
                 'EPSFModel', 'GriddedPSFModel', 'IntegratedGaussianPRF',
                 'PRFAdapter', 'BasicPSFPhotometry',
                 'IterativelySubtractedPSFPhotometry', 'DAOPhotPSFPhotometry',
                 'prepare_psf_model',
                 'get_grouped_psf_model', 'subtract_psf',
                 'resize_psf', 'create_matching_kernel',
                 'SplitCosineBellWindow', 'HanningWindow', 'TukeyWindow',
                 'CosineBellWindow', 'TopHatWindow')

__depr__[segmentation] = ('SourceCatalog', 'SegmentationImage', 'Segment',
                          'deblend_sources', 'detect_threshold',
                          'detect_sources', 'SourceFinder',
                          'make_2dgaussian_kernel')

__depr_mesg__ = ('`photutils.{attr}` is a deprecated alias for '
                 '`{module}.{attr}` and will be removed in the future. '
                 'Instead, please use `from {module} import {attr}` to '
                 'silence this warning.')

__depr_attrs__ = {}
for k, vals in __depr__.items():
    for val in vals:
        __depr_attrs__[val] = (getattr(k, val),
                               __depr_mesg__.format(module=k.__name__,
                                                    attr=val))
del k, val, vals

from photutils.detection.core import detect_threshold
from photutils.segmentation.detect import detect_sources
import numpy as np

def make_source_mask(data, nsigma, npixels, mask=None, mask_value=None,
                     filter_fwhm=None, filter_size=3, filter_kernel=None,
                     sigclip_sigma=3.0, sigclip_iters=5, dilate_size=11):
    """
    Make a source mask using source segmentation and binary dilation.

    Parameters
    ----------
    data : array_like
        The 2D array of the image.

    nsigma : float
        The number of standard deviations per pixel above the
        ``background`` for which to consider a pixel as possibly being
        part of a source.

    npixels : int
        The number of connected pixels, each greater than ``threshold``,
        that an object must have to be detected.  ``npixels`` must be a
        positive integer.

    mask : array_like, bool, optional
        A boolean mask with the same shape as ``data``, where a `True`
        value indicates the corresponding element of ``data`` is masked.
        Masked pixels are ignored when computing the image background
        statistics.

    mask_value : float, optional
        Deprecated.
        An image data value (e.g., ``0.0``) that is ignored when
        computing the image background statistics.  ``mask_value`` will
        be ignored if ``mask`` is input.

    filter_fwhm : float, optional
        The full-width at half-maximum (FWHM) of the Gaussian kernel to
        filter the image before thresholding.  ``filter_fwhm`` and
        ``filter_size`` are ignored if ``filter_kernel`` is defined.

    filter_size : float, optional
        The size of the square Gaussian kernel image.  Used only if
        ``filter_fwhm`` is defined.  ``filter_fwhm`` and ``filter_size``
        are ignored if ``filter_kernel`` is defined.

    filter_kernel : array-like (2D) or `~astropy.convolution.Kernel2D`, optional
        The 2D array of the kernel used to filter the image before
        thresholding.  Filtering the image will smooth the noise and
        maximize detectability of objects with a shape similar to the
        kernel.  ``filter_kernel`` overrides ``filter_fwhm`` and
        ``filter_size``.

    sigclip_sigma : float, optional
        The number of standard deviations to use as the clipping limit
        when calculating the image background statistics.

    sigclip_iters : int, optional
       The number of iterations to perform sigma clipping, or `None` to
       clip until convergence is achieved (i.e., continue until the last
       iteration clips nothing) when calculating the image background
       statistics.

    dilate_size : int, optional
        The size of the square array used to dilate the segmentation
        image.

    Returns
    -------
    mask : 2D bool `~numpy.ndarray`
        A 2D boolean image containing the source mask.
    """

    from scipy import ndimage

    threshold = detect_threshold(data, nsigma, background=None, error=None,
                                 mask=mask, mask_value=None,
                                 sigclip_sigma=sigclip_sigma,
                                 sigclip_iters=sigclip_iters)

    kernel = None
    if filter_kernel is not None:
        kernel = filter_kernel
    if filter_fwhm is not None:
        kernel_sigma = filter_fwhm * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(kernel_sigma, x_size=filter_size,
                                  y_size=filter_size)
    if kernel is not None:
        kernel.normalize()

    segm = detect_sources(data, threshold, npixels)
    if segm is None:
        return np.zeros(data.shape, dtype=bool)

    selem = np.ones((dilate_size, dilate_size))
    return ndimage.binary_dilation(segm.data.astype(bool), selem)


def __getattr__(attr):
    if attr in __depr_attrs__:
        obj, message = __depr_attrs__[attr]
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        return obj
    raise AttributeError(f'module {__name__!r} has no attribute {attr!r}')


# Set the bibtex entry to the article referenced in CITATION.rst.
def _get_bibtex():
    import os
    citation_file = os.path.join(os.path.dirname(__file__), 'CITATION.rst')

    with open(citation_file) as citation:
        refs = citation.read().split('@software')[1:]
        if len(refs) == 0:
            return ''
        bibtexreference = f'@software{refs[0]}'
    return bibtexreference


__citation__ = __bibtex__ = _get_bibtex()

del _get_bibtex
