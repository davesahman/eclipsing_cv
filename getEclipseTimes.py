from hipercam.hlog import Hlog
from astropy import coordinates as coord
from astropy.time import Time
from astropy.convolution import Box1DKernel, convolve
from scipy.signal import medfilt
import copy


def tcorrect(tseries, star, observatory, type='B'):
    """
    Correct for light travel time.

    Arguments:
    ----------
    tseries: hipercam.hlog.Tseries
        Time series object

    star: astropy.coordinate.SkyCoord
        Location of star on Sky

    observatory: string
        Observatory name. See coord.EarthLocation.get_site_names() for list

    type: string (default=B)
        Heliocentric (H) or Barcentric (B)

    Returns
    -------
    tseries_corr : hipercam.hlog.Tseries
        Time series object with corrected time axis
    """
    ts = copy.deepcopy(tseries)
    times = Time(tseries.t, format='mjd', scale='utc',
                 location=coord.EarthLocation.of_site(observatory))
    if type == 'B':
        corr = times.light_travel_time(star)
        corr = times.tdb + corr
    else:
        corr = times.light_travel_time(star, 'heliocentric')
        corr = times.utc + corr
    ts.t = corr.mjd
    return ts


def smooth_derivative(tseries, med_half_width, box_half_width):
    """
    Calculate a smoothed version of the lightcurve derivative

    First smooth lightcurve wiht a median filter, then smooth
    numerical derivative with boxcar convolution.

    Parameters
    -----------
    tseries: hipercam.hlog.Tseries
        Time series object

    med_half_width: int
        Half-width of median filter

    box_half_width: int
        Half-width of boxcar filter

    Returns
    --------
    x, dy: np.ndarray
        Locations and values of smoothed derivative
    """
    x = tseries.t.copy()
    y = tseries.y.copy()
    yf = medfilt(y, 2 * med_half_width + 1)
    deriv = (yf[1:] - yf[:-1]) / (x[1:] - x[:-1])
    locs = 0.5 * (x[1:] + x[:-1])
    kernel = Box1DKernel(2 * box_half_width + 1)
    return locs, convolve(deriv, kernel)


def get_tseries(logfile, ccdnam, ap_targ, ap_comp):
    log = Hlog.from_ulog(logfile)
    return log.tseries(ccdnam, ap_targ) / log.tseries(ccdnam, ap_comp)
