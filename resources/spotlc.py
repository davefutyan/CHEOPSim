import os
import numpy as np
import scipy.stats
import scipy.interpolate
# from sklearn.neighbors import KernelDensity

# Load data
path_to_datafiles = os.environ['CHEOPS_SW']+'/resources'

class ParameterSampler(object):
    def __init__(self):
        # Spectral type bins:
        # Teff < 4000: M
        # 4000 <= Teff < 5200: K
        # 5200 <= Teff < 6000: G
        # Teff >= 6000: F
        sptype_bins = [0, 4000, 5200, 6000, 10000]

        self.data = {}
        for rotper in ('10', '20'):
            decaytime_day, teff, amplitude_mag = np.loadtxt(os.path.join(
                path_to_datafiles, '{}d_data.dat'.format(rotper)), skiprows=1,
                                                            unpack=True)

            indexes = np.searchsorted(sptype_bins, teff)

            for i, sp in enumerate(['M', 'K', 'G', 'F']):
                cond = indexes == i+1
                dt, amp = decaytime_day[cond], amplitude_mag[cond]
                self.data[rotper+'_'+sp] = np.array([amp, dt])

        self.keys = self.data.keys()

class PickSampler(ParameterSampler):

    def __init__(self):
        super(PickSampler, self).__init__()
        
    def draw(self, rotper, sptype):
        if rotper < 15:
            prot = '10'
        else:
            prot = '20'
        sp = sptype[0]
        
        amp, dt = self.data[prot+'_'+sp]
        ind = np.random.choice(range(len(amp)))
        return np.array([amp[ind], dt[ind]])

class KdeSampler(ParameterSampler):

    def __init__(self):
        super(KdeSampler, self).__init__()

        raise NotImplementedError('KDE sampler not implemented')
        self.kernels = {}
        for key in self.data:
            self.kernels[key] = kde_sklearn(self.data[key].T)
            
    def draw(self, rotper, sptype):
        if rotper < 15:
            prot = '10'
        else:
            prot = '20'
        sp = sptype[0]
        return self.kernels[prot+'_'+sp].sample()[0]

def quasiper_kernel(alpha, dx):
    """
    The quasiperiodic kernel function. The difference matrix
    can be an arbitrarily shaped numpy array so make sure that you
    use functions like ``numpy.exp`` for exponentiation.
    
    :param alpha: ``(4,)`` The parameter vector ``(amplitude, decay time,
    period, structure param)``.
    
    :param dx: ``numpy.array`` The difference matrix. This can be
        a numpy array with arbitrary shape.
    
    :returns K: The kernel matrix (should be the same shape as the
        input ``dx``). 
    
    """
    return alpha[0]**2 * np.exp(-0.5 * dx**2 / alpha[1]**2 - 2.0 *
                                np.sin((np.pi * dx / alpha[2]))**2 /
                                alpha[3]**2)

def sample_gp(x, alpha, kernelfn=quasiper_kernel, size=1):
    """
    :param array x: position array
    :param int n: number of steps in lightcurve
    """
    dx = x[:, None] - x[None, :]
    
    K = kernelfn(alpha, dx)

    return x, np.random.multivariate_normal(np.zeros_like(x), K, size)
    

def kde_sklearn(x, bandwidth=0.2, **kwargs):
    """
    Kernel Density Estimation with Scikit-learn
    """
    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(x)
    return kde_skl
    # score_samples() returns the log-likelihood of the samples
    # log_pdf = kde_skl.score_samples(x_grid[:, np.newaxis])
    #return np.exp(log_pdf)


def draw_parameters(prot, sptype, method='pick'):
    if method=='pick':
        return picksampler.draw(prot, sptype)
    elif method=='kde':
        return kdesampler.draw(prot, sptype)
    
def spot_lc(prot, sptype, ttotal_minutes=2880, dt=1., ncurves=1, maxpoints=5000,
            forcedt=False, save=True, outtemplate='lc'):
    """
    Function that produces the spot light curve.

    :param float prot: stellar rotational period in days
    :param str sptype: stellar spectral type. Only first character is
    considered.
    :param float ttotal_minutes: total run duration in minutes.
    :param float dt: time sampling in minutes.
    :param int maxpoints: maximum number of points admitted.
    """

    # Prepare time array (sample at most 100 times per period)
    if not forcedt:
        timesample_minutes = float(np.max((dt, prot*24*0.6)))
    else:
        timesample_minutes = float(dt)
        
    # Check that total number of points do not exceed 2000
    npoints = np.floor(ttotal_minutes/timesample_minutes)+1
    if npoints > maxpoints:
        raise ValueError('Simulation too long. Try reducing the length of '
                         'the simulation by a factor of {:.2f}.'
                         ''.format(float(npoints)/ maxpoints))

    t = np.linspace(0.0, ttotal_minutes, npoints)

    # Prepare parameter vector
    # Randomly set structure parameter between 0.5 and 1.0
    period_minutes = prot * 24 * 60
    structure = np.random.random() * 0.5 + 0.5

    # some function of Spectral type has to be implemented
    amplitude_mag, decaytime_days = draw_parameters(prot, sptype,
                                                    method='pick')
    # Convert amplitude in mag to amplitude in relative flux
    amplitude_relflux = 1 - 10**(-0.4*amplitude_mag)

    print('Randomly drawn parameters are:')
    print('Amplitude: {:.2f} mmag.'.format(amplitude_mag*1e3))
    print('Amplitude: {:.2e} rel. flux.'.format(amplitude_relflux))
    print('Decay time: {:.2f} days.'.format(decaytime_days))
    print('Structure parameter: {:.2f}.'.format(structure))
    
    alpha = [amplitude_relflux, decaytime_days*24*60, period_minutes,
             structure]

    t, f = sample_gp(t, alpha, size=ncurves)

    # Check if interpolation is needed (if the sampling is worse than
    # requested).
    if timesample_minutes > dt:
        print('Resampling curve to requested resolution.')
        # Interpolate
        time = np.linspace(0.0, ttotal_minutes, np.floor(ttotal_minutes/dt))
        flux = np.empty((ncurves, len(time)))
        for i in range(ncurves):
            flux[i] = scipy.interpolate.interp1d(t, f[i])(time)
    else:
        time, flux = t, f
            
    if save:
        # Find i0
        i0 = 0
        while os.path.exists(outtemplate+'{:05d}_00001'.format(i0)):
            i0+=1
            
        for i in range(ncurves):
            fout = open(outtemplate+'{:05d}_{:05d}'.format(i0, i+1), 'w')
            for j in range(len(time)):
                fout.write('{:.12f}\t{:.12f}\n'.format(time[j],
                                                       flux[i][j]+1))
            fout.close()
        
    return time, flux

picksampler = PickSampler()
#kdesampler = KdeSampler()
