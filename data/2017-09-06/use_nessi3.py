import matplotlib.pyplot as plt
import numpy as np
# from nessi.tester import load_data
from nessi import integrator as nss
import matplotlib.colors as colors
import matplotlib.colors as mcolors
from astropy.io import fits as f
from scipy.interpolate import interp1d
import copy
import sunpy
import cocopy as cp
# from ISPy.io 
import solarnet # very use full to get time and to get spectral and coordinate positions
import sunpy.map
import astropy.units as u
from scipy.io import readsav as rs
from PIL import Image, ImageEnhance
from random import randint
import data_analysis as da
import datetime as dt
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.collections import PolyCollection
import os


class linestudier():
    '''
    linestudier() is a class of theoretical lineprofiles in a solar spectrum (quiet sun).
    It uses nessi to calculate starting of the disk centre (dc) profile and the centre to limb
    variations te overal spectrum.
    It has some functions to visualize an work with these profiles.

        param
        ------
        data: str
            filename of a .fits file for each wavelength that can be found here:
            https://cdsarc.cds.unistra.fr/ftp/J/A+A/671/A130/fits/ (This is Quiet sun data)
            (e.g. '7772_clv.fits')
        nr: int
            The 'nr' argument specifies the number of radial segments and azimuthal segments to use in the computation.
        map: Bool
            If a visual clearifying plot of the theoretical clv
            (=centre to limb variation) should be plotted
    '''
    def __init__(self, data_filename, atlas=None, nr=51, map=True, neglect_atlas=False):
        self.neglect_atlas = neglect_atlas
        if not neglect_atlas:
            if atlas is None:
                self.atlas = f.getdata(get_file_path_fits("E:/solar flares/data/2017-09-06/fits/solar_atlas_V1_405-1065.fits"))
            else:
                self.atlas = atlas
            atlas_w = np.arange(len(self.atlas)) * -0.003766534468 + 24700.0858041
            self.aw = 1e8 / atlas_w
        self.saas = nss.sun_as_a_star(nr=nr)

        self.theoretic_line(data_filename=data_filename, map=map) # load data for dc profile
        self.enroll_atlas() # load full disk spectrum atlas which will act as a control measure


    '''
    theoretic_line() takes data for the quiet sun from a website.
    This is then used to find a theoretical line profile

        param
        ------
        data: str
            filename of a .fits file for each wavelength that can be found here:
            https://cdsarc.cds.unistra.fr/ftp/J/A+A/671/A130/fits/ (This is Quiet sun data)
            (e.g. '7772_clv.fits')
        map: Bool
            If a visual clearifying plot of the theoretical clv
            (=centre to limb variation) should be plotted
    '''
    def theoretic_line(self, data_filename, map=True):
        self.cont_point=None
        if data_filename[-5:] == '.fits':
            self.has = f.open(get_file_path_fits(data_filename))

            self.sst_wav = self.has[1].data
            self.sst_mu = self.has[2].data
            self.sst_int = self.has[4].data.T
            self._extracted_from_theoretic_line_9()
        elif data_filename[-4:] == '.npy':
            self.has2 = np.load(get_file_path_line_data(data_filename), allow_pickle=True)
            self.sst_wav = self.has2[0]
            self.sst_mu  = self.has2[1]
            self.sst_int = self.has2[3].T
            self._extracted_from_theoretic_line_9()
        else:
            raise TypeError('Only .fits and .npy supported as nessi clv.')

        self.n = np.shape(self.sst_wav)

        # visual aid
        if map:
            self.map_expected_clv()

    # TODO Rename this here and in `theoretic_line`
    def _extracted_from_theoretic_line_9(self):
        self.sst_dc = self.sst_int[0]
        # CLV profiles should be normalized in such a way that at mu=1 it is is 1 all over, and the rest is lower.
        self.sst_clv = self.sst_int.copy()
        self.sst_clv /= self.sst_clv[0]
            


    def disgard_cont_point(self, spectrum, do=True):
        if do and (self.cont_point is not None):
            print("np.shape(spectrum)", np.shape(spectrum))
            if np.shape(spectrum) == self.n:
                return np.delete(spectrum, self.cont_point)
            else:
                raise("The point at continuum has already been disgarded")            
        return spectrum

    '''
    map_expected_clv() gives a visualisation of the dc-profile, the clv and the resulting
    sun-as-a-star (saas) profile.
    '''
    def map_expected_clv(self):

        # create_colormap
        cm = plt.cm.copper(np.linspace(0,1, 50))

        fig = plt.subplots(1, 3, figsize=(25,4))
        plt.subplot(131)
        plt.plot(self.sst_wav, self.sst_dc)
        plt.title('DC profile')
        plt.xlabel(r'Wavelength [$\rm\AA$]')
        plt.ylabel('Intensity')

        ax =plt.subplot(132)
        # Reverse the color map
        ax.set_prop_cycle('color', list(cm)[::-1])
        ax.plot(self.sst_wav,self.sst_clv.T)
        plt.title('CLV ')
        plt.xlabel(r'Wavelength [$\rm\AA$]')
        plt.ylabel('Intensity')
        norm = mcolors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.copper, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label(r'$\mu$')

        self.saas.update_clv(self.sst_mu,self.sst_wav,self.sst_clv,self.sst_wav,self.sst_dc)
        self.saas.update_vrot(0.,0.)
        self.saas_profile = self.saas.get_integration()

        plt.subplot(133)
        plt.plot(self.sst_wav, self.saas_profile)
        plt.title('SAAS Profile')
        plt.xlabel(r'Wavelength [$\rm\AA$]')
        plt.ylabel('Intensity')
        plt.show()
        
    def save_saas_profile(self, line_name):
        np.save(f'saas_profiles/{line_name}.npy')

    '''
    enroll_atlas()
        We include an atlas to compare the data with (always good to check if all works well!)
        Cut the atlas near where the line is (Keep in mind that we use air wavelenths,
        and they use vacuum, so youll need to shift it)
    '''
    def enroll_atlas(self):
        if self.neglect_atlas:
            return
        self.lw = self.sst_wav[0]
        self.rw = self.sst_wav[-1]

        if self.lw > self.aw[-1] or self.rw < self.aw[0]:
            print(f"The line of sst and atlas do not coincde! SST: min={self.lw}, max={self.rw}")
            print(f"and atlas: min={self.aw[0]}, max={self.aw[-1]}. The advise is taking a new consistent atlas.")
            raise WrongLineError("The line of sst and atlas do not coincde.")

        self.llw = np.where(self.aw > self.lw-3)[0][0]
        self.lrw = np.where(self.aw > self.rw+3)[0][0]

        self.fd = self.atlas[self.llw:self.lrw]
        self.fdw = self.aw[self.llw:self.lrw]

        self.ff = interp1d(self.fdw, self.fd, kind='linear', fill_value="extrapolate")
        # Use this function to compute the new values
        self.fdd = self.ff(self.sst_wav)



    '''
    saas_profile() calculates and plots the saas profile of a spectral line from its dc-profile
    and clv, integrating over the solar disk.

        param
        --------
        with_atlas: Bool
            If True the plots are accompagned by the atlass and the convolved atlas
            Default: True
        show_all: Bool
            Do we show all plots or only one. Default: False

        TODO:
            making an fitting algoritm for the scaling factors for the atlas
            making it usefull for all wavelengts.
    '''
    def saas_profile_atlas_check(self, show_all=False,initial_values_fit=None):
        if self.neglect_atlas:
            return
        test_si = self.saas_profile

        #test_si = test_si/np.min(test_si)*np.min(test_si)
        test_si = test_si/test_si[0]
        self.sst_dc = self.sst_dc/self.sst_dc[0]
        self.fdd = self.fdd/self.fdd[0]

        self.saas_profile = test_si


        import ISPy.spec.crisp as c
        dw = 0.07
        ntw = 59
        tw = (np.arange(ntw)-ntw//2)*dw
        # fpife = c.crisp(np.median(self.sst_wav))
        fc = c.crisp(np.median(self.sst_wav))
        tr = fc.dual_fpi(tw, erh = -0.022) #6301-2

        tr /= tr.sum()
        inst_prof = np.zeros((len(tr),2))
        inst_prof[:,0] = tw+np.abs(tw.min())
        inst_prof[:,1] = tr

        import ISPy.spec.calib as cb

        aa = cb.convolve_atlas(self.fdw,self.fd,inst_prof)
        self.aa2=aa
        self.fit_atlas_to_nessi(quality=show_all, initial_values=initial_values_fit)

        if show_all:
            print(len(self.fdw))

            plt.plot(self.fdw,self.fd, label="Atlas")
            plt.plot(self.fdw, aa, label="Atlas Convolved")
            plt.legend()
            plt.show()


        self.ff = interp1d(self.fdw+self.theta[0], self.aa2*self.theta[2]+self.theta[1], kind='linear', fill_value="extrapolate")
        # Use this function to compute the new values
        self.fdd = self.ff(self.sst_wav)


        fig = plt.subplots(1, 2, figsize=(15,7))
        plt.subplot(121)
        plt.plot(self.fdw + self.theta[0] ,self.fd*self.theta[2]+self.theta[1], color='black', linestyle='--', alpha=0.5, label='Atlas (Reiners et al. 2015)')
        plt.plot(self.fdw + self.theta[0],self.aa2*self.theta[2]+self.theta[1], color='black', lw=2, label='Atlas Convolved')
        plt.plot(self.sst_wav, self.sst_dc, color='blue', label='Disk Center')
        plt.plot(self.sst_wav, test_si, color='red', label='Full CLV')
        plt.xlabel(r'Wavelength [$\rm\AA$]')
        plt.ylabel('Intensity')
        plt.legend()
        # plt.xlim()

        self.atlas_convolved = [self.fdw + self.theta[0],self.aa2*self.theta[2]+self.theta[1]]

        plt.subplot(122)
        plt.plot(self.sst_wav, test_si/self.fdd, label='CLV', color='red')
        plt.plot(self.sst_wav, self.sst_dc/self.fdd, label='DC', color='blue')
        plt.xlabel(r'Wavelength [$\rm\AA$]')
        plt.ylabel('Offset [%]')
        plt.legend()

        # plt.ylim(0.93,1.025)
        plt.show()

    def fit_atlas_to_nessi(self, minimum=True, quality=True, initial_values=None):
        if self.neglect_atlas:
            return

        if minimum:
            # minimum determination if the lines are wells
            xatlas = np.where(self.fd == np.min(self.fd))
            # print(xatlas[0][0])
            x_si = np.where(self.saas_profile == np.min(self.saas_profile))
            # print(x_si[0][0])


            # print("the minima are at ", self.fdw[xatlas][0], '(atlas) and', self.sst_wav[x_si][0], "(nessi). Therefor they are ",
            #     self.fdw[xatlas][0]-self.sst_wav[x_si][0],"appart.")
            rv = - self.fdw[xatlas][0] + self.sst_wav[x_si][0]  #1.8431060192415316
        else:
            rv=0


        l = len(self.sst_wav)
        data = [self.sst_wav,  self.saas_profile ,np.zeros(l)+0.001,np.zeros(l)+0.001]

        multiplier = 1/0.85

        # theta = [horizontale translatie, verticale translatie, verticale schaalfactor]
        model_atlas = lambda theta : interp1d(self.fdw + theta[0], self.aa2*theta[2] + theta[1], kind='linear', fill_value="extrapolate")


        # self.fdw + rv,aa2*multiplier, color='black', lw=2, label='Atlas Convolved'
        # self.sst_wav, test_si, color='red', label='Full CLV'
        if initial_values is None:
            initial_values=np.array([rv, 0, multiplier])

        mini = da.optimalisatie(data, model=model_atlas, beginwaarden=initial_values, fout_model=None, plot=False)
        self.theta = mini['x']
        da.plot_fit(data, model=model_atlas, theta0=mini['x'], titel="Fitting atlas to nessi ",labelx=" $wavelength [\AA]$",
                    labely=" $Intensity$  [arbitrary units]", figname=None , error=False)
        if quality:
            print(mini)
            da.kwaliteit_fit(data, mini)
   

    def set_fov(self, X, Y, boundary, sr=959.63, reduction=False):
        boundary = convert_boundary_to_nan(boundary[:-1, :-1])
        X /= sr
        Y /= sr
        r = (X[1:,1:] + X[0:-1,0:-1])**2 / 4 + (Y[1:,1:] - Y[0:-1,0:-1])**2 / 4
        r = np.where(r >1, np.nan, 0)
        boundary = boundary * r
        self.fov = [X, Y, boundary]
        self.saas.update_clv(self.sst_mu,self.sst_wav,self.sst_clv,self.sst_wav,self.sst_dc)
        self.saas.update_vrot(0.,0.)
        self.saas_profile = self.saas.get_integration()
        dx = (X[1:,1:] - X[0:-1,0:-1])
        dy = (Y[1:,1:] - Y[0:-1,0:-1])
        self.fov_areafactor = np.nansum(dx*dy+boundary) / np.pi
        if reduction > 1 :
            X, Y, boundary = reduce_(X, Y, boundary, reduction)
        fov_spectra = np.array([boundary for _ in range(len(self.sst_wav))])
        # areafactor = 1 / np.pi * (X[-1,-1]-X[0,0]) * (Y[-1,-1]-Y[0,0])
        self.spectr_fov = -self.saas.get_diff_spectra_fov(X,Y,fov_spectra) / self.fov_areafactor
        plt.plot(self.sst_wav, self.saas_profile, label="SAAS profile NESSI")
        plt.plot(self.sst_wav, self.spectr_fov, label="FOV profile NESSI")
        plt.legend()
        plt.show()
        
    def set_quiet_sun(self, xlim, ylim, show=False, xyswitch=True, sr=959.63, reduction=False):
        if not hasattr(self, 'fov'):
            raise BufferError('self.fov is not yet defined. Please first run self.set_fov().')
        self.quiet_sun = [xlim, ylim]
        a, b = xlim
        c, d = ylim
        if xyswitch:
            a , b, c, d = c, d, a, b
        X = self.fov[0][a:b, c:d]
        Y = self.fov[1][a:b, c:d]
        boundary = self.fov[2][a:b-1, c:d-1]
        dx = (X[1:,1:] - X[0:-1,0:-1])
        dy = (Y[1:,1:] - Y[0:-1,0:-1])
        areafactor = np.nansum(dx*dy+boundary) /np.pi
        x0, y0 = np.nanmean(X), np.nanmean(Y)
        print(f"The quiet sun is pick at x {x0*sr, y0*sr} with average\
              mu {np.cos(np.arcsin((x0**2 +y0**2)**0.5))}")
        if reduction > 1 :
            X, Y, boundary = reduce_(X, Y, boundary, reduction)
        qs_spectra = np.array([boundary for _ in range(len(self.sst_wav))])
        if show:
            print(X,Y,qs_spectra, areafactor)
        self.spectr_qs = -self.saas.get_diff_spectra_fov(X,Y,qs_spectra) / areafactor
        
        plt.plot(self.sst_wav, self.saas_profile, label="SAAS profile NESSI")
        plt.plot(self.sst_wav, self.spectr_qs, label="QS profile NESSI")
        plt.plot(self.sst_wav, self.spectr_fov, label="FOV profile NESSI")
        plt.legend()
        plt.show()
        
    def adapted_fov_spectr(self):
        if not hasattr(self, "theta_nessi_to_quiet_sun"):
            raise ValueError("linestudier has no attribure theta_nessi_to_quiet_sun. First gauge to quiet sun!")
        theta = self.theta_nessi_to_quiet_sun
        return interp1d(self.sst_wav + theta[0], theta[2] * self.spectr_fov + theta[1]
                                    , kind='linear', fill_value="extrapolate")
        
    def adapted_qs_spectr(self):
        if not hasattr(self, "theta_nessi_to_quiet_sun"):
            raise ValueError("linestudier has no attribure theta_nessi_to_quiet_sun. First gauge to quiet sun!")
        theta = self.theta_nessi_to_quiet_sun
        return interp1d(self.sst_wav + theta[0], theta[2] * self.spectr_qs + theta[1]
                                    , kind='linear', fill_value="extrapolate")

def convert_boundary_to_nan(boundary):
    return np.where(boundary == 0, np.nan, 0)
     
         
def reduce_(X, Y, boundary, reduction=5):
    n = reduction
    X, Y, boundary = X[::n, ::n], Y[::n, ::n], boundary[::n,::n]
    assert(X.shape == Y.shape)
    if X.shape[0] == boundary.shape[0]:
        boundary = boundary[:-1]
    if X.shape[1] == boundary.shape[1]:
        boundary = boundary[:,:-1]
    print(X.shape, Y.shape, boundary.shape)
    return X, Y, boundary
    
        
def fix_mu_theor(theor_line, mu):
    theor_line.exact_mu = mu
    x = np.abs(theor_line.sst_mu-mu)
    index_mu = np.where(x == np.min(x))[0]
    theor_line.index_mu = index_mu
    theor_line.best_fit_clv = clv_fit(mu, theor_line)
    print('Mu also set to the theoretic nessi line.')


def clv_fit(mu, theor_line):
    return np.apply_along_axis(lambda arr: interp1d(theor_line.sst_mu, arr)(mu), axis=0, arr=theor_line.sst_clv)

def fit_qs_to_NESSI(theor_line, sst_data, initial_guess=None):
    import data_analysis as da
    # theta = [horizontale translatie, verticale translatie, verticale schaalfactor]
    # theta = [0.2, 0.3, 0.89]
    frame = sst_data.quiet_sun['frame']

    f_nessi_qs = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[1] * theor_line.spectr_qs 
                                    , kind='linear', fill_value="extrapolate")
    f_nessi_fov = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[1] * theor_line.spectr_fov 
                                    , kind='linear', fill_value="extrapolate")
    f_nessi_saas = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[1] * theor_line.saas_profile 
                                    , kind='linear', fill_value="extrapolate")
    g = len(sst_data._wavel)

    #To simulate a specific domain around the well we cam make the errors on the wings huge
    avs = sst_data.frame_integrated_spect(sst_data.quiet_sun['frame'], sst_data.quiet_sun['xlim'], sst_data.quiet_sun['ylim'], variation=True)
    stds = sst_data.var_spect
    print(sst_data._wavel, avs, f_nessi_qs([0,1])(sst_data._wavel))
    data = [sst_data._wavel,  avs ,stds,np.zeros(g)+0.01]
    if initial_guess is None:
        initial_guess = np.array([np.median(sst_data._wavel) - np.median(theor_line.sst_wav), np.average(avs)/np.average(theor_line.spectr_qs)])
    else:
        print(f'the given inital guess is {initial_guess}')
    print(initial_guess)

    mini = da.optimalisatie(data, model=f_nessi_qs, beginwaarden=initial_guess, fout_model=None, plot=False)
    print(mini)
    t0 = mini['x']
    theta = [t0[0], 0, t0[1]]
    sst_data.theta_nessi_to_quiet_sun = theta
    theor_line.theta_nessi_to_quiet_sun = theta
    da.plot_fit(
        data,
        model=f_nessi_qs,
        theta0=t0,
        titel="SST quiet sun patch vs nessi (fit)",
        labelx="wavelength $[\AA]$",
        labely=" Intensity []",
        error=False,
    )
    da.kwaliteit_fit(data, mini)

    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
    # fig.setitle("quiet sun determination, weiging constants")

    fig.suptitle('The gauge of the quiet sun spectrum', fontsize=16)

    sst_data.ccp_frame(frame,Show=False)

    ax[0].set_title(
        f"spectral lines of SST frame {str(frame)}  with normalized intensity"
    )
    sst_data.frame_integrated_spect(frame)
    ax[0].plot(sst_data._wavel, sst_data.av_spect, '--', label='sst FOV')
    ax[0].plot(sst_data._wavel, sst_data.quiet_spect, '--', label='sst quiet sun') #

    ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_saas(t0)(theor_line.sst_wav + theta[0]), label='NESSI saas')
    ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_qs(t0)(theor_line.sst_wav + theta[0]), label='NESSI QS')
    ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_fov(t0)(theor_line.sst_wav + theta[0]), label='NESSI FOV')
    ax[0].legend()
    ax[1].imshow(sst_data.current_ccp, origin='lower')
    xlim = sst_data.quiet_sun['xlim']
    ylim = sst_data.quiet_sun['ylim']
    color_wanted_patch = sst_data.quiet_sun['color']
    ax[1].plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], [ylim[0],ylim[0],ylim[1], ylim[1], ylim[0]], color=color_wanted_patch)
    ax[1].set_title(f"COCOplot of frame {str(frame)}")
    plt.show()

def fit_qs_to_NESSI_cont_point(theor_line, sst_data,  frame=0, xshift=None):
    xlim, ylim = theor_line.quiet_sun
    obs = sst_data.frame_integrated_spect(frame, xlim=xlim, ylim=ylim)
    theta = [(xshift if xshift is not None else sst_data._wavel[-1]-theor_line.sst_wav[-1]), 0,  obs[-1] /theor_line.spectr_qs[-1]]
    sst_data.theta_nessi_to_quiet_sun = theta
    print("the theta fit is ",theta)

    f_nessi_qs = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[2] * theor_line.spectr_qs + theta[1]
                                    , kind='linear', fill_value="extrapolate")
    f_nessi_fov = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[2] * theor_line.spectr_fov + theta[1]
                                    , kind='linear', fill_value="extrapolate")
    f_nessi_saas = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[2] * theor_line.saas_profile + theta[1]
                                    , kind='linear', fill_value="extrapolate")
    
    fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
    # fig.setitle("quiet sun determination, weiging constants")

    fig.suptitle('The gauge of the quiet sun spectrum', fontsize=16)

    sst_data.ccp_frame(frame,Show=False)

    ax[0].set_title(
        f"spectral lines of SST frame {str(frame)}  with normalized intensity"
    )
    sst_data.frame_integrated_spect(frame)
    ax[0].plot(sst_data._wavel, sst_data.av_spect, '--', label='sst FOV')
    ax[0].plot(sst_data._wavel, sst_data.quiet_spect, '--', label='sst quiet sun') #

    ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_saas(theta)(theor_line.sst_wav + theta[0]), label='NESSI saas')
    ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_qs(theta)(theor_line.sst_wav + theta[0]), label='NESSI QS')
    ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_fov(theta)(theor_line.sst_wav + theta[0]), label='NESSI FOV')
    ax[0].legend()
    ax[1].imshow(sst_data.current_ccp, origin='lower')
    color_wanted_patch = sst_data.quiet_sun['color']
    ax[1].plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], [ylim[0],ylim[0],ylim[1], ylim[1], ylim[0]], color=color_wanted_patch)
    ax[1].set_title(f"COCOplot of frame {str(frame)}")
    ax[0].set_xlim((np.min(sst_data._wavel), np.max(sst_data._wavel[:-1])))
    ax[0].set_ylim(0,10)
    plt.show()


























def intensity_match(I1, I2, error=0.001):
    G = (np.average(I1)+np.average(I2)) / 2
    D = np.max(np.abs(I1-I2))
    return D <= error * G
    
def guess_filters(n_wav):
    sd = n_wav/8
    return [[n_wav/6 - 0.5, sd], [3 * n_wav/6 - 0.5, sd], [5 * n_wav / 6 - 0.5, sd]]

def file_name_line_full_disk(name_of_line, full_path=None):
    filename = f'line_full_disk_{name_of_line}.npy'
    if full_path is not None:
        filename = os.path.join(full_path,filename)
    return filename

def get_file_path_line_data(filename, full_path=None):
    if filename[-4:] != ".npy":
        filename = f"{filename}.npy"
    filename = os.path.join('line_data',filename)
    if full_path is not None:
        filename = os.path.join(full_path,filename)
    return filename

def get_file_path_FOV(name_of_line, full_path=None):
    filename = f'FOV_spectrum_{name_of_line}.npy'
    filename = os.path.join('FOV', filename)
    if full_path is not None:
        filename = os.path.join(full_path,filename)
    return filename

def get_file_path_fits(name):
    if name[-5:] == '.fits':
        name = os.path.join('fits', name)
    else:
        raise TypeError('File should be .fits')
    return name


def official_start_flare(name, except_time=None):
    """Returns the indices over which to average to get a good quiet pattern to calculate the contrast profile from.

    Args:
        name (string): name of the line

    Returns:
        [t0, t1]: begin and end time of the quiet flare time
    """
    if "23a" in name:
        return "09:30:00"
    elif "23" in name:
        return "16:48:00"
    elif "15b" in name:
        return "11:55:00"
    elif "24a" in name:
        return "12:27:00"
    elif "22a" in name:
        return "08:00:00"
    elif "22" in name:
        return "07:35:00"
    elif "16" in name:
        return "11:28:00"
    elif "21a" in name:
        return "15:17:00"   
    elif "21" in name:
        return "09:46:00"
    for yy in ['19','13''9u', "17a","17","15a","15","14a","14"]:
        if yy in name:
            print(f'Not official start flare defined for {name} falling back to first of TIME[0], the execpt_time {except_time}.')
            return except_time
    raise NameError(f'WRONG NAME: the line {name} had no official start flare defined.')

class SST_data_from_multiple_fits_files():

    # timeframe_to_filename_fits is a function that sends the index of the frame to the associated filename
    # spectfilename is alike 'spectfile6563_93.idlsave'
    # timesfilename is alike "times6563_93_2017_09_06_11_55_47.idlsave"
    #
    def __init__(self, timeframe_to_filename_fits, spectfilename, number_of_frames, time, name_of_line,
                 thresh=[1e-10,2e-7], boundary_methode='search', alternative_datacube=None,
                 boundary_arguments=None, cont_point=None, with_stokes=False, with_time=False, shape=None):
        if alternative_datacube is not None:
            self.alternative_datacube = alternative_datacube
        self.with_stokes = with_stokes
        self.with_time = with_time
        self.ttff = timeframe_to_filename_fits
        self.filename = timeframe_to_filename_fits(0)
        print(self.filename)
        self.name_of_line = name_of_line
        self._number_of_frames = number_of_frames
        self.cont_point = None
        if shape is not None:
            self.shape = shape
        else:
            s = np.shape(self.datacube(0))
            print(f'{s=}')
            self.shape = (number_of_frames, 1 , s[-3], s[-2], s[-1]) if not with_stokes else (number_of_frames, '??' , s[0], s[1], s[2])
        print(f"the shape of the data is {self.shape}")
        
        
        # spectral positions 
        if spectfilename == 'use_solarnet':
            self._wavel = solarnet.get_wav(self.filename)
            print("the found wavelengts by solarnet where ", self._wavel)
        else:
            try:
                self._spectfile = rs(spectfilename)
                self._wavel=self._spectfile["spect_pos"]
            except:
                self._wavel=np.load(spectfilename)
        print(f"the wavelengths are {self._wavel = }")
        # Time file
        self._time = time
        assert(len(time) == self.shape[0])
        # cont_point is the index of the point in the continuum to exculde this in plots (will be mostly -1)
        if cont_point is not None:
            self.define_cont_point(cont_point)

        self._thresh=thresh
        # try:
        self.set_boundary_original(methode=boundary_methode,arguments=boundary_arguments)
        # except:
        #     print("IT WAS IMPOSSIBLE TO CONSTRUCT THE BOUNDARY AUTOMATICALLY!\nTRY TO DO THIS BY HAND.\n",
        #           "So try: \n\tself.set_boundary_original(methode=boundary_methode,arguments=boundary_arguments)",
        #           "\nWith other arguments.\n")
        print("The next thing to do is to initalise the filters. Use update_filters(self, MeanSd, form='normal')")
    
    def datacube(self, timeframe):
        if hasattr(self, "alternative_datacube"):
            return self.alternative_datacube(timeframe)
        if timeframe >= self._number_of_frames or timeframe < 0:
            raise IndexError(f"Timeframe {timeframe} is negative or exeeds numer of frames {self._number_of_frames}.")
        d = f.getdata(self.ttff(timeframe))
        if self.with_stokes:
            d = d[0]
        if self.with_time:
            d = d[0]
        return d
    
    # MeanSd: an array which gives mean and standard deviation for the RGB filters for the chosen form.
        # for a normal form for example: [[10,1.25], [6,1], [2,1.25]]
    # form: the form of filters. Normal distribution is Default.
    def update_filters(self, MeanSd=None, form="normal"):
        if MeanSd is None:
            MeanSd = guess_filters(len(self._wavel))
            print(f"The guessed filters are {MeanSd}")
        wavelengths = np.arange(self.shape[2])
        self._filt = cp.filter(wavelengths, form, MeanSd, plot=True)

    def define_cont_point(self, index):
        self.cont_point = index
        # line_lim wil give the xlim voor plotting a linesegment without the point in the continuum
        wav2 = self._wavel[:-1] if index==-1 else np.concatenate(self._wavel[:index], self._wavel[index+1:])
        min = np.min(wav2)
        max = np.max(wav2)
        D = max-min
        self.line_lim = [min-D/10 , max+D/10]
    
    def set_xlim_without_continuum(self, ax):
        ax.set_xlim(self.line_lim)

    def plot_filt(self, y,x):
        cp.filtplot(self.datacube(0)[:,y,x], self._filt)

    def disgard_cont_point(self, spectrum, do=True):
        if do and (self.cont_point is not None):

            print("np.shape(spectrum)", np.shape(spectrum), spectrum)
            if np.shape(spectrum)[0] == self.shape[2]:
                return np.delete(spectrum, self.cont_point)
            else:
                print(f"The point at continuum has already been disgarded. Got a {np.shape(spectrum)} array which should have been a {self.shape[0]}-array" )       
        return spectrum
    
    def set_mu(self, theor_line, number_of_last_frame=None, alternative_filename=None):
        filename = self.filename
        if alternative_filename is not None:
            filename = alternative_filename
        if number_of_last_frame is None:
            number_of_last_frame = self._number_of_frames - 1
        mu1 = give_mu_contourplot(filename, over=self, save=True, save_name=self.name_of_line)
        mu2 = give_mu_contourplot(filename, over=self,
                                    timeFrame=number_of_last_frame)
        self.mu = (mu1 + mu2) / 2
        print('We take average mu to be ', self.mu)
        fix_mu_theor(theor_line, self.mu)
        return self.mu

    def time_of_frame(self, frame):
        s = str(self._time[frame])
        if s[0] == 'b':
            s = s[2:]
        s = s.split('.')[0]
        tstr=("t="+s+" UT")
        self.current_time = tstr
        return tstr

    def ccp_frame(self, frame, Show=True):
        #coarse graining for memory reasons
        sx = 1 if (self.shape[3] < 2100) else self.shape[3] // 1000
        sy = 1 if (self.shape[4] < 2100) else self.shape[4] // 1000
        # Now view a datacube
        self.current_frame=frame
        cube = np.nan_to_num(self.datacube(frame)[:,::sx,::sy])
        cube[np.where(cube > self._thresh[1])] =self._thresh[1]
        cube[np.where(cube <self._thresh[0])] =self._thresh[0]
        tstr = self.time_of_frame(frame)
        if Show:
            print("COCOPLOT at ", tstr, "(frame number", frame,")")
        self.current_ccp = cp.plot(cube,self._filt, show=Show)
        # return self.current_ccp

    # frame is an integer, the frame number
    # pixels is an array listing [x,y] of the pixels of which the spectrum is to appear
    def interesting_pixels(self, frame, pixels, rand=False, numb=5):
        if rand:
            xmax=self.shape
            ymax=xmax[3]
            xmax = xmax[4]
            pixels = [[randint(0,xmax-1), randint(0, ymax-1)] for _ in range(numb)]
        else:
            numb = len(pixels)
        self.ccp_frame(frame,Show=False)
        colors=self.current_ccp

        fig, ax = plt.subplots(1)
        ax.set_title(f"spectral lines of some pixels: {self.time_of_frame(frame)}")
        if len(pixels)>0:
            corr = self.correction if hasattr(self, 'correction') else self._wavel*0
        for p in range(numb):
            ax.plot(
                self._wavel,
                self.datacube(frame)[:, pixels[p][1], pixels[p][0]] / self.scalar
                + corr,
                color=colors[pixels[p][1]][pixels[p][0]] / 255,
                label=f"pixel x={str(pixels[p][0])} y={str(pixels[p][1])}",
            )
        if hasattr(self, 'line_lim'):
            ax.set_xlim(self.line_lim)
        ax.set_xlabel(r"wavelength [$\rm\AA$]")
        ax.set_ylabel("intensity [units???]")
        if numb<6:
            ax.legend(fontsize=6)
        plt.show()

    def view_enhanced_frame(self, frame):
        # Now view a datacube
        self.ccp_frame(frame, Show=False)
        colplt=self.current_ccp
        img = Image.fromarray(colplt)
        new_img = ImageEnhance.Contrast(img).enhance(factor=2)
        print("Enchanced COCOPLOT at ", self.time_of_frame(frame), "(frame number", frame,")")
        plt.imshow(new_img, aspect='auto', origin='lower')
        plt.show()

    def calculate__FOV_spect_over_time(self, show_total_spectrum=False):
        if hasattr(self, 'FOV_spectrum'):
            if show_total_spectrum:
                self.plot_total_intensity()
            return
        filename = get_file_path_FOV(str(self.name_of_line))
        try :
            FOV_spectrum = np.load(filename)
        except FileNotFoundError:
            FOV_spectrum=self.frame_integrated_spect(0)[None]
        if np.shape(FOV_spectrum)[0] < self._number_of_frames:
            print(f'In total {self.shape[0]} frames.\nNow calculating frame:', )
            for frame in range(np.shape(FOV_spectrum)[0], self.shape[0]):
                s = ''.join('\r' for _ in range(len(str(frame))+1))
                print(s, end=str(frame))
                FOV_spectrum = np.concatenate((FOV_spectrum, self.frame_integrated_spect(frame)[None]) ) 
                np.save(filename, FOV_spectrum)

        self.FOV_spectrum = FOV_spectrum
        if show_total_spectrum:
            self.plot_total_intensity()

        
    def plot_total_intensity(self):
        Total_intensity = np.sum(self.FOV_spectrum, axis=1)
        # print(Total_intensity)
        frame_max = np.where(Total_intensity == np.max(Total_intensity))[0]
        print(f'The peak occurs at frame {frame_max} at time {[self.time_of_frame(i) for i in frame_max]}.')
        print(f"The SST data runs from {self.time_of_frame(0)} to {self.time_of_frame(-1)}.")

        fig, ax = plt.subplots()
        ax.plot(Total_intensity)
        ax.set_title("Spectral integrated intensity over time")
        ax.set_xlabel("number of frames")
        ax.set_ylabel('Intensity')
        plt.show



    def add_correction(self, theoretical_difference):
        if hasattr(self, 'correction'):
            print('the correction has already been done ')
            return self.FOV_spectrum, self.correction

        mean_difference_observation = np.mean(self.FOV_spectrum[:,-1]-self.FOV_spectrum[:,0])

        correction = -(mean_difference_observation - theoretical_difference) * (self._wavel - self._wavel[0])

        self.FOV_spectrum += correction
        self.correction = correction
        return self.FOV_spectrum, correction





    '''
    set_boundary_original() tries to make an boolean array for y,x which result False if the pixel is 'Outside'
    and True if the pixel actually contains data, so is 'inside'.

    param:
        methode: 'search', 'By_user'
            'search': this option let use_nessi try to calculate this itself.
            'By_user': Here the user gives the vertices of the quadrilateral that surrounds the boundary.

        arguments:
            'search': this option does not require any arguments
            'By_user': A list of tupples of the vertices of the quadrilateral that surrounds the boundary.

    It does it only for 1 color outside!!!
    TODO: go to two
    '''
    def set_boundary_original(self, methode = 'search', arguments=None):
        self.arguments = arguments
        if methode == 'search':
            self._search_boundary(arguments)
        elif methode == 'By_user':
            quadrilateral_vertices = arguments
            self.zeros=[]
            grid_shape = np.shape(self.datacube(0)[0,:,:])
            self.boundary = is_point_inside_quadrilateral_grid(grid_shape, quadrilateral_vertices)
            self._boundary_per_frame = False

        elif methode == 'No Boundary':
            self.zeros=[]
            print(f'defining no boundary for shape {self.shape }')
            self.boundary = np.ones(self.shape[-2:])
            print()
            self._boundary_per_frame = False

        elif methode == 'search_overlap':
            self._boundary_per_frame = True
            self._zero_per_frame = True
            print("Use arguments as dictionary with 'frames_of_overlap' to define a list of considerd frames to calculate overlap")
            f_o_o = arguments.get('frames_of_overlap', range(0,self._number_of_frames, self._number_of_frames//8))
            for i,frame in enumerate(f_o_o):
                self.zeros = self.calculate_zeros(frame=frame)
                B = self.calculate_boundary(frame=frame)
                self.boundary = B if i==0 else B * self.boundary
                print(f'calculated boundary at frame {frame}')
                plt.imshow(B)
                plt.show()
            print(f"calculated boundary over folowing frames {f_o_o}")
            self._boundary_per_frame = False
        else:
            raise ValueError(f'The provided Boundary methode {methode} is not implemented! only "search", "No Boundary", "By_user" and "search_overlap"')

        self.check_scalar_not_nan()


        self.plot_boundary()
        

    def _search_boundary(self, arguments):
        self.arguments=arguments
        self._boundary_per_frame = True
        self._zero_per_frame = 'zeros' not in arguments
        self.calculate_boundary()

    def check_scalar_not_nan(self):
        if hasattr(self, 'scalar'):
            if np.isnan(self.scalar):
                print('This is a problem. The self.scalar is nan.')
                # A scalar which will normalize the intensity
                self.scalar = 1
                self.scalar = np.average(self.frame_integrated_spect(0))
                if np.isnan(self.scalar):
                    print('The problem is not fixed by renormalization.\nMake sure no other constants are nan in the definition of the scalar')
        else:
            self.frame_integrated_spect(0)

    def calculate_zeros(self, zero="still to search", error=0.01, frame=0):
        if self.arguments is not None:
            if "zeros" in self.arguments:
                return self.arguments["zeros"]
        '''
        Here there are 2 posibilities for the 'outside' of te frame. Either these are denoted by NaN and then
        averaging is simpel by a numpy.mean() function wich only include the needed pixel

        If however by correcting for the stokes factors the 'outside' of the frame is turned into a spectrum as well than
        we need to find the zero/outside spectrum in order that we can leave the pixels with this spectum out. We will first
        suppose this collor is not different for any of the frames. So we find the color in frame = 0
        '''

        # OPTION 1 NAN


        if np.any(np.isnan(self.datacube(frame)[:,:,:])):
            return [np.nan]

        # OPTION 2 FAKE OUTSIDE DATA
        else:
            xmax=self.shape
            ymax=xmax[3]
            xmax = xmax[4]

            if zero=="still to search":
                # to find zeros. If not the param
                #
                # zero: gives an RGB value of the outside. any pixel that match it is excluded
                #
                # first we secure the "outside" of the frame so we can neglec it
                # most probable guess is on the outside lines in the middels. and check by looking to its neigbours.

                pixels_to_check = [[0,ymax//2],[xmax//2,0],[0+5,ymax//2],[xmax//2,0+5],[0,ymax//2+5],[xmax//2,0+5],
                                    [0,ymax//2-5],[xmax//2-5,0],[0+5,ymax//2+5],[xmax//2+5,0+5],[0+5,ymax//2-5],[xmax//2-1,0+5]]

                # as guess we take the first and check if they all match
                zero = self.datacube(frame)[:,ymax//2,0]
                possible = True

                for p in range(len(pixels_to_check)):
                    x,y=pixels_to_check[p]
                    # print("zero:", zero)
                    # print("vergelijk met:", self.datacube(frame)[:,y,x])

                    if not intensity_match(zero, self.datacube(frame)[:,y,x], error=error):
                        possible = False
                    # print("en ", possible, "bevonden")


                if not possible:
                    print("IMPORTANT MESSAGE\nThe middel of the outside lines is in the data.\n Give manualy the RGB value of the outside in.")
                    print("zero was:", zero)
                    print(f"The error for matching intensities is {error}")


            return [zero]

    def calculate_boundary(self, frame=0):
        
        error = self.arguments.get('error', 0.001) if self.arguments is not None else 0.001
        ind = self.arguments.get('index_of_zero', np.arange(len(self._wavel))) if self.arguments is not None else np.arange(len(self._wavel))
        self.zeros = self.calculate_zeros(error=error)
        assert(len(ind)>0)
        
        if not self._boundary_per_frame:
            return self.boundary
        if self._zero_per_frame:
            self.zeros = self.calculate_zeros(error=error, frame=frame)
        # print(f"calculation of the boundary and zero for frame {frame}")
        try:
            if np.any(np.isnan(self.zeros)):
                self.boundary = np.where(np.isnan(self.datacube(frame)[0,:,:]), 0, 1)
                return self.boundary
        except TypeError as e:
            print(f'{self.zeros = }')
            raise e

        shape=self.shape
        ymax=shape[3]
        xmax = shape[4]
        R0 = np.ones(self.datacube(frame).shape[1:])
        
        for zero in self.zeros:            
            R = np.ones(self.datacube(frame).shape[1:])
            for i in ind:
                    R *= np.where(np.abs(self.datacube(frame)[i,:,:]-zero[i])< error * np.average(self.datacube(frame)), 1, 0)
            R0 *= 1 - R    
        self.boundary = R0
        return R0


    def expect_pix(self, frame=0):
        y = self.shape[3]//2
        x = np.where(self.boundary[y,:])
        # print(x, x[0][0])
        return y, x[0][0]

    def second_boundary(self, plot=True, frame=0):
        '''
        This function has to be called if there is too little been excluded as a boundary.
        '''
        self._boundary_per_frame = False

        xmax=self.shape
        ymax=xmax[3]
        xmax = xmax[4]
        y, x= self.expect_pix(frame=frame)

        # as guess we take the first and check if they all match
        zero = self.datacube(frame)[:,y,x]
        # print(zero)
        # print(self.datacube[frame, 0,:,y,10:20])

        self.zeros.append(zero)
        # make now a weigthing array for averaging
        R = np.array([[(self.datacube(0)[:,i,j]==zero)[0]
                        for j in range(xmax)] for i in range(ymax)])
        R = np.where(R, 0, 1)

        self.boundary = (R + self.boundary)-1


        if plot:
            print(x, self.boundary[y, x-5:x+5])
            self.plot_boundary()

    def plot_boundary(self):
        print("Boundary\nBlue=Outside, yellow=Inside")
        plt.imshow(self.boundary,origin='lower')
        plt.title('Boundary')
        plt.show()
        print("If there is still boundary left to be excluded, call self.second_boundary()")

    def spectrum_in_posible_spectra(spectrum_given, posible_spectra):
        for i in range(len(posible_spectra)):
            if np.all(spectrum_given == posible_spectra[i]):
                return True
        return False

    # def set_boundary_frame(self, frame):
    #     #Option 1: NaN
    #     if np.any(np.isnan(self.datacube[frame, 0,:,:,:])):
    #         R = np.isnan(self.datacube[frame, 0,0,:,:])
    #         R = np.where(R, 0, 1)
    #         print("Original file: NaN are detected, exact boundary.")
    #     #Option 2: colored in boundaries
    #     else:
    #         xmax=np.shape(self.datacube)
    #         ymax=xmax[3]
    #         xmax = xmax[4]
    #         R = np.array([(self.datacube[frame,0,:,i,j]==self.zeros[0])[0] for i in range(ymax) for j in range(xmax)])
    #         for l in range(1,len(self.zeros),1):
    #             R += np.array([(self.datacube[frame,0,:,i,j]==self.zeros[l])[0] for i in range(ymax) for j in range(xmax)])
    #         R = np.where(R, 0, 1).reshape((ymax, xmax))

    #         self.boundary = R



    # Function gives the average/integrated spectrum. We normalize by setting first intensity to zero.
    # However there is a part of the frame which does not contain data so this has to be excluded.
    def frame_integrated_spect(self, frame, xlim=None, ylim=None, variation=False, without_cont_point=True):
        if self._boundary_per_frame:
            self.calculate_boundary(frame=frame)
        xmin = 0 if xlim is None else xlim[0]
        xmax = self.shape[4] if xlim is None else xlim[1]
        ymin = 0 if ylim is None else ylim[0]
        ymax = self.shape[3] if ylim is None else ylim[1]

        R = self.boundary
        corr = self.correction if hasattr(self, 'correction') else self._wavel*0
        if not hasattr(self, "scalar"):
            # A scalar which will normalize the intensity
            if np.any(np.isnan(self.zeros)):
                self.scalar =np.nanmean(self.datacube(0)[0,:,:])
            else:
                print(np.shape(self.datacube(0)[0,:,:]), np.shape(R))
                self.scalar = np.average([np.average(self.datacube(0)[i,:,:],weights=R) for i in range(len(self._wavel))])


        if np.any(np.isnan(self.zeros)):
            self.av_spect = self.disgard_cont_point(
                np.array([np.nanmean(self.datacube(frame)[i,ymin:ymax,xmin:xmax]) for i in range(len(self._wavel))]) / self.scalar + corr
                )
        else:
            self.av_spect = self.disgard_cont_point(
            np.array([np.average(self.datacube(frame)[i,ymin:ymax,xmin:xmax],weights=R[ymin:ymax,xmin:xmax]) for i in range(len(self._wavel))]) / self.scalar + corr
            )
        if variation:
            npxl = (ymax-ymin)*(xmax-xmin) #number of pixels
            if np.any(np.isnan(self.zeros)):
                self.var_spect = np.array([np.nanstd(self.datacube(frame)[i,ymin:ymax,xmin:xmax]/self.scalar) for i in range(len(self._wavel))])/npxl**0.5
            else:
                self.var_spect = np.array([np.std(self.datacube(frame)[i,ymin:ymax,xmin:xmax]/self.scalar)for i in range(len(self._wavel))])/npxl**0.5

        return self.av_spect




    # frame: frame number
    # zero: either "still to search" or an RGB array to excluded. (not important if .frame_integrated_spect has already runned.)
    # pixels : [[x,y], [x',y'],...] for values of pixels you wont to have plotted allong. Default [] leads to no extra lines
    def plot_integrate_spectr(self, frame, pixels=[]):
        self.frame_integrated_spect(frame)


        fig, ax = plt.subplots(1,1, figsize=(5,5), sharey=False)
        ax.set_title(
            f"spectral lines of SST frame {str(frame)} in relative intensity "
        )
        ax.plot(self._wavel, self.av_spect,color='red', label="sst averaged") # linewidth

        # to give relative intensities
        if len(pixels) > 0:
            corr = self.correction if hasattr(self, 'correction') else self._wavel*0
        for p in range(len(pixels)):
            x, y = pixels[p]
            col = self.current_ccp[y,x]
            ax.plot(
                self._wavel,
                self.datacube(frame)[:, y, x] / self.scalar + corr,
                color=col / 255,
                label=f"pixel@({str(x)},{str(y)})",
            )

        if hasattr(self, 'line_lim'):
            ax.set_xlim(self.line_lim)
        ax.set_xlabel(r"wavelength [$\rm\AA$]")
        ax.set_ylabel("intensity [units???]")
        ax.legend(fontsize=6, loc='lower left')
        plt.show()




    '''
    Set the quiet sun condition by assinging a patch on the frame to be quiet sun.
    '''
    def set_quiet_sun(self, frame, xlim, ylim, show=False, color='blue'):
        self.quiet_sun = {'frame':frame, 'xlim':xlim, 'ylim':ylim, 'color':color}
        self.quiet_spect=self.frame_integrated_spect(frame, xlim, ylim)
        self.ccp_frame(frame, Show=False)
        plt.imshow(self.current_ccp,origin='lower')
        plt.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], [ylim[0],ylim[0],ylim[1], ylim[1], ylim[0]], color=color)
        plt.show()

    def stand_dev_quiet_sun(self, option='wavl'):
        # option='area' or 'wavl'
        if option =='wavl':
            frame, xlim, ylim =  self.quiet_sun['frame'], self.quiet_sun['xlim'], self.quiet_sun['ylim']
            if hasattr(self, 'correction'):
                corr = self.correction
            else:
                corr = self._wavel*0
            self.std_quiet_sun = np.std(self.datacube(frame)[:,ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1], axis=(1,2))

    # def stand_dev_quiet_sun_interval(self, interval, option='wavl', to_print=False):
    #     # interval = [lower bound, upper bound] in angstroms
    #     wav = self._wavel

        # # option='area' or 'wavl'
        # if option=='wavl':
        #     if not hasattr(self, '.std_quiet_sun'):
        #         self.stand_dev_quiet_sun(option='wavl')
                
        #     try:
        #         stc = interp1d(wav, self.std_quiet_sun)
        #         integral, _ = scipy.integrate.quad(stc, interval[0], interval[1])
        #         l_intv = interval[1] - interval[0]
        #         appr_num = 
        #         std = integral / l_intv /

        #     stc = self.std_quiet_sun[(wav<= interval[1]) & (wav>= interval[0])]
        #     if to_print:
        #         print('take a look if those are not too much different?\n',stc)
        #     return np.mean(stc)

        # elif option=='area':
        #     frame, xlim, ylim =  self.quiet_sun['frame'], self.quiet_sun['xlim'], self.quiet_sun['ylim']
        #     return np.std(np.sum(self.datacube[frame,0,(wav<= interval[1]) & (wav>= interval[0]),ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1], axis=0))






    def fit_nessi_to_quiet_sun(self, theor_line, mu_data=None, initial_values=None):
        # theta = [horizontale translatie, verticale translatie, verticale schaalfactor]
        # theta = [0.2, 0.3, 0.89]

        if mu_data is None:
            mu_data = self.mu

        # to find the closed Mu:
        x = np.abs(theor_line.sst_mu-mu_data)
        index_mu  = np.where(x == np.min(x))[0][0]
        print("Best fitting mu is", theor_line.sst_mu[index_mu])

        f_nessi_clv = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[2] * theor_line.sst_dc*theor_line.sst_clv[index_mu]
                                             + theta[1], kind='linear', fill_value="extrapolate")
        f_nessi = lambda theta: interp1d(theor_line.sst_wav + theta[0], theta[2] * theor_line.sst_dc + theta[1]
                                        , kind='linear', fill_value="extrapolate")

        g = len(self._wavel)

        # dY = np.where(theor_Ha.sst_wav<6563.8, 0.01, 10) + np.where(6561.8<theor_Ha.sst_wav, 0.01, 10)
        #To simulate a specific domain around the well we cam make the errors on the wings huge

        data = [self._wavel,  self.quiet_spect, np.zeros(g)+0.001,np.zeros(g)+0.001]
        if initial_values is None:
            initial_values = np.array([-0.215, -0.111, 1.26])

        mini = da.optimalisatie(data, model=f_nessi_clv, beginwaarden=initial_values, fout_model=None, plot=False)
        theta = mini['x']
        self.theta_nessi_to_quiet_sun = theta
        da.plot_fit(data, model=f_nessi_clv, theta0=mini['x'], titel="fitting Nessi to sst data ",labelx=" $wavelength [\AA]$",
                    labely=" $intensity$  [arbitrary units]", figname=None , error=False)
        print(mini)
        da.kwaliteit_fit(data, mini)

        # self.f_nessi_clv = f_nessi_clv(theta)
        # self.f_nessi = f_nessi(theta)
        self.plot_fit_nessi_to_quiet_sun(self.quiet_sun['frame'], theor_line, f_nessi_theta=f_nessi(theta),
                                         f_nessi_clv_theta=f_nessi_clv(theta), xlim=self.quiet_sun['xlim'], ylim=self.quiet_sun['ylim'])
        return mini, f_nessi(theta), f_nessi_clv(theta)

    def plot_fit_nessi_to_quiet_sun(self, frame, theor_line, f_nessi_theta, f_nessi_clv_theta, xlim, ylim, restrict_to_line=False):

        fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))
        # fig.setitle("quiet sun determination, weiging constants")

        fig.suptitle('The gauge of the quiet sun spectrum', fontsize=16)

        self.ccp_frame(frame,Show=False)



        ax[0].set_title(
            f"spectral lines of SST frame {str(frame)}"
            + " H\u03B1 with normalized intensity"
        )
        self.frame_integrated_spect(frame)
        ax[0].plot(self._wavel, self.av_spect, label='sst data')
        ax[0].plot(self._wavel, self.quiet_spect, label='sst quiet sun') #

        theta = self.theta_nessi_to_quiet_sun

        ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_theta(theor_line.sst_wav + theta[0]), label='saas nessi')
        ax[0].plot(theor_line.sst_wav + theta[0], f_nessi_clv_theta(theor_line.sst_wav + theta[0]), label='nessi mu = 0.76')
        ax[0].legend()
        ax[1].imshow(self.current_ccp,origin='lower')
        ax[1].plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], [ylim[0],ylim[0],ylim[1], ylim[1], ylim[0]])
        ax[1].set_title(f"COCOplot of frame {str(frame)}")
        if hasattr(self, 'line_lim') and restrict_to_line:
            ax[0].set_xlim(self.line_lim)
        plt.show()


    '''
    Give a picture of the cocoplot with several posibbel area's which are quiet sun candidates.

        param
        ------
        x: list of two two tuple in a list which determines a rectangular patch by setting xlimits and ylimits.
            e.g. X = [[(750,940), (50,300)], [(650,975), (294,662)]] (2 patches)

    '''
    def possible_quiet_sun_patches(self, frame, theor_line, X, restrict_to_line=False):

        fig, ax = plt.subplots(nrows=1,ncols=2,figsize=(15, 7))

        fig.suptitle('Possible quiet sun patches', fontsize=16)

        self.ccp_frame(frame,Show=False)

        colors=['red', 'blue', 'yellow', 'orange', 'pink', 'purple', 'limegreen', 'darkgreen', 'gray']

        ax[0].set_title(
            f"spectral lines of SST frame {str(frame)}"
            + " H\u03B1 with normalized intensity"
        )
        self.frame_integrated_spect(frame)
        ax[0].plot(self._wavel, self.av_spect, '--', label='overal spectrum sst')
        theta = [0,0,1]
        ax[0].plot(theor_line.sst_wav + theta[0], theta[2] * theor_line.sst_dc*theor_line.sst_clv[12] + theta[1], '--', label='nessi mu = 0.76', color='black')

        ax[1].imshow(self.current_ccp,origin='lower')
        ax[1].set_title(f"COCOplot of frame {str(frame)}")
        for t, i in enumerate(X):
            xlim=i[0]
            ylim=i[1]
            color = (
                colors[t]
                if t < len(colors)
                else np.array(np.random.choice(range(256), size=3)) / 255
            )
            ax[0].plot(
                self._wavel,
                self.frame_integrated_spect(frame, xlim=xlim, ylim=ylim),
                color=color,
                label=f'{str(color)} area',
            )
            ax[1].plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], [ylim[0],ylim[0],ylim[1], ylim[1], ylim[0]], color=color)
        ax[0].legend()
        if hasattr(self, 'line_lim') and restrict_to_line:
            ax[0].set_xlim(self.line_lim)

        plt.show()

# ========================================================================
def get_extent(filename, timeFrame=0):
    """
    Read the coordinates of the corners to use them with imshow/extent

    Parameters
    filename : str
        name of the data cube
    timeFrame : int
        selected frame to get the coordinates for (Default value = 0)

    Returns
    -------
    extent_output : array_like
        1D array with solar coordinates in arcsec of the corners.

    Example
    -------
    >>> from ISPy.io import solarnet
    >>> extent = solarnet.get_extent(filename)
    >>> plt.imshow(data, extent=extent)

    :Author: 
        Carlos Diaz Baso (ISP/SU 2019)
    """
    io = f.open(filename)
    return [
        io[1].data[0][0][timeFrame, 0, 0, 0, 0],
        io[1].data[0][0][timeFrame, 0, 1, 1, 0],
        io[1].data[0][0][timeFrame, 0, 0, 0, 1],
        io[1].data[0][0][timeFrame, 0, 1, 1, 1],
    ]

def calculate_mu(x,y, radius=959.63):
    '''
    Calculates mu given x,y in arseconds the radius of the celestial object is option. 
    The value of the sun is 959.63 arcsec

    mu is the cos(theta) where theta is the angel between the zenit and the line of site of 
    a spot at the celestial objec.
    x,y are coordinates in arcsec (or any relative unit compared to the radius) from the middle
    of the celestial object.
    '''
    return (1-(x**2 + y**2)/radius**2)**0.5

def get_mu_mesh(filename, timeFrame=0, radius=959.63, shape=None):
    '''
    Creates a mesh of mu values with as many pixels as the observation. 
    '''
    io = f.open(filename)

    if shape is None:
        yp = io[0].data.shape[3]
        xp = io[0].data.shape[4]
    else:
        yp, xp = shape

    try:
        fxp = [io[1].data[0][0][timeFrame,0,0,0,0], io[1].data[0][0][timeFrame,0,1,1,0]]
        fyp = [io[1].data[0][0][timeFrame,0,0,0,1], io[1].data[0][0][timeFrame,0,1,1,1]]
    except IndexError as e:
        print('io:', io)
        print('np.shape(io[1].data):', np.shape(io[1].data))
        print('np.shape(io[1].data[0][0]):', np.shape(io[1].data[0][0]))
        raise e


    x = np.linspace(fxp[0], fxp[1], xp)
    y = np.linspace(fyp[0], fyp[1], yp)

    print(f"The frame is centered at {((fxp[0]+ fxp[1])/2,(fyp[0]+ fyp[1])/2)}")

    X, Y = np.meshgrid(x, y)

    MU = calculate_mu(X,Y, radius=radius)
    return MU, X, Y

from use_nessi import give_mu_contourplot
# def give_mu_contourplot(filename, timeFrame=0, over=None):
#     '''
#     Creates a contourplot of the mu values
#     if over is not None but a sst_data class object then the countour plot is shown over the specific frame
#     '''
#     shape = over.shape[3:5] if over is not None else None
#     MU, X, Y = get_mu_mesh(filename, timeFrame, shape=shape)
# #     extent = get_extent(filename, timeFrame)
#     if over is not None:
#         over.xx, over.yy = X, Y
# #     fig, ax = plt.subplots()
#     if over is not None:
#         over.frame_integrated_spect(timeFrame)
#         print(extent)
#         ax.imshow(Image.fromarray(over.current_ccp), origin='lower', extent=extent)
#         av_mu = np.average(MU, weights=over.boundary)
#         print(r'AVERAGE MU: The average $\mu$ wheigthed over the field of view is', av_mu, 'for timeframe', timeFrame)
#     CS = ax.contour(X, Y, MU)
#     ax.clabel(CS, inline=True, fontsize=10)
#     ax.set_title(r'Contour plot of the $\mu$ values.')
#     ax.set_xlabel('x-coordinate [arcsec]')
#     ax.set_ylabel('y-coordinate [arcsec]')

#     if over is not None:
#         return av_mu
    



# ========================================================================
def get_coord_creator(filename):
    """
    Converts pixel values to solar coordinates

    Parameters
    ----------
    filename : str
        name of the data cube
    pix_x, pix_y : int
        pixel location to convert
    timeFrame : int
        selected frame to the coordinates for (Default value = 0)

    Returns
    -------
    xy_output : list
        solar coordinates in arcsec

    Example
    -------
    >>> from ISPy.io import solarnet
    >>> [x_output, y_output] = solarnet.get_coord(filename, pix_x,pix_y)

    """
    io = f.open(filename)

    def coordinates(pix_x,pix_y, timeFrame=0):

        yp = [0, io[0].data.shape[3]]
        xp = [0, io[0].data.shape[4]]

        fxp = [io[1].data[0][0][timeFrame,0,0,0,0], io[1].data[0][0][timeFrame,0,1,1,0]]
        fyp = [io[1].data[0][0][timeFrame,0,0,0,1], io[1].data[0][0][timeFrame,0,1,1,1]]

        x_output = np.interp(pix_x, xp, fxp)
        y_output = np.interp(pix_y, yp, fyp) 

        return [x_output, y_output]
    
    return coordinates


def merge_wavelengths(wav1, wav2):
    wav1 = list(wav1)
    wav2 = list(wav2)
    return np.array(merge_sort_and_clip_lists(wav1, wav2))

def merge_sort_and_clip_lists(list1, list2):
    # Merge the two lists
    merged_list = list1 + list2

    # Sort the merged list
    return sorted(merged_list)

def is_point_inside_quadrilateral_grid(grid_shape, quadrilateral_vertices):
    def winding_number(n1, n2, vertices):
        wn = 0  # Winding number

        for i in range(len(vertices)):
            n1_1, n2_1 = vertices[i]
            n1_2, n2_2 = vertices[(i + 1) % len(vertices)]

            if n2_1 <= n2:
                if n2_2 > n2 and (n1_2 - n1_1) * (n2 - n2_1) - (n1 - n1_1) * (n2_2 - n2_1) > 0:
                    wn += 1
            else:
                if n2_2 <= n2 and (n1_2 - n1_1) * (n2 - n2_1) - (n1 - n1_1) * (n2_2 - n2_1) < 0:
                    wn -= 1

        return wn

    result = np.array([[winding_number(n1, n2, quadrilateral_vertices) != 0 for n2 in range(grid_shape[1])] for n1 in range(grid_shape[0])])
    return result


'''
This function appears to take an interval and an array x,
and then it returns a subset of x within the specified interval.
Additionally, it defines a nested function restx to restrict another
variable based on the same interval.
'''

def restrict_intervalx2(interval, x):
    # Find indices where x is within the specified interval
    l = np.where((x <= interval[1]) & (x >= interval[0]))[0]

    # Find the minimum and maximum indices within the interval
    m1 = np.min(l)
    m2 = np.max(l) + 1  # Adding 1 to include the upper bound

    # Extract the subset of x within the interval
    x2 = x[m1:m2]

    # Add the lower bound if it is not the first element
    if m1 > 0:
        x2 = np.insert(x2, 0, interval[0])

    # Add the upper bound if it is not the last element
    if m2 < len(x):
        x2 = np.append(x2, interval[1])

    # Define a nested function for restricting another variable based on the same interval
    def restx(y):
        # Extract the subset of y corresponding to the x interval
        y2 = y[m1:m2]

        # Add an extrapolated value at the upper bound if it is not the last element
        if m2 < len(x):
            t = y[m2 - 1] - (y[m2] - y[m2 - 1]) / (x[m2] - x[m2 - 1]) * (x[m2 - 1] - interval[1])
            y2 = np.append(y2, t)

        # Add an extrapolated value at the lower bound if it is not the first element
        if m1 > 0:
            t = y[m1] - (y[m1] - y[m1 - 1]) / (x[m1] - x[m1 - 1]) * (x[m1] - interval[0])
            y2 = np.insert(y2, 0, t)

        return y2

    # Return the restricted x and the function for restricting another variable
    return x2, restx


class WrongLineError(Exception):
    pass

####################################################################################
# OPTIMAL WIDTH PART
####################################################################################

def stand_dev_quiet_sun(self, option='wavl'):
    # option='area' or 'wavl'
    if option =='wavl':
        frame, xlim, ylim =  self.quiet_sun
        self.std_quiet_sun = np.std(self.datacube(frame)[:,ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1], axis=(1,2))

# def stand_dev_quiet_sun_interval(self, interval, option='wavl', to_print=False):
#     # interval = [lower bound, upper bound] in angstroms
#     wav = self._wavel

#     # option='area' or 'wavl'
#     if option=='wavl':
#         if not hasattr(self, '.std_quiet_sun'):
#             stand_dev_quiet_sun(self, option='wavl')

#         stc = self.std_quiet_sun[(wav<= interval[1]) & (wav>= interval[0])]
#         if to_print:
#             print('take a look if those are not too much different?\n',stc)
#         return np.mean(stc)/len(stc)**0.5

#     elif option=='area':
#         frame, xlim, ylim =  self.quiet_sun
#         return np.std(np.sum(self.datacube[frame,0,(wav<= interval[1]) &
#                                            (wav>= interval[0]),ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1],
#                                            axis=0))/(interval[1]-interval[0]) # here some squareroot?????


def Delta_I_over_I(self, interval, quiet_sun, flare_full_disk):
    # Here the thing is that the difference between the flare_full_disk and the quiet_sun is the same
    # as pure_flare - quiet_sun weighted with the area ratio.
    Delta_i = flare_full_disk



import scipy
def Delta_I_over_I2(wavl, quiet_sun, pure_flare, area_factor):
    x = wavl
    # print(quiet_sun, pure_flare)
    Delta_i = (pure_flare - quiet_sun)*area_factor
    Delta_i = scipy.integrate.simpson(y=Delta_i, x=x)

    I = scipy.integrate.simpson(y=quiet_sun, x=x, dx=1.0, axis=-1, even='avg')
    # print(Delta_i, I, Delta_i/I)
    return Delta_i / I


def get_noise(shape, mu=0, sigma=0.0001):
    return np.random.normal(mu, sigma, size=shape)

def ax_contourplot(fig, ax, X, Y, Z, x, line, decorations={}, seperate_colorbar=True, lim=0.001, logscale=False):
    if seperate_colorbar:
        pcm = ax.pcolormesh(X, Y, Z, cmap='RdBu_r',vmin=-np.max(Z),  shading='auto')
        fig.colorbar(pcm, ax=ax, extend='both')
    else:
        if logscale:
            pcm = ax.pcolormesh(X, Y, Z, cmap='RdBu_r', 
                            norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                              vmin=-lim, vmax=lim, base=10),
                                                shading='auto')
        else:
            pcm = ax.pcolormesh(X, Y, Z, cmap='RdBu_r', vmin=-lim, vmax=lim, shading='auto')


    if 'title' in decorations:
        ax.set_title(decorations['title'])
    # ax.set_ylabel('Time from '+str(decorations['time'])+' UT in [min]')
    # ax.set_xlabel(r'Wavelength [$\AA$]')

    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    if 'color' in decorations:
        color=decorations['color']
    else:
        color = 'black'
    # ax2.set_ylabel('Intensity []', color=color)  # we already handled the x-label with ax1
    ax2.plot(x,line, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    return pcm

def noise_plots(X,Y,Z, x, line, sigmas, seperate_colorbar=True, logscale=False):
    n = len(sigmas)
    rows = ((n-1)//3)+1
    cols = 3

    Z_noise=np.array([Z+get_noise(shape=np.shape(Z), mu=0, sigma=sigmas[i])  for i in range(n)])

    lim = np.max([-np.min(Z_noise), np.max(Z_noise)])


    fig, ax = plt.subplots(rows, cols, figsize=(5*cols,4*rows), constrained_layout=True) 
    fig.suptitle(r'Contrast profile of $H\alpha$ line, full disk', fontsize=20)
    for i in range(n):
        row = (i//3)
        col = i%3
        # noise = get_noise(shape=np.shape(Z), mu=0, sigma=sigmas[i]) 
        # print(np.shape(noise), np.shape(Z), np.shape(X), np.shape(Y))
        
        decorations = {'title':rf'noise: $\sigma=${sigmas[i]}.', 
               'time': '11:56:35', 'color':'black'}
        pcm = ax_contourplot(fig, ax[row, col], X,Y, Z_noise[i], x, line, decorations, 
                             seperate_colorbar=seperate_colorbar, lim=lim, logscale=logscale)

    if not seperate_colorbar:
        fig.colorbar(pcm, ax=ax[:, :], shrink=1/(rows+1))
    # fig.tight_layout()
    fig.supxlabel(r'Wavelength [$\AA$]', fontsize=16)
    fig.supylabel('Time from 11:56:35 UT in [min]', fontsize=16)
    # fig.supylabel('Intensity []', fontsize=16)
    plt.show()


######################################################################################################
#
# Here we do the analysis of the optimal width
#
######################################################################################################


from sklearn.isotonic import IsotonicRegression
import scipy
from matplotlib.collections import LineCollection



def stand_dev_quiet_sun(self, option='wavl', save=True):
    # option='area' or 'wavl'
    if option =='wavl':
        frame, xlim, ylim =  (self.quiet_sun['frame'], self.quiet_sun['xlim'], self.quiet_sun['ylim'])
        self.std_quiet_sun = np.nanstd(self.datacube(frame)[:,ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1], axis=(1,2))
    if save:
        filename = get_file_path_line_data(f'std_{self.name_of_line}')
        np.save(filename, self.std_quiet_sun)
        
def stand_dev_quiet_sun_interval(self, interval, option='wavl', to_print=False):
    # sourcery skip: extract-method, inline-immediately-returned-variable
    # interval = [lower bound, upper bound] in angstroms
    wav = self._wavel

    # option='area' or 'wavl'
    if option=='wavl':
        if not hasattr(self, '.std_quiet_sun'):
            stand_dev_quiet_sun(self, option='wavl')

        stc = interp1d(wav, self.std_quiet_sun)
        try:
            integral,_ = scipy.integrate.quad(stc, interval[0], interval[1])
            l_interv = interval[1]- interval[0]
            appr_num = len(self.std_quiet_sun) * l_interv / (wav[-1]-wav[0])
            std = integral / l_interv / appr_num**0.5
        except ValueError:
            stc = self.std_quiet_sun[np.where(wav <= interval[1]) & (wav >= interval[0])[0]]
            std = np.mean(stc) / len(stc) ** 0.5
        
        return std

    elif option=='area':
        frame, xlim, ylim =  (self.quiet_sun['frame'], self.quiet_sun['xlim'], self.quiet_sun['ylim'])
        return np.nanstd(np.sum(self.datacube(frame)[(wav<= interval[1]) &
                                           (wav>= interval[0]),ylim[0]:ylim[1]+1,xlim[0]:xlim[1]+1],
                                           axis=0))/(interval[1]-interval[0]) # here some squareroot?????
def candidate_intervals(sst_data, Deltas):
    peak = sst_data._wavel[sst_data.quiet_spect == np.min(sst_data.quiet_spect)][0] 
    intervals = np.array([peak-Deltas/2, peak+Deltas/2])
    return intervals.T

def STDs(sst_data, Deltas, intervals, show=True, inf_to=None):
    STD_Area = np.array([stand_dev_quiet_sun_interval(sst_data, intervals[i], option='area') for i,j in enumerate(Deltas)])
    STD_Wavl = np.array([stand_dev_quiet_sun_interval(sst_data, intervals[i], option='wavl') for i,j in enumerate(Deltas)])

    if show:
        plot_stds(Deltas, STD_Wavl, STD_Area)
    # print(f"{STD_Area[0] = }, {STD_Wavl[0] = }, {np.nan_to_num(STD_Area[0]) = }, {np.nan_to_num(STD_Wavl[0]) = }")

    if inf_to is None:
        inf_to = np.max(STD_Area[1]) * 2

    return np.nan_to_num(STD_Area, posinf=inf_to), np.nan_to_num(STD_Wavl, posinf=inf_to)


def plot_stds(Deltas, STD_Wavl, STD_Area):
    fig, ax = plt.subplots(ncols=2)
    ax[0].loglog(Deltas, STD_Wavl, label='option wavl')
    ax[0].loglog(Deltas, STD_Area, label='option area')
    ax[0].legend()

    ax[1].plot(Deltas, STD_Wavl, label='option wavl')
    ax[1].plot(Deltas, STD_Area, label='option area')
    ax[1].legend()
    plt.show()


def monotopic_STDs(STD_Area, STD_Wavl, Deltas, show=True):
    
    STD_Wavl2 = np.nan_to_num(STD_Wavl)
    STD_Area2 = np.nan_to_num(STD_Area)

    # Fit isotonic regression model
    ir = IsotonicRegression()
    STD_Wavl = -1*ir.fit_transform(Deltas, -1* STD_Wavl2)
    STD_Area = -1*ir.fit_transform(Deltas, -1*STD_Area2)

    if show:
        # Plotting the results
        plt.scatter(Deltas, STD_Wavl2, s=40, c="b", marker="x", label="data wavel")
        plt.plot(Deltas, STD_Wavl, label="Isotonic fit wavel", linewidth=3)
        plt.scatter(Deltas, STD_Area2, s=40, c="r", marker="x", label="data area")
        plt.plot(Deltas, STD_Area, label="Isotonic fit Area", linewidth=3)
        plt.xlabel("Number of measurements")
        plt.ylabel("Standard Deviation")
        plt.legend()
        plt.show()
    
    return STD_Area, STD_Wavl

def Delta_I_over_I2(wavl, quiet_sun, pure_flare, area_factor):
    x = wavl
    # print(quiet_sun, pure_flare)
    Delta_i = (pure_flare - quiet_sun)*area_factor
    Delta_i = scipy.integrate.simpson(y=Delta_i, x=x)

    I = scipy.integrate.simpson(y=quiet_sun, x=x, dx=1.0, axis=-1, even='avg')
    # print(Delta_i, I, Delta_i/I)
    return Delta_i / I

def restrict_intervalx2(interval, x):
    leftside = np.where(x<= interval[1])[0]
    rightside = np.where(x>= interval[0])[0]
     
    m1 = np.min(rightside)
    m2 = np.max(leftside)+1

    x2 = x[m1:m2]
    #adding the first
    if m1>0:
        x2 = np.insert(x2, 0,interval[0])
    #adding the last
    if m2<len(x):
        x2 =np.append(x2,interval[1])

    def restx(y):
        y2 = y[m1:m2]

        #adding the last
        if m2<len(x):
            t = y[m2-1] - (y[m2] - y[m2-1])/(x[m2] - x[m2-1]) * (x[m2-1]-interval[1])
            y2 = np.append(y2, t )

        #adding the first
        if m1>0:
            t = y[m1] - (y[m1] - y[m1-1])/(x[m1] - x[m1-1]) * (x[m1]-interval[0])
            y2 = np.insert(y2,0, t )
        return y2

    return x2, restx

def signal_of_interval(sst_data, FOV_spectrum, intervals, area_factor=(60**2/4/np.pi/959.63**2)):
    B =[] # signal

    for interval in intervals:
        print(f"Now calculating for interval {interval}", end='\r')

        wavl = sst_data._wavel
        quiet_sun = sst_data.quiet_spect
        w2, restx = restrict_intervalx2(interval, wavl)
        q2 = restx(quiet_sun)

        DI_I = []
        frames = range(0, sst_data._number_of_frames)

        for frame in frames:
            f2 = restx(FOV_spectrum[frame])
            DI_I.append(Delta_I_over_I2(w2, q2, f2, area_factor))

        B.append(np.array(DI_I))

    return np.array(B)

def signal_to_noise(signal, STD_Area, STD_Wavl):
    A_area = signal/STD_Area[:, np.newaxis]
    A_wavl = signal/STD_Wavl[:, np.newaxis]
    return A_area, A_wavl

def multiline(xs, ys, c, ax=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    c : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # find axes
    ax = plt.gca() if ax is None else ax

    # create LineCollection
    segments = [np.column_stack([x, y]) for x, y in zip(xs, ys)]
    lc = LineCollection(segments, **kwargs)

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    lc.set_array(np.asarray(c))

    # add lines to axes and rescale
    #    Note: adding a collection doesn't autoscalee xlim/ylim
    ax.add_collection(lc)
    ax.autoscale()
    return lc

def hulp_time(string):
    if len(string) < 6:
        string = string + '00000000000'
    return float(string[:2])*60 + float(string[3:5]) + float(string[6:8])/60

def time_to_minutes(_time):
    return np.array([ hulp_time(t) for i,t in enumerate(_time) ])

def get_TIME(sst_data):
    TIME = time_to_minutes(sst_data._time)
    TIME -= hulp_time(official_start_flare(sst_data.name_of_line, except_time=sst_data._time[0]))
    sst_data.TIME = TIME
    return TIME

def minutes_to_time_string(minutes, start=0):
  """Converts minutes in float to a string in 'hh:mm:ss' format.

  Args:
      minutes: A float representing the number of minutes.

  Returns:
      A string in 'hh:mm:ss' format.
  """

  # Ensure positive minutes (optional)
  minutes += start

  # Extract hours, minutes, and seconds
  hours = int(minutes // 60)
  minutes = int(minutes % 60)
  seconds = int((minutes % 60) * 60) % 60

  # Format time string with leading zeros
  time_string = f"{hours:02d}:{minutes:02d}:{seconds:02d}"

  return time_string

def optimal_ax(ax, A, STD, TIME, Deltas, criterion, timelim=30):
    MAX = np.nanmax(A)

    i = np.where(A==MAX)[0][0]
    j  = np.where(A==MAX)[1][0]
    print(f'The maximal signal to noise ratio for the {criterion} criterion is if the interval is of width $\Delta x={Deltas[i]}\AA$.')
    print(f'at an intesity of {A[i,j]}. For a standard deviation of {STD[i]}.\n')
    ax.plot([TIME[j], TIME[j]+5], [A[i,j],A[i,j]], color ='black')
    ax.text(TIME[j]+5, A[i,j], ''+str(round(Deltas[i], 2))+r' $\AA$', size=6)
    n_times_time = np.array([TIME[:np.shape(A)[1]] for i in Deltas])
    lc = multiline(n_times_time, A, Deltas, ax=ax, cmap='jet', lw=2)
    ax.set_title(f' STD by the {criterion} criterion ')
    ax.set_xlabel('t [minutes]')
    ax.set_xlim(0,timelim)
    return lc

def optimum_interval_width_plots(A_area, A_wavl, STD_area, STD_wavl, Deltas, TIME, timelim=30):
    fig, ax = plt.subplots(ncols=2, figsize=(14,5))
    fig.suptitle(r'Time evolution of the $\Delta I/\sigma I$ as function of width of the interval.')
    fig.supylabel(r'$\Delta I/\sigma I$ []')
    ld = optimal_ax(ax[0], A_area, STD_area, TIME, Deltas, criterion='area', timelim=timelim)
    lc = optimal_ax(ax[1], A_wavl, STD_wavl, TIME, Deltas, criterion='wavl', timelim=timelim)
    axcb = fig.colorbar(lc, ax=ax[:])
    axcb.set_label(r'Width of interval $\Delta x$')
    plt.show()
    
def plot_peak_enhancement(A_area, A_wavl, Deltas):
    fig, ax = plt.subplots(ncols=2, figsize=(14,5))
    fig.suptitle(r'peak signal to noice enhancement as a function of $\Delta W$.')
    fig.supylabel(r'$\Delta I/\sigma I$ []')
    ax_peak_enhancement(ax[0], A_area, Deltas, criterion='area')
    ax_peak_enhancement(ax[1], A_wavl, Deltas, criterion='wavl')
    plt.show()
    
    
def ax_peak_enhancement(ax, A, Deltas, criterion):
    MAX = np.nanmax(A, axis=1)
    
    ax.plot(Deltas, MAX)
    ax.set_title(f' (STD by the {criterion} criterion) ')
    ax.set_xlabel(r"Width of the interval $\Delta W$ $[\AA]$")
    
def get_file_path_opt_w_data(filename):
    if filename[-4:] != ".npy":
        filename = f"{filename}.npy"
    return os.path.join('line_data', 'optimal_width',filename)

def hulp_save_opt_w_data(name, array):
    filename = get_file_path_opt_w_data(name)
    np.save(filename,array)
    
def make_opt_w_dir():
    # Specify the directory path
    directory = "line_data/optimal_width"

    # Create the directory (raises an error if it already exists)
    try:
        os.mkdir(directory)
        print(f"Directory '{directory}' created successfully.")
    except FileExistsError:
        print(f"Directory '{directory}' already exists.")


def saves_As(name_of_line, A_area, A_wavl, Deltas):
    make_opt_w_dir()
    # save A_area, A_wavl, Deltas, ...
    hulp_save_opt_w_data("A_area_"+name_of_line, A_area)
    hulp_save_opt_w_data("A_wavl_"+name_of_line, A_wavl)
    hulp_save_opt_w_data("Deltas_"+name_of_line, Deltas)

def analyse_optimal_interval(sst_data, Deltas=np.arange(0, 3, 0.03), area_factor=(60**2/4/np.pi/959.63**2)):
    # if Deltas is None:
    #     Deltas = np.arange(0, 3, 0.03)
    intervals = candidate_intervals(sst_data, Deltas)

    STD_Area, STD_Wavl = STDs(sst_data, Deltas, intervals)
    print(f'{STD_Area = }, {STD_Wavl = } ')
    STD_Area, STD_Wavl = monotopic_STDs(STD_Area, STD_Wavl, Deltas)

    sst_data.calculate__FOV_spect_over_time()

    signal = signal_of_interval(sst_data, sst_data.FOV_spectrum, intervals, area_factor)
    A_area, A_wavl = signal_to_noise(signal, STD_Area, STD_Wavl)

    sst_data.TIME = get_TIME(sst_data)

    optimum_interval_width_plots(A_area, A_wavl, STD_Area, STD_Wavl, Deltas, sst_data.TIME, timelim=25)
    
    plot_peak_enhancement(A_area, A_wavl, Deltas)
    
    saves_As(sst_data.name_of_line, A_area, A_wavl, Deltas)
    
    return A_area, A_wavl, STD_Area, STD_Wavl, Deltas, sst_data.TIME 


from use_nessi import get_full_path, correct_flare_start, load_for_further_analysis, save_for_further_analysis
from use_nessi import add_enters, Movie_making, time_hulp