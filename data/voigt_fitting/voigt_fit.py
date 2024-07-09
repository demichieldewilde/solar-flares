from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
sys.path.append('E:/solar flares/data/2017-09-06')
import use_nessi3 as un
import data_analysis as da
import numpy as np
import matplotlib.pyplot as plt
import use_nessi as un
import use_nessi2 as un2
import os

def gaussian(x, amp_g, cen_g, sigma_g):
    """Gaussian function."""
    return amp_g * np.exp(-(x - cen_g)**2 / (2 * sigma_g**2))

def lorentzian(x, amp_l, cen_l, gamma_l):
    """Lorentzian function."""
    return amp_l * (gamma_l**2 / ((x - cen_l)**2 + gamma_l**2))

def voigt(x, amp_v, cen_v, sigma_v, gamma_v, offset=0):
    """Voigt profile: approximation by combining Gaussian and Lorentzian."""
    f_g = sigma_v**2 / (sigma_v**2 + gamma_v**2)
    return (f_g * gaussian(x, 1, cen_v, sigma_v) +
            (1 - f_g) * lorentzian(x, 1, cen_v, gamma_v)) * amp_v + offset

def voigt_plus_gaussian(x, amp_v, cen_v, sigma_v, gamma_v, amp_g, cen_g, sigma_g):
    """Combination of a Voigt profile and an additional Gaussian."""
    return voigt(x, amp_v, cen_v, sigma_v, gamma_v) + gaussian(x, amp_g, cen_g, sigma_g)

def fit_voigt_and_gaussian(x_data, y_data, initial_guess):
    """Fit a Voigt and Gaussian profile to the data."""
    popt, pcov = curve_fit(voigt_plus_gaussian, x_data, y_data, p0=initial_guess)
    return popt, np.sqrt(np.diag(pcov))

def fit_voigt(x_data, y_data, initial_guess, fix_ind=[]):
    """Fit a Voigt profile to the data."""
    def model(theta):
        for i in fix_ind:
            theta[i] = initial_guess[i]
        return voigt(theta)
    popt, pcov = curve_fit(voigt, x_data, y_data, p0=initial_guess)
    return popt, np.sqrt(np.diag(pcov))

def plot_data_fit_and_components(x_data, y_data, initial_guess, popt=None):
    """Plot the data, fit, and components."""
    if popt is None:
        popt, _ = fit_voigt_and_gaussian(x_data, y_data, initial_guess)
    y_fit = voigt_plus_gaussian(x_data, *popt)
    y_voigt = voigt(x_data, popt[0], popt[1], popt[2], popt[3])
    y_gaussian = gaussian(x_data, popt[4], popt[5], popt[6])
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_data, y_data, 'b.', label='Data')
    plt.plot(x_data, y_fit, 'r-', linewidth=2, label='Voigt + Gaussian Fit')
    plt.plot(x_data, y_voigt, 'g--', label='Voigt Component')
    plt.plot(x_data, y_gaussian, 'm--', label='Gaussian Component')
    plt.title('Data, Fit, and Components')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()
    
def exc_wings(left, rigth):
    return range(rigth, left+1)

def plot_data_fit_voigt(x_data, y_data, initial_guess, popt=None, neglect_points=[]):
    """Plot the data, fit, and components."""
    if popt is None:
        popt, _ = fit_voigt(x_data, y_data, initial_guess)
    y_voigt = voigt(x_data, popt[0], popt[1], popt[2], popt[3])
    y_diff = y_data - y_voigt
    plt.figure(figsize=(10, 6))
    plt.plot(x_data, y_data, 'y.', label='Data neglected')
    plt.plot(np.delete(x_data,neglect_points, axis=0), np.delete(y_data,neglect_points, axis=0), 'b.', label='Data for fit')
    plt.plot(x_data, y_voigt, 'b--', label='Voigt Component')
    a = np.average(y_data)
    plt.plot(x_data, y_diff + 1.01*a , 'r--', label='residue (shifted)')
    plt.title('Data, Fit, and Components')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.show()

def test_voigt_plus_gauss():
    # Generate synthetic data
    x = np.linspace(-5, 5, 100)
    y = voigt_plus_gaussian(x, 1, 0, 1, 0.5, 0.5, -2, 1) + np.random.normal(0, 0.1, len(x))

    # Initial guess: [amp_v, cen_v, sigma_v, gamma_v, amp_g, cen_g, sigma_g]
    initial_guess = [1, 0, 1, 1, 0.5, 0, 1]

    # Plot the data, the fit, and the individual components
    plot_data_fit_and_components(x, y, initial_guess)
    
def contrast_fit_voigt_gauss(wav, contrast, initial_guess, plot_rate=1000):
    param_fit = [] 
    residuals = []
    
    for i in range( np.shape(contrast)[0]):
        param_fit.append(fit_voigt_and_gaussian(wav, contrast[i], initial_guess))
        # print(param_fit[i])
        
        if i%plot_rate == plot_rate-1:
            plot_data_fit_and_components(wav, contrast[i], initial_guess)
            
        popt = param_fit[-1][0]
        y_data = contrast[i]
        y_voigt = voigt(wav, popt[0], popt[1], popt[2], popt[3])
        y_diff = y_data - y_voigt
        residuals.append(y_diff)
        
    return param_fit, residuals
        
def average_last_of_params(params, num=5, initial_guess=None, fix_ind=[]):
    if len(params) <= num:
        k =  [p[0] for p in params]
    else:
        k = [p[0] for p in params[-num : ]]
    if initial_guess is not None:
        k.append(initial_guess)
    guess= np.average(np.array(k), axis=0)
    for ind in fix_ind:
        guess[ind] = initial_guess[ind]
    return guess

def retreive_initial_guess(initial_guess, frame, many_guesses=True):
    """retreives the initial guess from the correct format
    when a list of four is given this is the initial guess for alle frames
    otherwise the format is 
    [[framenumber start, associated initial guess],
    [framenumber start, associated initial guess],
    ...
    """
    if np.shape(initial_guess) == (4,):
        many_guesses = False
    if not many_guesses:
         return initial_guess
    for i, j in enumerate(initial_guess):
        if frame < j[0]:
            return initial_guess[i-1][1]


def contrast_fit_voigt(wav, contrast, initial_guess, plot_rate=1000, neglect_points=[], 
                       fix_ind=[], exclude_frames=[], hard_fix=False):

    param_fit = []
    y_fits = []
    residuals = []
    n_frames = np.shape(contrast)[0]
    
    # parameter fitting
    for i in range(n_frames):
        i_g = retreive_initial_guess(initial_guess, i)

        if i in exclude_frames:
            k = i_g
            k[0] = 0
            popt = average_last_of_params(param_fit, num=5, initial_guess=k, fix_ind=fix_ind)
            param_fit.append([popt, popt])
        else:
            average_guess = average_last_of_params(param_fit, num=5, initial_guess=i_g, fix_ind=fix_ind)
            try:            
                w = np.delete(wav, neglect_points, axis=0)
                C = np.delete( contrast[i], neglect_points, axis=0)
                hard_fix_ind = fix_ind if hard_fix else []
                F = fit_voigt(w, C, average_guess, fix_ind)
                param_fit.append(F)

                if i%plot_rate == 0:
                    print(f'frame {i} with {param_fit[-1] = }. Here comes the plot:')
                    plot_data_fit_voigt(wav, contrast[i], average_guess, popt=F[0], neglect_points=neglect_points)

            except RuntimeError:
                if len(param_fit) != 0:
                    print(f"At frame {i} the voigt fitting was not succesfull and average guess ({average_guess}) is used as params and as std of params. ")
                    param_fit.append((average_guess, average_guess))
                else:
                    print(f"At frame {i} the voigt fitting was not succesfull and previous params are used. ")
                    param_fit.append(param_fit[-1])              
    
    # parameter smoothing
    param_fit = np.array(param_fit)
    param_fit = np.stack((complete_rolling_average(param_fit[:,0,:]), complete_rolling_average(param_fit[:,1,:])), axis=1)
    # print('it is smoothent')

    # residual calculation
    for i in range( n_frames):
        y_data = contrast[i]
        popt = param_fit[i][0]
        y_voigt = voigt(wav, popt[0], popt[1], popt[2], popt[3])
        y_fits.append(y_voigt)
        y_diff = y_data - y_voigt
        residuals.append(y_diff)
    try:
        return param_fit, np.array(y_fits) , np.array(residuals)
    except ValueError as error:
        _extracted_from_contrast_fit_voigt_35(param_fit, y_fits, residuals, error)

def _extracted_from_contrast_fit_voigt_35(param_fit, y_fits, residuals, error):
    print('The lenght of those arrays does not work out')
    t = [np.shape(k) for k in param_fit]
    print("the shapes of the param_fits:", set(t), t)
    t = [np.shape(k) for k in y_fits]
    print("the shapes of the y_fits:", set(t), t )
    t = [np.shape(k) for k in residuals]
    print("the shapes of the residuals:", set(t) , t)
    raise error
       
def complete_rolling_average(param, k=30):
    kernel = np.ones(k) / k
    pre = np.array([param[0] for _ in range((k-1)//2) ])
    aft = np.array([param[-1] for _ in range((k)//2)])
    param = np.concatenate( (pre ,  param, aft))
    param = np.array([np.convolve(param[:, i], kernel, mode='valid') for i in range(4) ]).T
    return param
 
def make_analysis(name, data, initial_guess, plot_rate=50, offset=0, neglect_points=None, 
                  fix_ind=None, hard_fix=False, exclude_frames=None):
    offset -= 1 
    if neglect_points is None:
        neglect_points = []
    if fix_ind is None:
        fix_ind = []
    wav, contr , time, line, std = un2.contrast_FD_data(name, data, quiet_sun_subtraction=False, num=40)
    print(f'The average is {np.average(contr)}')
    
    if exclude_frames is None:
        exclude_frames = np.where(un2.most_quiet_frames(name, time))[0] 
    # print(exclude_frames)
    # print(np.shape(wav), np.shape(DFOV), np.shape(time), np.shape(line))
    # Initial guess: [amp_v, cen_v, sigma_v, gamma_v, amp_g, cen_g, sigma_g]
    # print(np.shape(wav), np.shape(np.delete(wav, neglect_points, axis=0)), wav, np.delete(wav, neglect_points, axis=0))
    # print(np.shape(DFOV+offset), np.shape(np.delete(DFOV+offset, neglect_points, axis=1)), DFOV+offset, np.delete(DFOV+offset, neglect_points, axis=1))
    params, voigt, res = contrast_fit_voigt(wav, contr+offset, initial_guess, plot_rate, 
                                              neglect_points, fix_ind, exclude_frames=exclude_frames, hard_fix=hard_fix)
    visualize_analysis(res, voigt, wav, time, name)
    save_voigt_fits(name, params, res)
        
def display_OK():
    ok = [
        "  ____    _   __",
        " / __ \\  | | / /",
        "| |  | | | |/ /",
        "| |_ | | |   \\",
        " \____/  |_|\_\\"
    ]
    for line in ok:
        print(line * 5)  # Print 

    
    
def save_voigt_fits(name, params, res):
    fname = f"fit_data/voigt_data_{name}.npz"
    np.savez(fname, params, res)

def load_voigt_data(name):
    fname = f"D:/solar flares/data/voigt_fitting/fit_data/voigt_data_{name}.npz"
    data =  np.load(fname) 
    params = data['arr_0']
    res = data["arr_1"]
    return params, res    

def cut_off_data(data, up_lim=None, down_lim=None):
    if up_lim is None and down_lim is None:
        up_lim = max(-np.percentile(data, 3), np.percentile(data, 97))
        down_lim = - up_lim
    elif up_lim is None:
        up_lim = np.percentile(data, 97)
    elif down_lim is None:
        down_lim = np.percentile(data, 3)
    # Apply cutoff using boolean indexing and conditional assignment
    cutoff_data = np.where(data > up_lim, up_lim, data)
    cutoff_data = np.where(data < down_lim, down_lim, cutoff_data)
    return cutoff_data

def visualize_analysis(res, voigt, wav, time, name, non_centered=True):
    fig, ax = plt.subplots()

    vmax = np.percentile(res, 97)
    vmin = np.percentile(res, 3)

    print(f"{vmax = }, {vmin = }")
    c = ax.imshow(np.array(res), aspect="auto", cmap='RdBu_r', origin='lower', extent=(wav[0], wav[-1], time[0], time[-1]),
                vmax=vmax, vmin=vmin)
    #     # voigt = cut_off_data(voigt, up_lim=vmax, down_lim=vmin)
    # else:
    #     c = ax.imshow(np.array(res), aspect="auto", cmap='RdBu_r', origin='lower', extent=(wav[0], wav[-1], time[0], time[-1]))
    # # pcm = ax.pcolormesh(X, Y, Z, cmap='RdBu_r',vmin=-np.max(Z),  shading='auto')
    cb = fig.colorbar(c, ax=ax, extend='both')
    X, Y = np.meshgrid(wav, time)
    
    CS = ax.contour(X, Y, voigt, colors='black', alpha=0.5)
    ax.clabel(CS, inline=True, fontsize=10) 
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel('time from start of flare [min]')
    ax.set_title(f"voigt fit and Residue analysis for {name} flare")
    cb.set_label(r'the residues after Voigt fit [relative intensity]')
    plt.show()   
    
initial_guesses = {'Ha':[5.55265545e-01, 6.56293638e+03, 5.29257897e-01, 4.96399919e-01], 
                   "CaK":[4.5, 3933.675, -2.66732183e-01, 2.2284e-1], 
                   "CaIR":[1.5,  8.54203173e+03, -2.74932183e-01, -3.61301961e-04]}



element_from_name = un2.element_from_name

def list_or_dict(ob, i, key, default=0):
    if isinstance(ob, list):
        return ob[i]
    elif isinstance(ob, dict):
        return ob.get(key, default)
    elif isinstance(ob, int):
        return ob
    else:
        print("Could not get a good type retruning Deafault")
        return default

def make_full_analysis(data, names_list, offsets=None, plot_rate={}, init_guesses=initial_guesses, list_neglect_points=None, fix_ind={}):
    if list_neglect_points is None:
        list_neglect_points = [[] for _ in range(len(names_list))]
    if offsets is None:
        offsets = np.zeros(len(names_list))
    for i, name in enumerate(names_list):
        el = element_from_name(name)
        make_analysis(name, data, init_guesses[el], list_or_dict(plot_rate, i, el, 500), 
                      offset=list_or_dict(offsets, i, el, 500), neglect_points=list_or_dict(list_neglect_points, i, el, []), 
                      fix_ind=list_or_dict(fix_ind, i, el, []))