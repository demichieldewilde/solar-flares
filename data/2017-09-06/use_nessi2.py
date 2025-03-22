import numpy as np
import matplotlib.pyplot as plt
import use_nessi as un
import os
from importlib import reload
from matplotlib import cm
import matplotlib.cbook as cbook
import matplotlib.colors as mplcolors
from scipy.interpolate import interp1d
from scipy.signal import convolve2d


sr = solar_radius = 959.63
area_factor = 60**2/np.pi/sr**2

def element_from_name(name):
    if 'CaK' in name:
        return 'CaK(2)' if '(2)' in name else 'CaK'
    lines = ['Ha', 'CaI', "Fe6173", "Hbeta", "He"]
    for line in lines:
        if line in name:
            return "CaIR" if line=="CaI" else line
    if "Fe" in name:
        print(f"Assuming that {name} is spectral line Fe6173.")
        return "Fe6173"
    raise ValueError(f'The given name {name} is not known as a spectral name (yet).')

def disgard_cont_point(name, data):
    # FOV_spectrum
    data[f'FOV_{name}'] = data[f'FOV_{name}'][:, :-1]

    # load quiet_sun profile
    data[f'quiet_sun_{name}'] = data[f'quiet_sun_{name}'][:, :-1]
    
    # correcting nessi profile
    data[f'nessi_{name}'] = data[f'nessi_{name}'][:,:-1]

def split_data_in_two_lines(name, data, wav_split, lines):
    lines.append(f"{name}(2)")
    
    wav_sst = data[f'quiet_sun_{name}'][0]
    ind_sst = np.where(wav_sst > wav_split)[0][0]
    
    # splitting FOV_spectrum
    data[f'FOV_{name}(2)'] = data[f'FOV_{name}'][:, ind_sst:]
    data[f'FOV_{name}'] = data[f'FOV_{name}'][:, :ind_sst]

    # splitting quiet_sun profile
    data[f'quiet_sun_{name}(2)'] = data[f'quiet_sun_{name}'][:, ind_sst:]
    data[f'quiet_sun_{name}'] = data[f'quiet_sun_{name}'][:, :ind_sst]
    
    wav_nessi = data[f'nessi_{name}'][0]
    ind_nessi = np.where(wav_nessi > wav_split)[0][0]

    # splitting nessi profile
    data[f'nessi_{name}(2)'] = data[f'nessi_{name}'][:,ind_nessi:]
    data[f'nessi_{name}'] = data[f'nessi_{name}'][:,:ind_nessi]

def name_no2(name):
    return name[:-3] if name[-3:] == '(2)' else name
        
sr = solar_radius = 959.63
area_factor = 60**2/np.pi/sr**2    

def test_contrast(data, name_of_line):
    print(f'Test of contrast profile and associated curves of {name_of_line}')
    wav, DFD, time, line, std = difference_FD_data(name_of_line, data, False)
    plt.plot(wav, DFD[0]  , label='DFD_0')
    plt.plot(wav, DFD[-1] , label='DFD_-1')
    plt.legend()
    plt.show()
    FOV = data[f"FOV_{name_of_line}"]
    wav_qs, qs_spec, std_qs = data[f"quiet_sun_{name_of_line}"]
    time = data[f"TIME_{name_no2(name_of_line)}"]
    wav_nessi, fov_nessi, saas_nessi = data[f"nessi_{name_of_line}"]
    plt.plot(wav_qs, FOV[0], label='FOV 0')
    plt.plot(wav_qs, FOV[-1], label='FOV -1')
    plt.plot(wav_qs, qs_spec, label='quiet sun')
    plt.plot(wav_nessi, fov_nessi, label='FOV nessi')
    plt.plot(wav_nessi, saas_nessi, label='saas nessi')
    plt.legend()
    plt.show()
    wav, contr, time, line, std = contrast_FD_data(name_of_line, data, quiet_sun_subtraction=False, num=100,area_factor=60**2/np.pi/959.63**2, add_noise=False)
    plt.plot(wav, contr[0], label='contr 0')
    plt.plot(wav, contr[-1], label='contr -1')
    plt.legend()
    plt.show()

def most_quiet_flare_time(name):
    """Returns the indices over which to average to get a good quiet pattern to calculate the contrast profile from.

    Args:
        name (string): name of the line

    Returns:
        [t0, t1]: begin and end time of the quiet flare time
    """
    if '19'in name:
        return [-10,-5]
    elif '13' in name:
        return [60,65]
    elif '9u' in name:
        return [50,61]
    elif "17a" in name:
        return [120,135]
    elif "17" in name:
        return [43,50]
    elif "15b" in name:
        return [5,10]    
    elif "15a" in name:
        return [35,40]
    elif "15" in name:
        return [70,80]
    elif "14a" in name:
        return [-90, -75]
    elif "14" in name:
        return [4.5,6]
    elif "23a" in name:
        return [10,13]
    elif "23" in name:
        return [55,60]
    elif "16" in name:
        return [55, 60]
    elif "21a" in name:
        return [45, 52]    
    elif "21" in name:
        return [25, 27]
    elif "24a" in name:
        return [45,55]
    elif "22a" in name:
        return [47, 57] 
    elif "22" in name:
        return [55,65]
    else:
        raise NameError(f'WRONG NAME: the line {name} had no most quiet flare time defined.')

def most_quiet_frames(name_of_line, time):
    T = most_quiet_flare_time(name_of_line)
    try:
        assert(T[-1] > time[0])
        assert(T[0]< time[-1])
    except AssertionError as e:
        print("ERROR Wrong quiet times")
        print("time is ", time)
        print("quiet frames are ", T)
        raise(e)
    return np.where(time >= T[0], 1,0) * np.where( time <= T[1], 1, 0)
    
def most_quiet_FD_spectr(FD_spec, name_of_line, time):
    R = most_quiet_frames(name_of_line, time)
    return np.average(FD_spec, axis=0, weights=R)
    
def contrast_FD_data(name_of_line, data, quiet_sun_subtraction=False, num=100,area_factor=60**2/np.pi/959.63**2, add_noise=False, theoretical=False): 
    wav, DFD, time, line, std = difference_FD_data(name_of_line, data, quiet_sun_subtraction, num, area_factor, add_noise)
    saas = line
    
    FD = DFD + line
    if theoretical:
        mq_FD = line
    else:
        mq_FD = most_quiet_FD_spectr(FD, name_of_line, time)
    
    # plt.plot(wav, line, label="NESSI")
    # plt.plot(wav, most_quiet_FD_spectr(FD, name_of_line, time), '--', label="Most quiet FD")
    # plt.legend()
    # plt.show()
    
    # plt.plot(wav, most_quiet_FD_spectr(FD, name_of_line, time) - line, '--', label="Difference")
    # plt.legend()
    # plt.show()
    
    contr_prof = FD / mq_FD
    
    
    return wav, contr_prof, time, line, (std/mq_FD**0.5 if std is not None else None)

def normalized_spectral_change(name_of_line, data, quiet_sun_subtraction=False, num=100,area_factor=60**2/np.pi/959.63**2, add_noise=False, theoretical=False): 
    wav, DFD, time, line, std = difference_FD_data(name_of_line, data, quiet_sun_subtraction, num, area_factor, add_noise)
    DS = DFD / line if theoretical else DFD / most_quiet_FD_spectr(line + DFD, name_of_line, time)
    return wav, DS, time, line, std/line

def difference_FD_data(name_of_line, data, quiet_sun_subtraction=False, num=100,area_factor=60**2/np.pi/959.63**2, add_noise=False):            
    wav, DFOV, time, line, std = difference_FOV_data(name_of_line, data, quiet_sun_subtraction, num)

    # Correct normalization for area and mu-value (all intensities are normalized on the first wavelength)
    #   -correction for area is area_factor
    # thereby we have area_factor * clv(at wav_qs[0])/1 thus
    DFD = area_factor * DFOV
    
    if add_noise:
        noise = np.random.normal(loc=0, scale=std[0], size=(DFD.shape))
        DFD += noise

    return wav, DFD, time, line, std * area_factor**0.5

'''
The scale_pix_to_saas is the number by which the standard deviation has to be multiplied to get to the the std of the saas observation
it should be the inverse of square root of (number of pixels in sst (1920×1200 CHROMIS or 1k × 1k CRISP) * areafactor (60**2/np.pi/sr**2) ) 
hence 1/53.544360373477076
'''
def difference_FOV_data(name_of_line, data, quiet_sun_subtraction=False, num=100, scale_pix_to_saas=1/53.4):
    FOV = data[f"FOV_{name_of_line}"]
    wav_qs, qs_spec, std_qs = data[f"quiet_sun_{name_of_line}"]
    time = data[f"TIME_{name_no2(name_of_line)}"]
    wav_nessi, fov_nessi, saas_nessi = data[f"nessi_{name_of_line}"]
    std_qs *= scale_pix_to_saas

    if quiet_sun_subtraction:
        qs_spectc =  qs_spec
        wav_qsc = wav_qs
    else:
        qs_spectc = fov_nessi
        wav_qsc = wav_nessi

    wav = smooth_wavelengths(wav_qs, wav_nessi, num)
    
    DFOV = np.array([interp1d(wav_qs, FOV[i,:])(wav)  for i in range(np.shape(FOV)[0])]) - interp1d(wav_qsc, qs_spectc)(wav)

    line = interp1d(wav_nessi, saas_nessi)(wav)
    std = interp1d(wav_qs, std_qs)(wav)
    # Correct normalization for mu-value has already been done by using the gauge to NESSI
    # Area normalization should still happen though!
    return wav, DFOV , time, line, std

def smooth_wavelengths(wav1, wav2, num=100):
    if num <= 0:
        num = len(wav1)
    start = max(np.min(wav1), np.min(wav2))
    stop  = min(np.max(wav1), np.max(wav2))
    return np.linspace(start, stop, num)
    
def get_harps_std(name, backup=1):
    harps = get_Harps(name)
    if harps[0] is None:
        return backup
    wingl = harps[0][:, :200]
    wingr = harps[0][:, 800:]
    stdl = np.std(wingl) 
    stdr = np.std(wingr)
    return np.mean([stdl, stdr])

def noise_alike_harps(name, shape, backup=1, x_steps=1, scale_noise=1):
    std = get_harps_std(name, backup) * scale_noise
    return np.array(
        [
            generate_random_array(shape[1], block_size=x_steps, std=std)
            for _ in range(shape[0])
        ]
    )

def generate_random_array(size, block_size, std):
    """
    Generates a random array with size 'size' where values are repeated in blocks 
    of size 'block_size' and randomly sampled from a normal distribution with standard deviation 'std'.

    Args:
        size: The total size of the desired array.
        block_size: The size of each block where values will be the same.
        std: The standard deviation of the normal distribution for random sampling.

    Returns:
        A NumPy array with the desired properties.
    """
    # Ensure block_size is a divisor of size
    # if size % block_size != 0:
    #     raise ValueError("Size must be a multiple of block_size.")

    # Calculate the number of blocks
    num_blocks = size // block_size

    # Generate random values for each block
    block_values = np.random.normal(loc=0, scale=std, size=num_blocks)

    # Repeat block values to create the final array
    return np.repeat(block_values, block_size)

def time_hark(time, arr2D, cad):
    t2 = []
    a2D2 = []
    t=time[0]
    while t<= time[-1]:
        i = np.where(time >= t)[0][0]
        t2.append(time[i])
        a2D2.append(arr2D[i])
        t += cad
    return np.array(t2), np.array(a2D2)       
    
def degenerate_contrast_as_Harps(name_of_line, data, quiet_sun_subtraction=True, area_factor=60**2/np.pi/959.63**2,
                                 normal=True, add_noise=False, scale_noise=1):
    wav, contr, time, line, std = contrast_FD_data(name_of_line, data, quiet_sun_subtraction, area_factor=area_factor)
    
    time, contr = time_hark(time, contr, cad=4.5)
    
    if add_noise:
        noise = noise_alike_harps(name_of_line, contr.shape, backup=std[0], scale_noise=scale_noise)
        contr += noise

    return wav, smooth(contr), time, line, (std if add_noise else None)
    
def get_Harps(name_of_line):
    # Harps starts at 8:57 and end at 14:33 UT
    # first flare starts exactly at 8:57 (X2.2) the second at 11:53 (X9.3)
    folder = "E:/solar flares/data/2017-09-06/Harps/"
    flare = np.load(f'{folder}Flux_corrected.npy')[:63]

    wav = np.load(f'{folder}Wavelength.npy')

    summary = np.loadtxt('E:\solar flares\data\\2017-09-06\Harps\dates.txt', dtype=str) 
    time = np.array([un.hulp_time(t[11:123]) for t in summary])

    quiet_min = [43,58]
    
    if '9u' in name_of_line :
        quiet_min = [50, 55]
        time -= un.hulp_time("08:57:00")

    else:
        quiet_min = [43,55]
        time -= un.hulp_time("11:53:00")
    
    quiet_indeces = [np.where(time>= quiet_min[0])[0][0], np.where(time>= quiet_min[1])[0][0]]
    assert(quiet_indeces[0] < quiet_indeces[1])

    timeavg = np.median(flare[quiet_indeces[0]:quiet_indeces[1]], axis=0)
    timeavg = np.median(flare[10:30], axis=0)
    
    if '9u' in name_of_line :
        time = time[:5]  
        flarerange = flare[:5]
    else:
        time = time[30:63] 
        flarerange = flare[30:63]
        
    flarerange = (flarerange/timeavg) 

    if 'Ha' in name_of_line : 
        #Pick out window
        linecore = 6562.8
        flare_win = flarerange[:,265800:266800]
        wav_win = wav[265800:266800]

    elif 'CaK' in name_of_line : #Cak
        #Pick out window
        linecore = 3933.66
        flare_win = flarerange[:,2870:3870]
        wav_win = wav[2870:3870]

    elif 'CaH' in name_of_line : #CaH
        #Pick out window
        linecore = 3968.47
        flare_win = flarerange[:,6300:7300]
        wav_win = wav[6300:7300]
    else:
        print(f"Line could not be determined. only Ha, CaK and CaH possible. Got {name_of_line}.")
        return None, None, None

    return smooth(flare_win), wav_win, time



def create_gaussian_kernel(n, m, sigma):
    """
  Creates a 2D Gaussian distribution over an n x m matrix.

  Args:
      n: Number of rows in the matrix.
      m: Number of columns in the matrix.
      sigma: Standard deviation of the Gaussian distribution.

  Returns:
      A numpy array representing the Gaussian distribution.
  """
    x, y = np.meshgrid(np.arange(m), np.arange(n))
    mu_x = m // 2
    mu_y = n // 2
    return np.exp(-((x - mu_x) ** 2 + (y - mu_y) ** 2) / (2 * sigma**2))

from scipy.signal import savgol_filter

def smooth(data):
    if np.any(np.isnan(data)):
        raise ValueError("Data contains NaN. This is incompatible with the filter!")
    windowlength = min(101, np.shape(data)[1]//10)
    return savgol_filter(data, windowlength, 3)

def smooth2(data, n_wav=500, n_time=1, mode='same'):
    """
    Smooths a NumPy array by averaging over n neighbors.

    Args:
        data: A NumPy array of any shape.
        n: The number of neighbors to average over.

    Returns:
        A new NumPy array with the same shape as the input data,
        containing the smoothed values.
    """
    sigma = 0.4*n_wav
    kernel = create_gaussian_kernel(n_time, n_wav, sigma)
    kernel /= np.sum(kernel)
    sum_kernel = np.sum(kernel)
    return convolve2d(data, kernel, mode=mode)

def ax_contrastplot(fig, ax, X, Y, Z, x, line, decorations={}, seperate_colorbar=True, vlim=None, 
                    vlimscale=1, logscale=False, xlim=None, cmap='RdBu_r', lambda_0=None, non_centered=False):
    if vlim is None:   
        if non_centered:
            vmax = np.percentile(Z, 97)
            vmin = np.percentile(Z, 3)
        else:
            vmax = max(2-np.percentile(Z, 9), np.percentile(Z, 91))
            vmin =2-vmax
            print(f"Centerd contrast plot {vmin=}, {vmax=}.")
        # vmin = np.percentile(Z, 3)
        # vmax = np.percentile(Z, 97)
    else:
        vmin = 1-vlim
        vmax = vlim+1
        
    if logscale:
            pcm = ax.pcolormesh(X, Y, Z, cmap=cmap, norm=mplcolors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                            vmin=vmin, vmax=vmax, base=10), shading='auto')
    else:
        # norm = mplcolors.Normalize(vmin=vmin, vmax=vmax)
        pcm = ax.pcolormesh(X, Y, np.clip(Z, vmin, vmax), vmax=vmax, vmin=vmin, cmap=cmap, shading='auto', label=f'$\Delta I / \sigma$ []') #norm=norm, vmin=vmin, vmax=vmax,
        # pcm = ax.imshow(np.clip(Z, vmin, vmax), aspect="auto", cmap=cmap, origin='lower', 
        #                 extent=(np.min(X), np.max(X), np.min(Y), np.max(Y)), vmax=vmax, vmin=vmin) #
    if seperate_colorbar:
        cb = fig.colorbar(pcm, ax=ax, extend='both')#, format=cbformat)
        cb.ax.yaxis.set_offset_position('left')  

        # colorbar = plt.colorbar(pcm, label=f'$\Delta I / \sigma$ []')
        # colorbar.set_label('Z-Values')  # Replace with your desired label

    if 'title' in decorations:
        ax.set_title(decorations['title'], y=1.10)
    if "ylabel" in decorations:
        ax.set_ylabel(decorations['ylabel'])
    if "xlabel" in decorations:
        ax.set_ylabel(decorations['xlabel'])
        
    if xlim is not None:
        ax.set_xlim(xlim)

    ax2 = ax.twinx()  # instantiate a second axes that shares the same x-axis

    color = decorations['color'] if 'color' in decorations else 'black'
    # ax2.set_ylabel('Intensity []', color=color)  # we already handled the x-label with ax1
    ax2.plot(x,line, color=color)
    ax2.tick_params(axis='y', left=False, right=False, which="both", labelright=False)
    ax3 = ax.secondary_xaxis('top', functions=(wav_2_doppler(lambda_0), doppler_2_wav(lambda_0)))

    return pcm

def wav_2_doppler(lambda_0):
    c = 299792.458 # speed of light in km/s
    return lambda wavelengths :  c * ((wavelengths-lambda_0) / lambda_0) # Correct Doppler formula

def doppler_2_wav(lambda_0):
    c = 299792.458 # speed of light in km/s
    return lambda doppler_shifts : ((doppler_shifts * lambda_0) / c) + lambda_0

def acx_coord(ax, row, col):
    return ax[row+col] if ax.ndim == 1 else ax[row, col]

def Create_flare_contrast_plots(names_of_lines_list, data, quiet_sun_subtraction_list, long_names, Harps=False, 
                                scale_up=False, title='Contrast profiles', start_of_flare='', shared_colors_row=False, normal=True, vlims=None, add_noise=True, scale=10):
    
    rows = 2 + Harps + scale_up
    cols = len(names_of_lines_list)

    # make figure
    fig, ax = plt.subplots(rows, cols, figsize=(5*cols,4*rows), constrained_layout=True)
    fig.suptitle(title, fontsize=20)
    fig.supylabel(f"Minutes from start of flare {start_of_flare}")
    fig.supxlabel(r"Wavelength [$\AA$]")

    # row 1: sst FOV contrastplots
    for i, name in enumerate(names_of_lines_list):
        print("FOV: Line", name)
        wav, DFOV , time, line, std = contrast_FOV_data(name, data, quiet_sun_subtraction_list[i], normal=normal)
        W, T = np.meshgrid(wav, time)

        decorations={"title":long_names[name]}
        if i==0:
            decorations["ylabel"] = "FOV contrast profile"

        ax_contrastplot(fig, acx_coord(ax, 0, i), W, T, DFOV, wav, line, decorations , seperate_colorbar=True, logscale=False)

    # row 2: nessi + FOV = Full disk contrastplots
    for i, name in enumerate(names_of_lines_list):
        vlim = vlims[name] if vlims is not None else None
        print("FD: Line", name, vlim)
        wav, DFD , time, line, std = contrast_FD_data(name, data, quiet_sun_subtraction_list[i], normal=normal, add_noise=add_noise)
        W, T = np.meshgrid(wav, time)

        decorations={}
        if i==0:
            decorations["ylabel"] = "Full disk (Nessi + flare) contrast profile"

        ax_contrastplot(fig, acx_coord(ax, 1, i), W, T, DFD, wav, line, decorations, seperate_colorbar=True, vlim=vlim)

    # row 3: scaled up flare
    if scale_up:
        for i, name in enumerate(names_of_lines_list):
            
            wav, DFD , time, line, std = contrast_FD_data(name,data, quiet_sun_subtraction_list[i], area_factor=area_factor*scale, normal=normal, add_noise=add_noise)
            W, T = np.meshgrid(wav, time)

            decorations={}
            if i==0:
                decorations["ylabel"] = f"{scale}x Scaled flare profile"
            vlim = vlims[name] if vlims is not None else None
            print("scaled flare: Line", name, vlim)
            ax_contrastplot(fig, acx_coord(ax, 2, i), W, T, DFD, wav, line, decorations, seperate_colorbar=True, vlim=vlim)


    # row 4: Harps data
    if Harps:
        acx_coord(ax, 3, 0).set_ylabel("Harps contrast profile")
        for i, name in enumerate(names_of_lines_list):
            
            flare_win, wav_win, time = get_Harps(name)
            if flare_win is not None:
                wav_nessi, dc_nessi, clv_nessi = data[f"nessi_{name}"]
                if "CaK" in name:
                    wav_nessi, dc_nessi = wav_nessi[:-1], dc_nessi[:-1]
                W, T = np.meshgrid(wav_win, time)

                decorations={}
                if i==0:
                    decorations["ylabel"] = "Harps contrast profile"
                vlim = vlims[name] if vlims is not None else None
                print("Harps: Line", name, vlim)
                ax_contrastplot(fig, acx_coord(ax, 3, i), W, T, flare_win, wav_nessi, dc_nessi, decorations , vlimscale=2/3,
                                seperate_colorbar=True, xlim=(wav_nessi[0], wav_nessi[-1]), vlim=vlim)

    figname="contrastplot"
    for name in names_of_lines_list:
        figname += name
    fig.savefig(f'D:\solar flares\data\plots\{figname}.png')
  



def Harps_contrast_plots(names_of_lines_list, data, quiet_sun_subtraction_list, long_names, Harps=False, 
                                scale_up=False, title='Contrast profiles Harps', start_of_flare='', shared_colors_row=False):
    
    rows = 1
    cols = len(names_of_lines_list)

    # make figure
    fig, ax = plt.subplots(rows, cols, figsize=(5*cols,4*rows), constrained_layout=True)
    fig.suptitle(title, fontsize=20)
    fig.supylabel(f"Minutes from start of flare {start_of_flare}")
    fig.supxlabel(r"Wavelength [$\AA$]")

    # row 1: sst Harps contrastplots
    if Harps:
        for i, name in enumerate(names_of_lines_list):
            print("Harps: Line", name)
            flare_win, wav_win, time = get_Harps(name)
            if flare_win is not None:
                wav_nessi, dc_nessi, clv_nessi = data[f"nessi_{name}"]
                if "CaK" in name:
                    wav_nessi, dc_nessi = wav_nessi[:-1], dc_nessi[:-1]
                W, T = np.meshgrid(wav_win, time)

                decorations={"title":long_names[name]}
                if i==0:
                    decorations["ylabel"] = "Harps contrast profile"

                ax_contrastplot(fig, ax[i], W, T, flare_win, wav_nessi, dc_nessi, decorations , 
                                seperate_colorbar=True, vlimscale=2/3, logscale=False, xlim=(wav_nessi[0], wav_nessi[-1]))
                


import matplotlib.patches as patches

# Function to plot rectangles representing the field of view
def plot_fov(ax, x, y, width, height, color='red', label='', angle=0):
    rect = patches.Rectangle((x - width / 2, y - height / 2), width, height,
                             linewidth=1, edgecolor=color, facecolor='none', label=label, angle=angle)
    ax.add_patch(rect)
    ax.text(x, y+height, label, ha='center', va='center', fontsize=8)

def plot_centrical_circ(ax, mu, x, y):
    r = (1-mu**2)**0.5
    theta = np.arange(0, 2*np.pi, 0.01)
    circx = r*np.cos(theta) + x
    circy = r*np.sin(theta) + y
    ax.plot(circx, circy, linestyle='--', alpha=0.9, color='gray', linewidth=0.7)
    if mu>=0.4 or mu==0:
        ax.text(0, r, fr'$\mu$={mu:.1f}', ha='center', va='bottom', fontsize=8, color='gray')

def plot_mu_grid(ax, mu):
    for i in mu:
        plot_centrical_circ(ax, i, 0, 0)

# Function to create the sun plot
def create_sun_plot(centers, flarelabels, z=60, sr=959.63, angels=None):
    fig, ax = plt.subplots(figsize=(5,5))

    # Plot the yellow disk for the sun
    sun = plt.Circle((0, 0), 1, color='yellow', ec='black')
    ax.add_patch(sun)

    # Plot rectangles for the field of view at various locations
    for i, center in enumerate(centers):
        angle = angels[i] if angels is not None else 0
        plot_fov(ax, x=center[0]/sr, y=center[1]/sr, width=z/sr, height=z/sr, label=flarelabels[i], angle=angle)

    # plot mu grid
    plot_mu_grid(ax, mu = np.arange(0, 1.05, 0.1))

    # Customize the plot
    ax.set_aspect('equal')  # Ensure aspect ratio is equal
    ax.set_title('Sun with Field of View of observations', fontsize=15)
    ax.axis('off')  # Hide the axes

    plt.show()
