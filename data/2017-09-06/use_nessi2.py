import numpy as np
import matplotlib.pyplot as plt
import use_nessi as un
import os
from importlib import reload
from matplotlib import cm
import matplotlib.cbook as cbook
import matplotlib.colors as colors
from scipy.interpolate import interp1d
from scipy.signal import convolve2d


sr = solar_radius = 959.63
area_factor = 60**2/np.pi/sr**2


def disgard_cont_point(name, data):
    # FOV_spectrum
    data[f'FOV_{name}'] = data[f'FOV_{name}'][:, :-1]

    # load quiet_sun profile
    data[f'quiet_sun_{name}'] = data[f'quiet_sun_{name}'][:, :-1]


sr = solar_radius = 959.63
area_factor = 60**2/np.pi/sr**2

'''
The scale_pix_to_saas is the number by which the standard deviation has to be multiplied to get to the the std of the saas observation
it should be the inverse of square root of (number of pixels in sst (1920×1200 CHROMIS or 1k × 1k CRISP) * areafactor (60**2/np.pi/sr**2) ) 
hence 1/53.544360373477076
'''
def contrast_FOV_data(name_of_line, data, quiet_sun_subtraction=True, num=100, normal=True, scale_pix_to_saas=1/53.4):
    FOV = data[f"FOV_{name_of_line}"]
    wav_qs, qs_spec, std_qs = data[f"quiet_sun_{name_of_line}"]
    time = data[f"TIME_{name_of_line}"]
    wav_nessi, dc_nessi, clv_nessi = data[f"nessi_{name_of_line}"]
    std_qs *= scale_pix_to_saas

    if quiet_sun_subtraction:
        qs_spectc =  qs_spec
        wav_qsc = wav_qs
    else:
        qs_spectc = dc_nessi * clv_nessi
        wav_qsc = wav_nessi

    wav = smooth_wavelengths(wav_qs, wav_nessi, num)
    
    if normal:  
        DFOV = np.array([interp1d(wav_qs, FOV[i,:])(wav)  for i in range(np.shape(FOV)[0])]) - interp1d(wav_qsc, qs_spectc)(wav)
    else:
        DFOV = np.array([interp1d(wav_qs, FOV[i,:])(wav)  for i in range(np.shape(FOV)[0])]) / interp1d(wav_qsc, qs_spectc)(wav) -1

    line = interp1d(wav_qsc, qs_spectc)(wav)
    std = interp1d(wav_qs, std_qs)(wav)

    return wav, DFOV , time, line, std

def smooth_wavelengths(wav1, wav2, num=100):
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

def noise_alike_harps(name, shape, backup=1, x_steps=1):
    std = get_harps_std(name, backup)
    noise = np.array([generate_random_array(shape[1], block_size=x_steps, std=std)  for j in range(shape[0])]) #np.random.normal(loc=0, scale=std, size=shape)
    return noise

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
    
def degenerate_contrast_as_Harps(name_of_line, data, quiet_sun_subtraction=True, area_factor=60**2/np.pi/959.63**2,normal=True, add_noise=False):
    wav, DFOV, time, line, std = contrast_FOV_data(name_of_line, data, quiet_sun_subtraction, normal=normal)
    wav_nessi, dc_nessi, clv_nessi = data[f"nessi_{name_of_line}"]

    line = interp1d(wav_nessi, dc_nessi)(wav)
    DFD = area_factor * DFOV 
    
    time, DFD = time_hark(time, DFD, cad=5)
    
    if add_noise:
        noise = noise_alike_harps(name_of_line, DFD.shape, backup=std[0])
        DFD += noise

    return wav, smooth(DFD), time, line, (std if add_noise else None)
    

def contrast_FD_data(name_of_line, data, quiet_sun_subtraction=True, area_factor=60**2/np.pi/959.63**2,normal=True, add_noise=False):            
    wav, DFOV, time, line, std = contrast_FOV_data(name_of_line, data, quiet_sun_subtraction, normal=normal)
    wav_nessi, dc_nessi, clv_nessi = data[f"nessi_{name_of_line}"]

    line = interp1d(wav_nessi, dc_nessi)(wav)
    # Correct normalization for area and mu-value (all intensities are normalized on the first wavelength)
    #   -correction for area is area_factor
    #   -correction factor for average vs quiet sun normalization is 1/(dc_nessi[0]*clv_nessi[0])
    #   -correction factor for mu value is clv[0]/1
    # thereby we have area_factor * 1/(dc_nessi[0]*clv_nessi[0]) * clv[0]/1 thus
    DFD = area_factor * DFOV / dc_nessi[0]
    print('the correction for the mu_value normalization is ', 1/dc_nessi[0])
    
    if add_noise:
        noise = np.random.normal(loc=0, scale=std[0], size=(DFD.shape))
        DFD += noise

    return wav, DFD, time, line, (std if add_noise else None)
    
def get_Harps(name_of_line):
    # Harps starts at 8:57 and and at 14:33 UT
    # first flare starts exactly at 8:57 (X2.2) the second at 11:53 (X9.3)
    folder = "D:/solar flares/data/2017-09-06\Harps/"
    flare = np.load(f'{folder}Flux_corrected.npy')

    wav = np.load(f'{folder}Wavelength.npy')

    timeavg = np.median(flare[10:30], axis=0)

    flarerange = (flare[:63] / timeavg) - 1
    cadence =  5 # Harps cadence in minutes
    time = np.arange(np.shape(flarerange)[0]) * cadence

    line = 6563

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
        print(f"Line could not be detirmend. only Ha, CaK and CaH possible. Got {name_of_line}.")
        return None, None, None

    flare_win = flare_win[:5] if '9u' in name_of_line else flare_win[30:]
    time = time[:5] if '9u' in name_of_line else time[30:63]- (11*60+53) + (8*60+57)
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

def ax_contrastplot(fig, ax, X, Y, Z, x, line, decorations={}, seperate_colorbar=True, vlim=None, vlimscale=1, logscale=False, xlim=None):
    if vlim is None:   
        vlim = np.max(np.abs(Z)) * vlimscale
        
    if logscale:
            pcm = ax.pcolormesh(X, Y, Z, cmap='RdBu_r', norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03,
                                            vmin=-vlim, vmax=vlim, base=10), shading='auto')
    else:
        pcm = ax.pcolormesh(X, Y, Z, cmap='RdBu_r',vmin=-vlim, vmax=vlim, shading='auto', label=f'$\Delta I / \sigma$ []')

    if seperate_colorbar:
        # print('X', X, 'Y', Y, 'Z', Z)
        fig.colorbar(pcm, ax=ax, extend='both')
        # colorbar = plt.colorbar(pcm, label=f'$\Delta I / \sigma$ []')
        # colorbar.set_label('Z-Values')  # Replace with your desired label

    if 'title' in decorations:
        ax.set_title(decorations['title'])
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
    ax2.tick_params(axis='y', labelcolor=color)

    return pcm


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
def plot_fov(ax, x, y, width, height, color='red', label=''):
    rect = patches.Rectangle((x - width / 2, y - height / 2), width, height,
                             linewidth=1, edgecolor=color, facecolor='none', label=label)
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
def create_sun_plot(centers, flarelabels, z=60/959.63):
    fig, ax = plt.subplots(figsize=(5,5))

    # Plot the yellow disk for the sun
    sun = plt.Circle((0, 0), 1, color='yellow', ec='black')
    ax.add_patch(sun)

    # Plot rectangles for the field of view at various locations
    for i, center in enumerate(centers):
        plot_fov(ax, x=center[0], y=center[1], width=z, height=z, label=flarelabels[i])

    # plot mu grid
    plot_mu_grid(ax, mu = np.arange(0, 1.05, 0.1))

    # Customize the plot
    ax.set_aspect('equal')  # Ensure aspect ratio is equal
    ax.set_title('Sun with Field of View of observations', fontsize=15)
    ax.axis('off')  # Hide the axes

    plt.show()
