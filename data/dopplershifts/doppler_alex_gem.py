##All spec

spec2d = 0 # global spec2d is now initialized here, in global scope.

def dopplerck(wavelengths):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 3934 # central wavelength in Angstrom
    if spec2d:
        doppler_shifts = c * ((wavelengths) / lambda_0)
    else:
        doppler_shifts = c * ((wavelengths-lambda_0) / lambda_0) # Correct Doppler formula
    return doppler_shifts

def idopplerck(doppler_shifts):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 3934 # central wavelength in Angstrom
    if spec2d:
        wavelengths = (lambda_0 * doppler_shifts) / c
    else:
        wavelengths = ((doppler_shifts * lambda_0) / c) + lambda_0 # Correct inverse Doppler formula
    return wavelengths

def dopplerha(wavelengths):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 6563 # central wavelength in Angstrom
    if spec2d:
        doppler_shifts = c * ((wavelengths) / lambda_0)
    else:
        doppler_shifts = c * ((wavelengths-lambda_0) / lambda_0) # Correct Doppler formula
    return doppler_shifts

def idopplerha(doppler_shifts):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 6563 # central wavelength in Angstrom
    if spec2d:
        wavelengths = (lambda_0 * doppler_shifts) / c
    else:
        wavelengths = ((doppler_shifts * lambda_0) / c) + lambda_0 # Correct inverse Doppler formula
    return wavelengths

def dopplerc8(wavelengths):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 8542 # central wavelength in Angstrom
    if spec2d:
        doppler_shifts = c * ((wavelengths) / lambda_0)
    else:
        doppler_shifts = c * ((wavelengths-lambda_0) / lambda_0) # Correct Doppler formula
    return doppler_shifts

def idopplerc8(doppler_shifts):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 8542 # central wavelength in Angstrom
    if spec2d:
        wavelengths = (lambda_0 * doppler_shifts) / c
    else:
        wavelengths = ((doppler_shifts * lambda_0) / c) + lambda_0 # Correct inverse Doppler formula
    return wavelengths


import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


#plot lavenders
fig = plt.figure(figsize=(10, 7))


plt.subplot(351)
c = np.arange(25)
wavck = np.arange(-1.05,1.05,0.1)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Reds)

# Create dummy data for the first row (Ca II K)
num_spectra_ck = 20 # Based on loop range
reg1ik = np.random.rand(num_spectra_ck, len(wavck)) + 1 # Intensity data, add 1 to make it positive
reg3ik = np.random.rand(num_spectra_ck, len(wavck)) + 1
reg5ik = np.random.rand(num_spectra_ck, len(wavck)) + 1
regwick = np.random.rand(num_spectra_ck, len(wavck)) + 1
wavelck = wavck + 3934 # dummy wavelength array centered at 3934
selected_pixels_ck = np.random.rand(len(wavelck) -1, 327) # Shape inferred from plotting
wwck = wavelck[:-1]-3933.66+3934 # dummy wavelenghts for gray subplot


for i in range(num_spectra_ck):
    plt.plot(wavck+3934,reg1ik[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.xticks([3933,3934,3935])
plt.yticks([])

ax1 = plt.gca()  # Get the current axis
ax1a = ax1.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax1a.set_xticks([-50,0,50])
ax1.tick_params(axis='x', labelrotation=45)

plt.subplot(352)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)

for i in range(num_spectra_ck):
    plt.plot(wavck+3934,reg3ik[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.xticks([3933,3934,3935])
plt.yticks([])

ax2 = plt.gca()
ax2a = ax2.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax2a.set_xticks([-50,0,50])
ax2.tick_params(axis='x', labelrotation=45)

plt.subplot(353)
c=int(327)
norm = mpl.colors.Normalize(vmin=-50, vmax=c)
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.gray_r)
wwck = wavelck[:-1]-3933.66+3934
for i in range(50,327,5):
    plt.plot(wwck,selected_pixels_ck[:,i], c=cmap.to_rgba(i + 1))
c=int(162)
plt.xticks([3933,3934,3935])
plt.yticks([])

ax3 = plt.gca()
ax3a = ax3.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax3a.set_xticks([-50,0,50])
ax3a.set_xlabel('V [km/s]')
ax3.tick_params(axis='x', labelrotation=45)

plt.subplot(354)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
for i in range(num_spectra_ck):
    plt.plot(wavck+3934,reg5ik[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.xticks([3933,3934,3935])
plt.yticks([])

ax4 = plt.gca()
ax4a = ax4.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax4a.set_xticks([-50,0,50])
ax4.tick_params(axis='x', labelrotation=45)

plt.subplot(355)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Purples)
for i in range(num_spectra_ck):
    plt.plot(wavck+3934,regwick[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.xticks([3933,3934,3935])
plt.yticks([])

linecore = 3933.66
ax5 = plt.gca()
ax5a = ax5.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax5a.set_xticks([-50,0,50])
ax5.tick_params(axis='x', labelrotation=45)

plt.subplot(356)
c = np.arange(25)
wavh = np.arange(-1.5,1.5,0.2)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Reds)

# Create dummy data for the second row (H-alpha)
num_spectra_ha = 20 # Based on loop range
reg1ih = np.random.rand(num_spectra_ha, len(wavh)) + 1
reg3ih = np.random.rand(num_spectra_ha, len(wavh)) + 1
reg5ih = np.random.rand(num_spectra_ha, len(wavh)) + 1
regwiha = np.random.rand(num_spectra_ha, len(wavh)) + 1
wavelha = wavh + 6563 # dummy wavelength array centered at 6563
selected_pixels_ha = np.random.rand(13, 160) # Shape inferred from plotting


plt.yticks([])
for i in range(num_spectra_ha):
    plt.plot(wavh+6563,reg1ih[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.xticks([6562,6563,6564])
plt.ylabel('Intensity') # Only Y label is on the first subplot of second row

ax6 = plt.gca()
ax6a = ax6.secondary_xaxis('top', functions=(dopplerha, idopplerha))
ax6a.set_xticks([-50,0,50])
ax6.tick_params(axis='x', labelrotation=45)


plt.subplot(357)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)
for i in range(num_spectra_ha):
    plt.plot(wavh+6563,reg3ih[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.xticks([6562,6563,6564])
plt.yticks([])

ax7 = plt.gca()
ax7a = ax7.secondary_xaxis('top', functions=(dopplerha, idopplerha))
ax7a.set_xticks([-50,0,50])
ax7.tick_params(axis='x', labelrotation=45)
ax7.tick_params(axis='x', labelrotation=45)


#HIEEEEEEEEERRRRRRRRRRRR
plt.subplot(358)
c=28#int(175)
norm = mpl.colors.Normalize(vmin=-50, vmax=c)
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.gray_r)
wwha = wavelha +0.2

listor = np.zeros([15,28])
ii = 0
for i in (np.arange(20,160,5)):

    listor[:-2,ii] = np.flip(selected_pixels_ha, axis=1)[:,i]
    ii+=1
listor = np.flip(listor, axis=1)

for i in np.arange(28):
    plt.plot(wwha,listor[:,i], alpha=1, c=cmap.to_rgba(i))
plt.yticks([])

ax8 = plt.gca() # Fixed variable name ax7->ax8 to avoid overwrite
ax8a = ax8.secondary_xaxis('top', functions=(dopplerha, idopplerha)) # Fixed variable name ax7->ax8
ax8a.set_xticks([-50,0,50])
plt.xticks([6562,6563,6564])
ax8.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax8

plt.subplot(359)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
for i in range(num_spectra_ha):
    plt.plot(wavh+6563,reg5ih[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.yticks([])
plt.xticks([6562,6563,6564])

ax9 = plt.gca() # Fixed variable name ax7->ax9
ax9a = ax9.secondary_xaxis('top', functions=(dopplerha, idopplerha)) # Fixed variable name ax7->ax9
ax9a.set_xticks([-50,0,50])
ax9.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax9

plt.subplot(3,5,10)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Purples)
for i in range(num_spectra_ha):
    plt.plot(wavh+6563,regwiha[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.yticks([])
plt.xticks([6562,6563,6564])

linecore = 6562.8
ax10 = plt.gca() # Fixed variable name ax7->ax10
ax10a = ax10.secondary_xaxis('top', functions=(dopplerha, idopplerha)) # Fixed variable name ax7->ax10
ax10a.set_xticks([-50,0,50])
ax10.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax10


plt.subplot(3,5,11)
c = np.arange(20)
wav8 = np.arange(-0.7,0.8,0.1)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Reds)

# Create dummy data for the third row (Ca II IR)
num_spectra_c8 = 15 # Based on loop range
reg1ic = np.random.rand(num_spectra_c8, len(wav8)) + 1
reg3ic = np.random.rand(num_spectra_c8, len(wav8)) + 1
reg5ic = np.random.rand(num_spectra_c8, len(wav8)) + 1
regwic8 = np.random.rand(num_spectra_c8, len(wav8)) + 1
wavel8542 = wav8 + 8542 # dummy wavelength array centered at 8542
selected_pixels = np.random.rand(11, 160) # Shape inferred from plotting

for i in range(num_spectra_c8):
    plt.plot(wav8+8542,reg1ic[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.yticks([])
plt.xticks([8541.5,8542,8542.5])

linecore = 8542
ax11 = plt.gca()
ax11.tick_params(axis='x', labelrotation=45)
ax11a = ax11.secondary_xaxis('top', functions=(dopplerc8, idopplerc8))
ax11a.set_xticks([-25,0,25]) # Doppler range is smaller for Ca II IR
ax11.tick_params(axis='x', labelrotation=45)

plt.subplot(3,5,12)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)
plt.xticks([8541.5,8542,8542.5])
for i in range(num_spectra_c8):
    plt.plot(wav8+8542,reg3ic[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.yticks([])

linecore = 8542
ax12 = plt.gca() # Fixed variable name ax7->ax12
ax12.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax12
ax12a = ax12.secondary_xaxis('top', functions=(dopplerc8, idopplerc8)) # Fixed variable name ax7->ax12
ax12a.set_xticks([-25,0,25])

plt.subplot(3,5,13)
c=28#int(175)
norm = mpl.colors.Normalize(vmin=-50, vmax=c)
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.gray_r)

listor2 = np.zeros([15,28])
ii = 0
for i in (np.arange(20,160,5)):

    listor2[:-4,ii] = selected_pixels[:,i]
    ii+=1

for i in np.arange(28):
    plt.plot(wavel8542,listor2[:,i], alpha=1, c=cmap.to_rgba(i))
plt.yticks([])
plt.xlabel(r'$\lambda$ [$\rm{\AA}$]') # Only X label is on the first subplot of third row

linecore = 8542
ax13 = plt.gca() # Fixed variable name ax7->ax13
ax13.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax13
ax13a = ax13a = ax13.secondary_xaxis('top', functions=(dopplerc8, idopplerc8)) # Fixed variable name ax7->ax13
ax13a.set_xticks([-25,0,25])

c = np.arange(20)
plt.subplot(3,5,14)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
for i in range(num_spectra_c8):
    plt.plot(wav8+8542,reg5ic[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.yticks([])
plt.xticks([8541.5,8542,8542.5])

linecore = 8542
ax14 = plt.gca() # Fixed variable name ax7->ax14
ax14.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax14
ax14a = ax14.secondary_xaxis('top', functions=(dopplerc8, idopplerc8)) # Fixed variable name ax7->ax14
ax14a.set_xticks([-25,0,25])

plt.subplot(3,5,15)
norm = mpl.colors.Normalize(vmin=np.min(c), vmax=np.max(c))
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Purples)
for i in range(num_spectra_c8):
    plt.plot(wav8+8542,regwic8[i], c=cmap.to_rgba(i + 1)) # Corrected indexing
plt.yticks([])
plt.xticks([8541.5,8542,8542.5])

linecore = 8542
ax15 = plt.gca() # Fixed variable name ax7->ax15
ax15.tick_params(axis='x', labelrotation=45) # Fixed variable name ax7->ax15
ax15a = ax15.secondary_xaxis('top', functions=(dopplerc8, idopplerc8)) # Fixed variable name ax7->ax15
ax15a.set_xticks([-25,0,25])

plt.tight_layout()
plt.savefig('spectraall.pdf')
plt.show()
### Diagonal