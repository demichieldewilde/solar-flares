##All spec

global spec2d
spec2d = 0


def dopplerck(wavelengths):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 3934 # central wavelength in Angstrom
    if spec2d:
        doppler_shifts = c * ((wavelengths) / lambda_0)
    else:
        doppler_shifts = c * ((wavelengths-lambda_0) / lambda_0)
    return doppler_shifts

def idopplerck(doppler_shifts):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 3934 # central wavelength in Angstrom
    if spec2d:
        wavelengths = (lambda_0 * doppler_shifts) / c
    else:
        wavelengths = (lambda_0 * doppler_shifts) / c - lambda_0
    return wavelengths

def dopplerha(wavelengths):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 6563 # central wavelength in Angstrom
    if spec2d:
        doppler_shifts = c * ((wavelengths) / lambda_0)
    else:
        doppler_shifts = c * ((wavelengths-lambda_0) / lambda_0)
    return doppler_shifts

def idopplerha(doppler_shifts):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 6563 # central wavelength in Angstrom
    if spec2d:
        wavelengths = (lambda_0 * doppler_shifts) / c
    else:
        wavelengths = (lambda_0 * doppler_shifts) / c - lambda_0
    return wavelengths

def dopplerc8(wavelengths):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 8542 # central wavelength in Angstrom
    if spec2d:
        doppler_shifts = c * ((wavelengths) / lambda_0)
    else:
        doppler_shifts = c * ((wavelengths-lambda_0) / lambda_0)
    return doppler_shifts

def idopplerc8(doppler_shifts):
    c = 299792.458 # speed of light in km/s
    lambda_0 = 8542 # central wavelength in Angstrom
    if spec2d:
        wavelengths = (lambda_0 * doppler_shifts) / c
    else:
        wavelengths = (lambda_0 * doppler_shifts) / c - lambda_0
    return wavelengths

#todo
#check if red bit of arcade leads to shifting out
#f70 follow fibril arcade thing
#f103, 206 follow arcade


import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mpl



#plot lavenders
fig = plt.figure(figsize=(10, 7))


plt.subplot(351)
c = np.arange(25)
wavck = np.arange(-1.05,1.05,0.1)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Reds)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavck+3934,(reg1ik[0:])[i], c=cmap.to_rgba(i + 1))
#plt.colorbar(cmap)
plt.xticks([3933,3934,3935])
#plt.ylabel('Intensity')
plt.yticks([])

ax1 = plt.gca()  # Get the current axis (i.e., the one just created)
ax1a = ax1.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax1a.set_xticks([-50,0,50])
ax1.tick_params(axis='x', labelrotation=45)

plt.subplot(352)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavck+3934,(reg3ik[6:])[i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap)
plt.xticks([3933,3934,3935])
#plt.colorbar(cmap, ticks=c)
plt.yticks([])

ax2 = plt.gca()  # Get the current axis (i.e., the one just created)
ax2a = ax2.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax2a.set_xticks([-50,0,50])
ax2.tick_params(axis='x', labelrotation=45)

plt.subplot(353)
c=int(327)
norm = mpl.colors.Normalize(vmin=-50, vmax=c)
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.gray_r)
wwck = wavelck[:-1]-3933.66+3934
for i in range(50,327,5):
    plt.plot(wwck,selected_pixels_ck[:-1,i], c=cmap.to_rgba(i + 1))
c=int(162)
plt.xticks([3933,3934,3935])
plt.yticks([])


ax3 = plt.gca()  # Get the current axis (i.e., the one just created)
ax3a = ax3.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax3a.set_xticks([-50,0,50])
ax3a.set_xlabel('V [km/s]')
ax3.tick_params(axis='x', labelrotation=45)

plt.subplot(354)
c = np.arange(25)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavck+3934,(reg5ik[:])[i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap, ticks=c)
plt.xticks([3933,3934,3935])
#plt.colorbar(cmap)
plt.yticks([])

ax4 = plt.gca()  # Get the current axis (i.e., the one just created)
ax4a = ax4.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax4a.set_xticks([-50,0,50])
ax4.tick_params(axis='x', labelrotation=45)

plt.subplot(355)
c = np.arange(25)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Purples)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavck+3934,regwick[0:][i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap, ticks=c)
plt.xticks([3933,3934,3935])
#plt.colorbar(cmap)
plt.yticks([])

linecore = 3933.66
ax5 = plt.gca()  # Get the current axis (i.e., the one just created)
ax5a = ax5.secondary_xaxis('top', functions=(dopplerck, idopplerck))
ax5a.set_xticks([-50,0,50])
ax5.tick_params(axis='x', labelrotation=45)

plt.subplot(356)
c = np.arange(25)
wavh = np.arange(-1.5,1.5,0.2)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Reds)
#cmap.set_array([])
plt.yticks([])
for i in range(20):
    plt.plot(wavh+6563,(reg1ih[0:])[i], c=cmap.to_rgba(i + 1))
#plt.colorbar(cmap)
plt.xticks([6562,6563,6564])
plt.ylabel('Intensity')

ax6 = plt.gca()  # Get the current axis (i.e., the one just created)
ax6a = ax6.secondary_xaxis('top', functions=(dopplerha, idopplerha))
ax6a.set_xticks([-50,0,50])
ax6.tick_params(axis='x', labelrotation=45)


plt.subplot(357)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavh+6563,(reg3ih[6:])[i], c=cmap.to_rgba(i + 1))
#plt.colorbar(cmap)
plt.xticks([6562,6563,6564])
#plt.colorbar(cmap, ticks=c)
plt.yticks([])

ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
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

listor = np.zeros([13,28])
ii = 0
for i in (np.arange(20,160,5)):

    listor[:,ii] = np.flip(selected_pixels_ha, axis=1)[:,i]
    ii+=1
listor = np.flip(listor, axis=1)

for i in np.arange(28):
    plt.plot(wwha,listor[:,i], alpha=1, c=cmap.to_rgba(i))
plt.yticks([])

ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerha, idopplerha))
ax7a.set_xticks([-50,0,50])
plt.xticks([6562,6563,6564])
ax7.tick_params(axis='x', labelrotation=45)

plt.subplot(359)
c = np.arange(25)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavh+6563,(reg5ih[10:])[i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap, ticks=c)
#plt.colorbar(cmap)
plt.yticks([])
plt.xticks([6562,6563,6564])

ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerha, idopplerha))
ax7a.set_xticks([-50,0,50])
ax7.tick_params(axis='x', labelrotation=45)

plt.subplot(3,5,10)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Purples)
#cmap.set_array([])
for i in range(20):
    plt.plot(wavh+6563,regwiha[:][i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap, ticks=c)
#plt.colorbar(cmap)
plt.yticks([])
plt.xticks([6562,6563,6564])

linecore = 6562.8
ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerha, idopplerha))
ax7a.set_xticks([-50,0,50])
ax7.tick_params(axis='x', labelrotation=45)


plt.subplot(3,5,11)
c = np.arange(20)
wav8 = np.arange(-0.7,0.8,0.1)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Reds)
#cmap.set_array([])

for i in range(15):
    plt.plot(wav8+8542,(reg1ic[0:])[i], c=cmap.to_rgba(i + 1))
#plt.colorbar(cmap)
#plt.ylabel('Intensity')
#plt.xlabel(r'$\lambda$ [$\rm{\AA}$]')
plt.yticks([])
plt.xticks([8541.5,8542,8542.5])

linecore = 8542
ax11 = plt.gca()  # Get the current axis (i.e., the one just created)
ax11.tick_params(axis='x', labelrotation=45)
ax11a = ax11.secondary_xaxis('top', functions=(dopplerc8, idopplerc8))
ax11a.set_xticks([-25,0,25])
ax11.tick_params(axis='x', labelrotation=45)

plt.subplot(3,5,12)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Oranges)
#cmap.set_array([])
plt.xticks([8541.5,8542,8542.5])
for i in range(15):
    plt.plot(wav8+8542,(reg3ic[10:])[i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap)
#plt.xlabel(r'$\lambda$ [$\rm{\AA}$]')
plt.yticks([])

linecore = 8542
ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7.tick_params(axis='x', labelrotation=45)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerc8, idopplerc8))
ax7a.set_xticks([-25,0,25])

plt.subplot(3,5,13)
c=28#int(175)
norm = mpl.colors.Normalize(vmin=-50, vmax=c)
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.gray_r)

listor2 = np.zeros([11,28])
ii = 0
for i in (np.arange(20,160,5)):

    listor2[:,ii] = selected_pixels[:,i]
    ii+=1
#listor2 = np.flip(listor2, axis=1)

for i in np.arange(28):
    plt.plot(wavel8542,listor2[:,i], alpha=1, c=cmap.to_rgba(i))
plt.yticks([])


#for i in range(0,150,5):
#    plt.plot(wavel8542,selected_pixels[:,i], c=cmap.to_rgba(i + 1))
#plt.yticks([])
plt.xlabel(r'$\lambda$ [$\rm{\AA}$]')

linecore = 8542
ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7.tick_params(axis='x', labelrotation=45)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerc8, idopplerc8))
ax7a.set_xticks([-25,0,25])

#plt.colorbar(cmap, ticks=c)
c = np.arange(20)
plt.subplot(3,5,14)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Blues)
#cmap.set_array([])
for i in range(15):
    plt.plot(wav8+8542,(reg5ic[7:])[i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap, ticks=c)
#plt.xlabel(r'$\lambda$ [$\rm{\AA}$]')
plt.yticks([])
plt.xticks([8541.5,8542,8542.5])

linecore = 8542
ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7.tick_params(axis='x', labelrotation=45)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerc8, idopplerc8))
ax7a.set_xticks([-25,0,25])
#plt.colorbar(cmap)

plt.subplot(3,5,15)
norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.Purples)
#cmap.set_array([])
for i in range(15):
    plt.plot(wav8+8542,regwic8[5:][i], c=cmap.to_rgba(i + 1))
#fig.colorbar(cmap, ticks=c)
#plt.xlabel(r'$\lambda$ [$\rm{\AA}$]')
plt.yticks([])
plt.xticks([8541.5,8542,8542.5])
#plt.colorbar(cmap)

linecore = 8542
ax7 = plt.gca()  # Get the current axis (i.e., the one just created)
ax7.tick_params(axis='x', labelrotation=45)
ax7a = ax7.secondary_xaxis('top', functions=(dopplerc8, idopplerc8))
ax7a.set_xticks([-25,0,25])

plt.tight_layout()
plt.savefig('spectraall.pdf')
plt.show()
### Diagonal