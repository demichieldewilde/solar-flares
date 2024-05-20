#!/usr/bin/env python3
# coding: utf-8
# https://docs.astropy.org/en/stable/generated/examples/io/skip_create-large-fits.html

import numpy as np

from astropy.io import fits


def create_nd_fits_file(fname, dims, dtype=np.float32):

  n = len(dims)
  fdims = tuple([1,]*n)
  fake_data = np.zeros(fdims, dtype=dtype)

  #
  # Allow header stuff:
  nblocks = 5#10


  #
  # Create fits header with the fake (allocatable and proper rank) data
  hdu = fits.PrimaryHDU(data=fake_data)
  header = hdu.header
  while (len(header) < (36 * nblocks - 1)):
      header.append()  # Adds a blank card to the end

  #
  # Modify header to accomodate the full array (reverse dims order as python's fastest axis is the latest):
  for i in range(header['NAXIS']):
    header['NAXIS%i' % (i+1,)] = dims[::-1][i]
  header.tofile(fname, overwrite=True)

  shape = tuple(header[f'NAXIS{ii}'] for ii in range(1, header['NAXIS']+1))
  with open(fname, 'rb+') as fobj:
    p = len(header.tostring()) + (np.prod(shape) * np.abs(header['BITPIX']//8)) - 1
    print(p)
    fobj.seek(p)
    fobj.write(b'\0')

  return


dims = (512, 1, 15, 1550, 1498)
fname = "test_torem.fits"
create_nd_fits_file(fname, dims)


with fits.open(fname, mode="update") as hdul:
    hdul.info()
    nn,_,_,_,_ = hdul[0].data.shape
    for i in range(nn):
      print(i)
      hdul[0].data[i,:,:,:] = hdul[0].data[i,:,:,:] + (i+1)
    
