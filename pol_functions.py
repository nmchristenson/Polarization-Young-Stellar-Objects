#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 16:58:30 2022

Functions used for the polarization cameras

"""
import os
import sys
import numpy as np
from astropy.io import fits

def determine_filetype_opfilename (input_file):
    split_tup = os.path.splitext(input_file) 
    file_extension = split_tup[1]
    output_file = input_file + '.pdf'
    return file_extension, output_file

# This function reads an input FITS image from the QHY550 camera 
# and returns a 2D array with the pixel values in the image. While saving 
# make sure the fits image is not stretched to 16 bits/pixel. Rather, the
# fits should retain 12 bits/pixel from the QHY550 camera. It also determines
# which software created it, and crops the bad edges appropriately
def read_fits(fitsfile):
    img_array_raw = fits.getdata(fitsfile, ext=0)
    bits_per_pixel = 12
    
    # Determine whether the image was created using NINA or SharpCap
    sw = fits.getval(fitsfile, 'SWCREATE', ext=0)
    if (sw.find('SharpCap') != -1):
        print('Reading FITS image taken with SharpCap')
        img_data = img_array_raw[0:-22, 0:-12]
    elif (sw.find('N.I.N.A.') != -1):
        print('Reading FITS image taken with NINA')
        img_data = img_array_raw[22:, 12:]
    else:
        print('Neither SharpCap nor NINA! Quitting ...')
        sys.exit(0)

    return bits_per_pixel, img_data

def read_fits_save_scaled_split_fits (fitsfile):
    img_array_raw, header = fits.getdata(fitsfile, ext=0, header=True)
    
    # Create the names of the 4 output files
    split_tup = os.path.splitext(fitsfile)
    op_000 = split_tup[0] + '_p000.fits'
    op_045 = split_tup[0] + '_p045.fits'
    op_090 = split_tup[0] + '_p090.fits'
    op_135 = split_tup[0] + '_p135.fits'
    
    
    
    
    # Determine whether the image was created using NINA or SharpCap
    sw = fits.getval(fitsfile, 'SWCREATE', ext=0)
    if (sw.find('SharpCap') != -1):
        print('Reading FITS image taken with SharpCap')
        
        # Crop image to keep good region
        img_data = img_array_raw[0:-22, 0:-12]
        
        # Get rid of some junk header keywords inserted by SharpCap
        header.pop('COLORTYP', None)
        header.pop('XBAYROFF', None)
        header.pop('YBAYROFF', None)
        header.pop('YBAYROFF', None)
        header.pop('BAYOFFX', None)
        header.pop('BAYOFFY', None)
        header.pop('BAYERPAT', None)
    

    elif (sw.find('N.I.N.A.') != -1):
        print('Reading FITS image taken with NINA')
        
        # Crop image to keep good region
        img_data = img_array_raw[22:, 12:]
        
    else:
        print('Neither SharpCap nor NINA! Quitting ...')
        sys.exit(0)
        
    # Divide by 16 (= 2^16 / 2^12) because the FITS files are stretched!
    img_data_scaled = img_data/16
    img_data_scaled = img_data_scaled.astype('int16')
    
    # Extract the subimages
    orient_000 = img_data_scaled[1::2, 1::2]
    orient_045 = img_data_scaled[::2,  1::2]
    orient_090 = img_data_scaled[::2,   ::2]
    orient_135 = img_data_scaled[1::2,  ::2]
    
    # Write the subimages as FITS
    header['POL_ANG'] = 0
    fits.writeto(op_000, orient_000, header, overwrite=True)
    header['POL_ANG'] = 45
    fits.writeto(op_045, orient_045, header, overwrite=True)
    header['POL_ANG'] = 90
    fits.writeto(op_090, orient_090, header, overwrite=True)
    header['POL_ANG'] = 135
    fits.writeto(op_135, orient_135, header, overwrite=True)

    return 0


# This function reads raw images from PHX050 camera, decides the bit dept,
# and returns a 2D array with the pixel values in the image.
def read_raw(rawfile):
    img_height, img_width = 2048, 2448
    NPIX = img_height * img_width           # Total number of pixels

    # Read binary data from raw file into 8-bit unsigned integer numpy array
    binary_array = np.fromfile(rawfile, dtype='uint8')

    num_bytes = len(binary_array)
    if (num_bytes == NPIX):          # Mono8 format
    
        bpp = 8   # Bits/pixel
        print('Decoding 8-bit image:', rawfile)
        img_arr_1d = binary_array
        
    elif (num_bytes == 2*NPIX):      # Mono12 format

        bpp = 12   # Bits/pixel
        print('Decoding 12-bit image:', rawfile)
        
        # These are the even bytes, i.e. the 0th, 2nd, 4th, 6th, 8th, 10th, ...
        evn_arr = binary_array[0::2].copy()

        # These are the odd bytes, i.e. the 1st, 3rd, 5th, 7th, 9th, 11th, ...
        odd_arr = binary_array[1::2].copy()


        # Recast these arrays as 16-bit integer arrays because an 8-bit integer
        # array can store from 0 to 255 only whereas we need to store bigger 
        # numbers when working with 12-bit integers.
        evn_arr = evn_arr.astype('int16')
        odd_arr = odd_arr.astype('int16')

        # Shift the odd bytes 8 binary digits and then add to the even byte to 
        # get the 12-bit pixel value.
        pix_val = (odd_arr << 8) + evn_arr
        img_arr_1d = pix_val
        
    else:
        
        print('Unrecognized image format')
        sys.exit(0)
        
    # Reshape 1D array img_arr_1d to 2D image array data
    img_data = img_arr_1d.reshape(img_height, img_width)
    
    return bpp, img_data


# give a raw from PHX050, or fits from QHY550, returns bits per
# pixel, array of pixel values, and the name of output file.
def read_image_get_bpp (input_file):
    fileType ,opfile = determine_filetype_opfilename (input_file)
    
    if fileType == ".fits" :    # Images taken with QHY550 camera
        bpp, img_array = read_fits(input_file)
    elif fileType == ".raw" :   # Images taken with PHX050 camera
        bpp, img_array = read_raw(input_file)
        
    return bpp, img_array, opfile


    
# This function extracts the subimages corresponding to different
# polarizer orientations. 
def extract_subimages(full_img_array):
    orient_090 = full_img_array[::2,   ::2].astype('float64')
    orient_045 = full_img_array[::2,  1::2].astype('float64')
    orient_135 = full_img_array[1::2,  ::2].astype('float64')
    orient_000 = full_img_array[1::2, 1::2].astype('float64') 
    return orient_000, orient_045, orient_090, orient_135


# Given an array of pixel values and bit depth, return some pixel stats
def get_arrray_stats(image_array):
    avg_pix = np.average(image_array)
    med_pix = np.median(image_array)
    q1, q99 = np.percentile(image_array, [1, 99])
    minpix  = np.amin(image_array)
    maxpix  = np.amax(image_array)
    
    print(' mean and median pixel values (ADU):', f'{avg_pix:.3}', med_pix)
    print(' min, 1-pctile, 99-pctile,  max pixel values (ADU):', 
          minpix, q1, q99, maxpix)
    return avg_pix, med_pix, q1, q99, minpix, maxpix

