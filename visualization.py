#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import os
import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D

#Input: four FITS files 
#Output: PDF of various visualization plots and three FITS files of the AoLP, DoLP, and PA calculations

#directory pointing to input files
datadir = "Example/Directory" 


#input filenames in fits format - one for each angle

ipfile_p000 = 'p000.fits'   
ipfile_p045 = 'p045.fits'
ipfile_p090 = 'p090.fits'
ipfile_p135 = 'p135.fits'

# joining directory and input file
file_p000 = os.path.join(datadir, ipfile_p000)  
file_p045 = os.path.join(datadir, ipfile_p045)
file_p090 = os.path.join(datadir, ipfile_p090)
file_p135 = os.path.join(datadir, ipfile_p135)


# output file name and directory
opfnam = 'EggNebula_Results_Github' + '.pdf'
opfile = os.path.join(datadir, opfnam)


# This function reads an input FITS image from the QHY550 camera 
# and returns a 2D array with the pixel values in the image. While saving 
# make sure the fits image is not stretched to 16 bits/pixel. Rather, the
# fits should retain 12 bits/pixel from the QHY550 camera. 

def read_fits(fitsfile):
    
    img_array_raw = fits.getdata(fitsfile, ext=0)
    img_data = img_array_raw
    
    return img_data



# getting FITS image data for all 4 orientations
o_000 = read_fits(file_p000)
o_045 = read_fits(file_p045)
o_090 = read_fits(file_p090)
o_135 = read_fits(file_p135)

# crop input files using array cropping to zoom in on object

#size of cropped image
sizex = 201
sizey = 201
size = (sizex, sizey)

#center position 
position = (675,516)

# makes cutout object
o_000_crop = Cutout2D(o_000, position, size)   
o_045_crop = Cutout2D(o_045, position, size)
o_090_crop = Cutout2D(o_090, position, size)
o_135_crop = Cutout2D(o_135, position, size)

# issue- Cutout2D not recognizing WCS in input file header when ip files are plate solved

#re-define data array in relation to cutout
o_000 = o_000_crop.data
o_045 = o_045_crop.data
o_090 = o_090_crop.data
o_135 = o_135_crop.data


#bits per pixel
bpp = 12

#full well pixel value
saturation = 2**bpp - 1

# angle to rotate image is specific to authors egg nebula imaging session
 
plate_scale = 0.000174 #degree/pix 
angle_rot = 168.5 #angle in degrees

# Get overall pixel stats for each orientation
print("Stats: 0 degrees")

avg_pixval0 = np.average(o_000)
minpixval0 = np.amin(o_000)
maxpixval0 = np.amax(o_000)
maxpixpct0 = 100.0*maxpixval0/saturation

print(' Average pixel value (ADU):', f'{avg_pixval0:.5}')
print(' Minimum pixel value (ADU):', minpixval0)
print(' Maximum pixel value (ADU):', maxpixval0,
      '(', f'{maxpixpct0:.3}','% of saturation)')


print("Stats: 45 degrees")

avg_pixval45 = np.average(o_045)
minpixval45 = np.amin(o_045)
maxpixval45 = np.amax(o_045)
maxpixpct45 = 100.0*maxpixval45/saturation

print(' Average pixel value (ADU):', f'{avg_pixval45:.5}')
print(' Minimum pixel value (ADU):', minpixval45)
print(' Maximum pixel value (ADU):', maxpixval45,
      '(', f'{maxpixpct45:.3}','% of saturation)')

print("Stats: 90 degrees")

avg_pixval90 = np.average(o_090)
minpixval90 = np.amin(o_090)
maxpixval90 = np.amax(o_090)
maxpixpct90 = 100.0*maxpixval90/saturation

print(' Average pixel value (ADU):', f'{avg_pixval90:.5}')
print(' Minimum pixel value (ADU):', minpixval90)
print(' Maximum pixel value (ADU):', maxpixval90,
      '(', f'{maxpixpct90:.3}','% of saturation)')

print("Stats: 135 degrees")

avg_pixval135 = np.average(o_135)
minpixval135 = np.amin(o_135)
maxpixval135 = np.amax(o_135)
maxpixpct135 = 100.0*maxpixval135/saturation

print(' Average pixel value (ADU):', f'{avg_pixval135:.5}')
print(' Minimum pixel value (ADU):', minpixval135)
print(' Maximum pixel value (ADU):', maxpixval135,
      '(', f'{maxpixpct135:.3}','% of saturation)')


# Ignore brightnesses beyond linearity
nonlinear = 0.95*saturation


######################################################################
# DoLP and AoLP computations
S0 = (o_000 + o_045 + o_090 + o_135)/2
S1 = o_000 - o_090
S2 = o_045 - o_135
DoLP = 100.0*(np.sqrt((S1*S1) + (S2*S2)))/S0
AoLP = np.rad2deg( 0.5 * np.arctan2(S2, S1) )

# PA computation

pos_angle = np.empty([sizex,sizey])   

for j in range(len(AoLP)):
    
    for i in range(len(AoLP)):
       
        
        if ((AoLP[i][j] >= -11.5).any() and (AoLP[i][j]<=90).any()): 
            
            y = 90 - (AoLP[i][j] -angle_rot) -180
            pos_angle[i][j] = y
            
        elif ((AoLP[i][j] >= -90).any() and (AoLP[i][j] <= -11.51).any()):
            y = 90-(AoLP[i][j] - angle_rot) -360
            pos_angle[i][j] = y

######################################################################
# Visualization
from copy import copy
import matplotlib.pyplot as plt
from astropy.visualization import (ManualInterval, MinMaxInterval, ZScaleInterval,
                                    PercentileInterval, LinearStretch, ImageNormalize)
from matplotlib.backends.backend_pdf import PdfPages
import cmocean
from matplotlib.colors import hsv_to_rgb

# Set up a colormap so that saturated/nonlinear pixels are red:
# use copy so that we do not mutate the global colormap instance
mycmap = copy(plt.cm.gray)
mycmap.set_bad('r', 1.0)

# Location and size of the main image [xmin, ymin, dx, dy]   
main_axis = [0.05, 0.05, 0.75, 0.85]


# Radial grid
rmin, rmax, rpts = 0.7, 1, 100
radii = np.linspace(rmin, rmax, rpts)

# theta values on the right side of the color circle
thpts = 500 
azimuthsR = np.linspace(-90, 91, thpts)
valuesR =  azimuthsR * np.ones((rpts, thpts))
###############




# Different normalizations
# myzsnorm = ImageNormalize(img_array, interval=ZScaleInterval(),
#                       stretch=LinearStretch())
# mymmnorm = ImageNormalize(img_array, interval=MinMaxInterval(), 
#                       stretch=LinearStretch())
# mypctnorm = ImageNormalize(img_array, interval=PercentileInterval(99), 
#                       stretch=LinearStretch())
# mynorm = ImageNormalize(img_array, interval=ManualInterval(0, nonlinear), 
#                       stretch=LinearStretch())

print('Creating plots...')


with PdfPages(opfile) as pdf:
   

    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    cbar_label = 'Counts (ADU)'
    ax0.set_title('Subimage with pixels oriented $0^\circ$, measured CCW from horizontal')
    #mynorm = ImageNormalize(o_000, interval=MinMaxInterval(),
                        # stretch=LinearStretch())
    myzsnorm = ImageNormalize(o_000, interval=ZScaleInterval(),
                          stretch=LinearStretch())
    im = ax0.imshow(o_000, origin='upper', norm=myzsnorm, cmap=mycmap)
   
    ax1 = fig.add_axes([0.85,0.10, 0.04,0.75])
    cbar = fig.colorbar(im, cax=ax1)
    cbar.set_label(cbar_label)
    pdf.savefig(fig)
    plt.close()


    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    cbar_label = 'Counts (ADU)'
    ax0.set_title('Subimage with pixels oriented $90^\circ$, measured CCW from horizontal')
    myzsnorm = ImageNormalize(o_090, interval=ZScaleInterval(),
                          stretch=LinearStretch())
    im = ax0.imshow(o_090, origin='upper', norm=myzsnorm, cmap=mycmap)
   
    ax1 = fig.add_axes([0.85,0.10, 0.04,0.75])
    cbar = fig.colorbar(im, cax=ax1)
    cbar.set_label(cbar_label)
    pdf.savefig(fig)
    plt.close()


    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    cbar_label = 'Counts (ADU)'
    ax0.set_title('Subimage with pixels oriented $45^\circ$, measured CCW from horizontal')
    myzsnorm = ImageNormalize(o_045, interval=ZScaleInterval(),
                          stretch=LinearStretch())
    im = ax0.imshow(o_045, origin='upper', norm=myzsnorm, cmap=mycmap)
   
    ax1 = fig.add_axes([0.85,0.10, 0.04,0.75])
    cbar = fig.colorbar(im, cax=ax1)
    cbar.set_label(cbar_label)
    pdf.savefig(fig)
    plt.close()


    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    cbar_label = 'Counts (ADU)'
    ax0.set_title('Subimage with pixels oriented $-45^\circ$, measured CCW from horizontal')
    myzsnorm = ImageNormalize(o_090, interval=ZScaleInterval(),
                      stretch=LinearStretch())
    im = ax0.imshow(o_135, origin='upper', norm=myzsnorm, cmap=mycmap)
   
    ax1 = fig.add_axes([0.85,0.10, 0.04,0.75])
    cbar = fig.colorbar(im, cax=ax1)
    cbar.set_label(cbar_label)
    pdf.savefig(fig)
    plt.close()


   
    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Histogram of pixel intensities')
    kwargs = dict(histtype='step', alpha=0.99, bins=100, log='True')
    ax.hist(o_000.ravel(), **kwargs, color='r', label='0 deg')
    ax.hist(o_045.ravel(), **kwargs, color='g', label='45 deg')
    ax.hist(o_090.ravel(), **kwargs, color='b', label='90 deg')
    ax.hist(o_135.ravel(), **kwargs, color='y', label='135 deg')
    ax.legend()
    #ax.set_xlim([0, 42])
    ax.set_xlabel('Pixel value (ADU)')
    ax.set_ylabel('Number of pixels')
    pdf.savefig(fig)
    plt.close()
   
   
    mycmap = copy(plt.cm.viridis)
    mycmap.set_bad('w', 1.0)   
    fig = plt.figure(figsize=(10, 8), dpi=200)
   
    dolpmin, dolpmax = np.amin(DoLP), np.amax(DoLP)
    #mynorm = ImageNormalize(DoLP, interval=ManualInterval(dolpmin, dolpmax),
    #                      stretch=LinearStretch())
    myzsnorm = ImageNormalize(DoLP, interval=ZScaleInterval(),
                          stretch=LinearStretch())
    ax0 = fig.add_axes(main_axis)
    ax0.set_title('Degree of Linear Polarization (white pixels were nonlinear and excluded)')
    im = ax0.imshow(DoLP, origin='upper', norm=myzsnorm, cmap=mycmap)

    ax1 = fig.add_axes([0.85,0.10, 0.04,0.75])
    cbar = fig.colorbar(im, cax=ax1)
    cbar.set_label('Degree of linear polarization (percent)')
    pdf.savefig(fig)
    plt.close()
   
   
    mycmap = copy(plt.cm.hsv)
    mycmap.set_bad('w', 1.0)   
    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    ax0.set_title('Angle of Linear Polarization (white pixels were nonlinear and excluded)')
    im = ax0.imshow(AoLP, origin='upper', cmap=mycmap)

    ax1 = fig.add_axes([0.72,0.45, 0.25,0.25], projection='polar')
    ax1.grid(False)
    ax1.axis('off')
    ax1.pcolormesh(azimuthsR*np.pi/180.0, radii, valuesR, cmap=mycmap)

    # Label AoLP angles
    for ii in np.arange(-90, 91, 30):
        iirad = ii*np.pi/180
        ax1.plot( (iirad, iirad), (rmax-0.03, rmax+0.00), color='k', ls='-')
        ax1.plot( (iirad, iirad), (rmin-0.00, rmin+0.03), color='k', ls='-')
        labl = str(ii) + "$^\circ$"
        if np.absolute(ii)==90:
            labl = "$\pm 90^\circ$"
        ax1.text(iirad, 1.20, labl, style='italic', fontsize=12, rotation=0, 
                horizontalalignment='center', verticalalignment='center')
       
    pdf.savefig(fig)
    plt.close()
   
    ##################################### PA #####################################
   
    mycmap = copy(plt.cm.hsv)
    mycmap.set_bad('w', 1.0)   
    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    ax0.set_title('Position Angle of Linear Polarization (white pixels were nonlinear and excluded)')
    im = ax0.imshow(pos_angle, origin='upper', cmap=mycmap)

    ax1 = fig.add_axes([0.72,0.45, 0.25,0.25], projection='polar')
    ax1.grid(False)
    ax1.axis('off')
    ax1.pcolormesh(azimuthsR*np.pi/180.0, radii, valuesR, cmap=mycmap)

    # Label PA angle
    for ii in np.arange(-90, 91, 30):
        iirad = ii*np.pi/180
        ax1.plot( (iirad, iirad), (rmax-0.03, rmax+0.00), color='k', ls='-')
        ax1.plot( (iirad, iirad), (rmin-0.00, rmin+0.03), color='k', ls='-')
        labl = str(ii) + "$^\circ$"
        if ii==90:
            labl = "$90^\circ(East) $"
           
        if ii==-90:
            labl = "$-90^\circ(West) $"
        ax1.text(iirad, 1.20, labl, style='italic', fontsize=12, rotation=0, 
                horizontalalignment='center', verticalalignment='center')
       
    pdf.savefig(fig)
    plt.close()
   


    H = 2*pos_angle/360 + 0.5   # AOLP/PA values are between [-90, 90] and hue needs to be [0,1]   
    V = (DoLP - dolpmin) / (dolpmax - dolpmin)
    S = np.ones_like(V)
   
    HSV = np.dstack((H,S,V))
    RGB = hsv_to_rgb(HSV)
    to_mask = np.maximum.reduce([o_000, o_045, o_090, o_135])
    to_mask3 = np.dstack((to_mask, to_mask, to_mask))
    RGB3 = np.where(to_mask3 >= nonlinear, 1, RGB)
    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    ax0.set_title('PA and DoLP combined')
    im = ax0.imshow(RGB3, origin='upper')
   
    hh, vv = np.mgrid[0:1:300j, 0:1:100j]
    ss = np.ones_like(vv)
    myHSV = np.dstack((hh,ss,vv))
    myRGB = hsv_to_rgb(myHSV)
    ax1 = fig.add_axes([0.87,0.10, 0.12,0.75])
    im = ax1.imshow(myRGB, origin="lower", extent=[dolpmin, dolpmax, -90, 90], aspect='auto')
    ax1.set_xlabel('DoLP (%)')
    ax1.set_ylabel('Position Angle (degrees)')
    pdf.savefig(fig)
    plt.close()
   
   
   
    mycmap = copy(cmocean.cm.phase)
    mycmap.set_bad('w', 1.0)
    ocean_colors = mycmap(H)
    oRGB = ocean_colors[:,:,:3]
    dolp3 = np.dstack((V, V, V))
    oRGB3 = np.multiply(oRGB, dolp3)
    oRGB3_masked_white = np.where(to_mask3 >= nonlinear, 1, oRGB3)
    fig = plt.figure(figsize=(10, 8), dpi=200)
    ax0 = fig.add_axes(main_axis)
    ax0.set_title('PA and DoLP combined - Different colormap')
    im = ax0.imshow(oRGB3_masked_white, origin='upper')
       
    ohh = cmocean.cm.phase(hh)[:,:,:3]
    vv3 = np.dstack((vv,vv,vv))
    ohh = np.multiply(vv3, ohh)
    ax1 = fig.add_axes([0.87,0.10, 0.12,0.75])
    im = ax1.imshow(ohh, origin="lower", extent=[dolpmin, dolpmax, -90, 90], aspect='auto')
    ax1.set_xlabel('DoLP (%)')
    ax1.set_ylabel('PA (degrees)')
    pdf.savefig(fig)
    plt.close()
   
   
   
    d = pdf.infodict()
    d['Title'] = 'Polarization results'
    d['Author'] = 'Wheaton Physics and Astronomy'   

# save AoLP, DoLP, PA as new fits files

print("Creating AoLP, DoLP, and PA FITS files...")
hdu_aolp = fits.PrimaryHDU(AoLP)
hdul_aolp = fits.HDUList([hdu_aolp])
hdul_aolp.writeto('AoLP_EggNebula_ccdproc.fits', overwrite= True)

hdu_dolp = fits.PrimaryHDU(DoLP)
hdul_dolp = fits.HDUList([hdu_dolp])
hdul_dolp.writeto('DoLP_EggNebula_ccdproc.fits', overwrite = True)

hdu_PA = fits.PrimaryHDU(pos_angle)
hdul_PA = fits.HDUList([hdu_PA])
hdul_PA.writeto('PA_EggNebula_Github.fits', overwrite = True)


print('Finished successfully.')





