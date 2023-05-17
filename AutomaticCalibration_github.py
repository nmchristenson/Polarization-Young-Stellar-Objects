# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:50:34 2023

@author: ninac
"""
#Automating Polarization Image Calibration 
#Spring 2023 QHY550P Research
#Author: Nina Christenson

#Input: the angle for which you want calibrated science images (ex: p000) AND exposure time for dark images in line 86
# "combinedDark.header['EXPTIME']"

#*** main directory and all folders are specific to author - change to your filename convention***

#Output: Combined Bias, Bias-Subtracted Darks AND Combined Dark, Bias+Dark Subtracted Flat Fields 
#and Combined Flat Field (not normalized), calibrated science images (flat normalization occurs in this step)

import os
from astropy.nddata import CCDData
from astropy import units as u
from ccdproc import ImageFileCollection
from ccdproc import combine
import ccdproc 


#main directory
mainDir = "Example/Parent/Directory" 
angle = input("Angle: ")


#reading in bias folder
biasFolder = mainDir + "bias/" + angle +"_bias/"
#location to write combined bias file to
biasWriteTo = mainDir + "bias/"

#creating image file collection of bias files
bias_images = ImageFileCollection(biasFolder)

#using a median method to be consistent with AiJ calibration
combinedBias = combine(bias_images.files_filtered(include_path = True), method = 'median', unit = "adu")

#add header to indicate its a combined file
combinedBias.header['Combined'] = True

combinedBias.write(biasWriteTo + "CombinedBias_" + angle + ".fits", overwrite = True)

print("Combined Bias Created...")


#reading in dark images and making a combined dark image

#filepath to folder with dark images
darkFolder = mainDir +  "darks/" + angle + "_dark" 
darkWriteTo = mainDir + "darks/" + angle + "_biasreduced"

#this creates the folder to write the dark images if it does not exist already
if not os.path.exists(darkWriteTo):
    os.makedirs(darkWriteTo)
    
#create image collection of dark files
darkImCollection = ImageFileCollection(darkFolder)

#define header to indicate they are bias reduced 
darkMeta = dict()
darkMeta = {'Bias_Reduced':"True"}

#subtract bias from dark frames

for i, filename in darkImCollection.data(return_fname=True):
    dark = CCDData.read(darkFolder+ "/" + filename, unit = "adu")
    darkSubtract = ccdproc.subtract_bias(dark, combinedBias)
    darkSubtract.header = darkMeta
    darkSubtract.write(darkWriteTo + "/biasReduced_"+ filename, overwrite = True)
print("Dark Images Bias Subtracted...")

#make master dark
reducedDarkImCollection = ImageFileCollection(darkWriteTo)
combinedDark = combine(reducedDarkImCollection.files_filtered(include_path = True), method = 'median', unit = 'adu')    
combinedDark.header['Combined']= True

# USER SUPPLIED DARK EXPOSURE TIME !! Was previously having issues with getting ccdproc.subtract_dark to recognize 
# the exposure time keyword so it's currently user supplied 

combinedDark.header['EXPTIME'] = 60
combinedDark.write(darkWriteTo + "_CombinedDark_"+ angle+".fits", overwrite = True)

print("Combined Dark Created...")

# bias/dark subtract from flats

flatFolder = mainDir + "flats/"+angle+"_flats"
flatWriteTo = mainDir + "flats/"+angle+"_flats/bias_subtracted/"

if not os.path.exists(flatWriteTo):
    os.makedirs(flatWriteTo)

flatCollection = ImageFileCollection(flatFolder)

for j,filename in flatCollection.data(return_fname=True):
    flat = CCDData.read(flatFolder + "/" +filename, unit = "adu")
    flatBiasSubtract = ccdproc.subtract_bias(flat, combinedBias)
    flatDarkSubtract = ccdproc.subtract_dark(flatBiasSubtract, combinedDark, 
                                            exposure_time = 'EXPTIME', exposure_unit = u.second, scale = True)
    flatDarkSubtract.write(flatWriteTo + "calibrated_"+filename, overwrite = True)
    
print("Flat Images Calibrated...")

flatReduced = ImageFileCollection(flatWriteTo)
combinedFlat = combine(flatReduced.files_filtered(include_path = True), method = 'median', unit = 'adu')
combinedFlat.write(flatWriteTo + "CombinedFlat_"+angle+".fits", overwrite = True)


print("Combined Flat Created...")


# calibrate individual science images

scienceFolder = mainDir + "science_images/" + angle
scienceWriteTo = mainDir + "science_images/"+angle+"_calibrated"
if not os.path.exists(scienceWriteTo):
    os.makedirs(scienceWriteTo)

scienceImageCollection = ImageFileCollection(scienceFolder)

for k,filename in scienceImageCollection.data(return_fname = True):
    scienceImage = CCDData.read(scienceFolder + "/" + filename, unit = "adu")
    bSub = ccdproc.subtract_bias(scienceImage, combinedBias)
    dSub = ccdproc.subtract_dark(bSub, combinedDark, exposure_time = "EXPTIME", exposure_unit = "seconds")
    final_calibrated = ccdproc.flat_correct(dSub, combinedFlat)
    final_calibrated.write(scienceWriteTo + "/calibrated_"+filename, overwrite = True)

print("Science Images Calibrated...")



# #Appendix: Extra Code 


# #summary attribute is an astropy table that displays FITS header data as columns
# #this shows all FITS header keywords
# #print(darkImCollection.summary.colnames)

# #basic syntax for iterating over the methods of the image collection
# #for i,fname in darkImCollection.hdus(return_fname=True):
#     #print(i,i.header['blklevel'])

#different method of calibrating science images- this doesn't normalize flat

# for k,filename in scienceImageCollection.data(return_fname = True):
#     scienceImage = CCDData.read(scienceFolder + "/" + filename, unit ="adu")
#     processed = ccdproc.ccd_process(scienceImage,master_bias=combinedBias,
#                                     dark_frame = combinedDark,
#                                     master_flat = combinedFlat, exposure_key = "EXPTIME",
#                                     exposure_unit = u.second)
#     processed.write(scienceWriteTo + "/processed_"+filename, overwrite = True)

# print("Science Images Calibrated...")