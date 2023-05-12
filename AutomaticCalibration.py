# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 13:50:34 2023

@author: ninac
"""
#Automating Polarization Image Calibration 
#Spring 2023 QHY550P Research
#Author: Nina Christenson
#Version 1 Goals: to automatically dark subtract and flat field 
#the science images produced by the QHY550P Camera (if dark images have same exposure time as science images)

#packages
from pathlib import Path
import numpy as np
from astropy.nddata import CCDData
from astropy import units as u
import ccdproc 
import os

datadir = "C:/Users/ninac/spring23research/footprint_nebula" 

# working so that this is a function that will do p000,p045,p090, and p135 separately

angle = "p000" # this will be an argument you can pass in
biasFolder = datadir + "/bias/" + angle +"_bias"
combinedFolder = Path('.', 'combined')
# biasFilenameList = os.listdir(biasFolder)

# numBiasFiles = len(biasFilenameList)

# listCCDDataObjects = [] #all fits files must be read in as CCDData types

# for i in range(numBiasFiles):
#     f = CCDData.read(biasFolder + "/"+biasFilenameList[i], unit = "adu")
#     listCCDDataObjects.append(f)


biasCollection = ccdproc.ImageFileCollection(location = biasFolder,find_fits_by_reading = True)
biasData = biasCollection.data()

#c = ccdproc.Combiner(biasCollection.data())

c = ccdproc.combine(biasData)



# Combiner objects requires CCDData objects, loop produces a list of CCDData objects, but I don't know how 
#to put those objects in a folder in my directory to then put them in an ImageFIleCollection to eventually
# use the combiner function


