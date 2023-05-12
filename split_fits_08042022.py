#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pre-reqs:
Requires functions defined in the file pol_functions.py
Save pol_functions.py in the same folder where you save this file.

Purpose:
Splits a QHY550PM FITS image into 4 subimages after cropping out the 
boundaries and selecting the four polarization orientations.

Usage:
    python split_fits.py path-to-input-fits-file

Output:
    Four fits subimages with _pXXX suffixes where XXX are 
    000, 045, 090, 135
"""

# updated by Nina 8/4/2022
'''
loops through file names to split all files in a given folder 
'''


import os
# import all functions in file pol_functions.py
import pol_functions as pf

dirName = "Example/Directory"

fileList = os.listdir(dirName)  #does not check if its a fits file

for i in range(len(fileList)):
    ip = dirName + fileList[i]
    print(ip)
    ret = pf.read_fits_save_scaled_split_fits (ip)
    
    

    