# Polarization-Young-Stellar-Objects
Wheaton College Massachusetts - Advisor: Professor Dipankar Maitra (Physics and Astronomy Department)

Observational astronomy research project looking at the polarization of young stellar objects using the QHY550PM Polarization camera with the SONY IMX250MZR imaging sensor. This sensor has an integrated four-directional repeating polarizer array (0°, 45°, 90°, -45°); once an image has been taken, one can extract four subarrays that are all pixels with the same polarization filter. By use of the Stokes parameters, we are able to calculate the degree and angle of linear polarization for a given image. 

I worked on this project for my senior capstone as a Physics major at Wheaton College Massachusetts, completing 8 weeks of work in the summer before senior year, and continuing to work on the project until graduation. All images were taken using a Meade LX600 12-inch telescope at the Wheaton College Observatory.

Major goals:
  - Explore and quantify all camera settings (gain scale)
  - Identify targets for imaging (based on magnitude and known polarization of the object, as well as visibility from our location)
  - Complete successful imaging sessions (including science and calibration images)
  - Analysis: calibrate images in AstroImageJ and run Python visualization of the degree and angle of linear polarization of the object
  - Use Python ccdproc package to automate calibration process (bias subtracting, dark subtracting, flat fielding)


Included in this repository are the visualization and automatic calibration Python scripts and all of our defined convenience functions. 

Acknowledgements: We would like to thank NASA/Rhode Island Space Grant and Wheaton College for funding this research.
