# coding: utf-8
#Photometry is the procedure of measuring the brightness of objects (usually stars or galaxies) in the nighttime sky. There are a number of possible programs and procedures to use, such as IDL, but we'll use python. In this case, we are comparing the flux of a star with that of a reference star whose magnitude we already know. Therefore, we will be taking the 'fluxes' of these two stars and using the ratio etween them to find out the magnitude of the target star. You will see them represented as the reference star (or ref) and the target star (tar).
#First we want to start by importing our necessary libraries
import numpy as np
from matplotlib import pyplot
from matplotlib import patches
from astropy.io import fits

#Now, we want to define a few variables 
#We want to define the hdulist as the opened fits file. The fits file is what image format pictures from telescopes come from. Just insert the file name of whatever the fits file is
hdulist = fits.open('filename')

#The header of the fits file contains important infomration like what telescope this file was from, what the date was, temp, airmass, exposure time etc. For now we are only interested in the date in which the image was taken. Usually the dates are given in as a Mean Julian Date of Observation
header = hdulist[0].header
date = header['MJD-OBS']

#Depending what kind of photometry you are doing, you might have a reference star. Here, insert the RA and Dec of the reference star in decimal degrees
decRAref  = 16.*15+14./4. + 20.912/240.
decDECref = -(19.+6./60.+5.5/3600.)

#Now, the RA and Dec of the target star, again in decimal degrees
decRA  = 16.*15+14./4. + 20.3/240.
decDEC = -(19.+6./60.+48.1/3600.)

#Here's where it gets tricky. So, obviously we are going to need to get the location of the star on the image itself. Now, the image itself is composed of pixels, each with it's own intensity read out. We are now going to define a functions that take the actual RA and DEC of the star and turn it into coordinates on the image. Thankfully, the header contains information we need to do such a thing
def RA2pix(RA, head):
    returnValue = head['CRPIX1']+(RA-head['CRVAL1'])/head['CD1_1']
    returnValue = int(returnValue)
    return returnValue

def DEC2pix(DEC, head):
    returnValue = head['CRPIX2']+(DEC-head['CRVAL2'])/head['CD2_2'] 
    returnValue = int(returnValue)
    return returnValue
#Now what we want to do is 'fix' the image. We define our function to do a few things. First, since our image from the telescope isn't perfect, what we want to do is make sure any overexposed parts of the image (pixels that have greater than a 65000 DN), are made irrelevant. The easiest way to do that is to assign any pixel with greater than 65,000 = 1. 
#added another variable for the background "XY"
def imagePlot(rawImage, X, Y, X1, Y1, XY):
    viewIm = rawImage-np.min(rawImage)+1
    viewIm[viewIm>65000] = 1

#After everything is nice and reasonable, we are going to go ahead and scale the image by log 10. Now, what this does is this actually makes the 'data' image dimmer than the actual one. But, ultimatey this won't matter, since we are looking at the ratios of brightness between stars. 
    viewIm = np.log10(viewIm)
#The folling command just tells pyplot we want to see the image in gray.
    pyplot.imshow(viewIm, 'gray')
    ax = pyplot.gca()

#What we are doing now is creating a 25 x 25 pixel aperture around each of the stars. More specifically, we are drawing the boxes around the area in which we are reading the DN of the pixels. 
    ax.add_patch(patches.Rectangle((X-12,Y-12),25,25,fill=False,color='red'))
    ax.add_patch(patches.Rectangle((X1-12,Y1-12),25,25,fill=False,color='cyan'))

#Of course, we create the extra box for the background of sky that we are using so that we can visually see that we don't have stars in the background
    ax.add_patch(patches.Rectangle((X-12+XY,Y-12+XY),25,25,fill=False,color='blue'))

#The following commands just creates a legend in the plot
    red = patches.Patch(color='red', label='Target Star')
    cyan = patches.Patch(color='cyan', label='Reference Star')
    blue = patches.Patch(color='blue', label='Background')
    pyplot.legend(handles=[red, cyan, blue])
    pyplot.show()
    pyplot.close()

#The region about bg pixels that are to the left and below our aperture
bg = 47

image = hdulist[0].data
pixRA = RA2pix(decRA, header)
pixDEC = DEC2pix(decDEC, header)
pixRAref = RA2pix(decRAref, header)
pixDECref = DEC2pix(decDECref, header)
imagePlot(image, pixRA, pixDEC, pixRAref, pixDECref, bg)
decRA, decDEC

#replaced the number of pixels moved in x and y with bg
bkgAp = image[pixDEC+bg-12:pixDEC+bg+12,pixRA+bg-12:pixRA+bg+12]
meanBkg = np.mean(bkgAp)
medianBkg = np.median(bkgAp)

#Compare both the mean and median to see which gives the better result. 
meanBkg, medianBkg

#Target star's aperture

rawAp = image[pixDEC-12:pixDEC+12, pixRA-12:pixRA+12]

# Subtract out the background

corAp = rawAp-medianBkg
totalCount = np.sum(corAp)

#To get the electron count from the DN, in the header there's a 'gain' which gives you the conversion between DN and electron number. 
gain = header['GAIN']
totalElectrons = totalCount*gain

#Using Poisson statistics we know that errors in measurements are equal to sqrt(n)
Error = np.sqrt(totalElectrons)

print('For target {0:8.2f}+/-{1:6.2f} electrons on MJD:{2:9.3f}'.format(totalElectrons, Error, date))

#now we wish to generate the electron count from the reference star
RefAp = image[pixDECref -12:pixDECref + 12, pixRAref - 12:pixRAref+12]
RefcorAp = RefAp-medianBkg
totalcountref = np.sum(RefcorAp)

totalelectronsRef = totalcountref*gain
ErrorRef = np.sqrt(totalelectronsRef)
print('For reference {0:8.2f} +/- {1:6.2f} electrons on MJD:{2:9.3f}'.format(totalelectronsRef,ErrorRef, date)) 

#Now we want to find the magnitude of the star.
m_ref = 13.5
m = -2.5*np.log10(totalElectrons/totalelectronsRef) + m_ref

#To find the error in this star, we use the following formula
#We first define the following varables to make it a bit easier
Uncertainty_tar = Error/totalElectrons
Uncertainty_ref = ErrorRef/totalelectronsRef

delm = (-2.5/np.log(10))*np.sqrt(Uncertainty_tar**2 + Uncertainty_ref**2)
print('Magnitude of the target star is {0:5.5f} +/- {1:5.5f}'.format(m, delm))

