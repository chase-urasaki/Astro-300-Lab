#To find the flux for the individual spectra lines
#Import the necessary libraries and packages
#We have to integrate the area under each of the plots that we are after. 
import numpy as np
from matplotlib import pyplot
from scipy.optimize import curve_fit
from scipy.integrate import quad

data = np.genfromtxt('DEIMOS_05.txt', delimiter = ',')

#Here we'll normalize the plot
Xs = data[:,0]
ymin = np.min(data[:,1])
Ys = (data[:,1]-ymin)/(np.max(data[:,1])-ymin)
#pyplot.plot(Xs,Ys)

#Insert redshift
rs = 1.+0.28

#Rest frame wavelength
wl = 6563.*rs

#Now, we want to define a Gaussian function to fit for 
def Gaussian(x, mu, sigma, b, sc):
    return b+sc*(0.39894/sigma)*np.exp(-(x-mu)**2/(2*sigma**2))
pyplot.plot(Xs,Ys)
pyplot.xlim(wl-5., wl+5.)
#pyplot.show()

#We want to give our fitting function a little help, so we use mu = wl*(1+rs)
popt, pcov = curve_fit(Gaussian, Xs, Ys, p0=[wl,0.5,0.1, 5])

pyplot.plot(Xs, Gaussian(Xs, popt[0], popt[1], popt[2], popt[3]), linestyle = ':', color = 'r')
pyplot.show()

Flux = np.absolute(popt[3])
error = np.sqrt(np.diagonal(pcov))
dFlux = error[3]

print('{0:4.2f} +/- {1:4.2f}'.format(Flux, dFlux))
