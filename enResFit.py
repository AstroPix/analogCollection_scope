import scipy
from scipy.optimize import curve_fit
import pylab 
import numpy as np
import matplotlib.pyplot as plt
import h5py


# Define the fit function (a Gaussian)
def Gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))

def GaussFit(data, scale, savePlots, lowFit=0, highFit=np.inf):
    _, xbins0 = np.histogram(data)
    xmax=xbins0.max()
    #bins set by scope resolution
    bins=int(xmax/scale)
    histRange=np.linspace(0,xmax,bins)
    ydata,xbins1=np.histogram(data,bins=histRange)
    xbins=xbins1[:-1]
    xspace = np.linspace(0, xmax, 10000) # This creates a smoother plot when plotting the fit
    #AMANDA - Instead of defining number of bins, keep bins defined by scope resolution

    # Guesses for p01: [Amplitude, Mu, Sigma]
    muGuess = xbins[ydata.argmax()]
    if (highFit<np.inf): 
        muGuess = (highFit-lowFit)/2 + lowFit #midpoint of range
    elif (lowFit>0):
        muGuess = (xmax-lowFit)/2 + lowFit #can't algebra with infinity
    p01 = [ydata.max(), muGuess, muGuess/2.]

    #AMANDA - include check that fit converges
    #AMANDA - print chi2/pval from fit
    # popt returns the best fit values for amplitude, mean, and sigma. pcov returns a covariance matrix; the diagonal of this matrix returns the errors associated with the three returned values, which is used to determine the error in the energy resolution.
    popt, pcov = curve_fit(Gauss, xdata=xbins, ydata=ydata, p0=p01, maxfev=5000, bounds=([0, lowFit, 0], [np.inf, highFit, np.inf]))

    plt.figure(figsize = (8, 6))

    # Plot the data
    plt.bar(xbins, ydata, width=xbins[1] - xbins[0], color='blue', label=r'Data')
    # Plot the fit
    plt.plot(xspace, Gauss(xspace, *popt), 'r-', label='Fit')
    (Amp, Mu, Sigma) = popt
    # Print the outputs 
    print("Amplitude = %0.4f, Mu = %0.4f, Sigma = %0.4f" %(Amp, Mu, Sigma))
    energy_res = (2.355*Sigma*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
    print("The energy resolution is approximately %0.2f percent." %(abs(energy_res)))
    plt.xlim([0,xmax])
    plt.legend()
    plt.xlabel('Energy (V)')
    plt.ylabel('Counts')


    if savePlots:
        energy_res_str = str(abs(round(energy_res,2)))
        plt.savefig(path+outName+'_'+dsName+'_'+energy_res_str+'enRes.png')
    
    plt.show()    
    
    return popt, energy_res, pcov
    
    