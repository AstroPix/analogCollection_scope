import scipy
from scipy.optimize import curve_fit
import pylab 
import numpy as np
import matplotlib.pyplot as plt


# Define the fit function (a Gaussian)
def Gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))

def SetArray(xmax):
    #AMANDA - optimize x axis range and bin number (so not hardcoded to 1 and 1000)
    bins = np.linspace(0,1,500) # Here, the number of bins has been set to 100 and the length of the x axis has been set to one. This can be adjusted; inspect the matplotlib histogram to see what the binning / x axis should be.
    fit_hist, bins_1 = np.histogram(peaks, bins=bins)
    bins_2 = np.array([bins[i] for i in range(len(bins)-1)])
    #bins_2 = 0.5*(bins[0:-1]+bins[1:])

def GaussFit(data, savePlots, bins=500, xmaxin=1, lowFit=0, highFit=1):
    _, xbins0 = np.histogram(data)
    xmax=xbins0.max()
    histRange=np.linspace(0,xmaxin,bins)
    ydata,xbins1=np.histogram(data,bins=histRange)
    xbins=xbins1[:-1]
    #AMANDA - WIP user defined fit window
    #fitRange=np.linspace(lowFit,highFit,1000)
    xspace = np.linspace(0, xmax, 10000) # This creates a smoother plot when plotting the fit

    # Guesses for p01: [Amplitude, Mu, Sigma]. These guesses must be reasonable.
    mean=np.mean(ydata)
    p01 = [ydata.max(), mean, mean/2.]

    #AMANDA - include check that fit converges
    # popt returns the best fit values for amplitude, mean, and sigma. pcov returns a covariance matrix; the diagonal of this matrix returns the errors associated with the three returned values, which is used to determine the error in the energy resolution.
    popt, pcov = curve_fit(Gauss, xdata=xbins, ydata=ydata, p0=p01, maxfev=5000)

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
    
    