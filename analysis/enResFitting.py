import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from scipy import special

#################
#
# Functions to help with fitting and calculating energy resolution / energy calibration
#
###################

def Gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))
    
def Edge(x,A,B,mu,sigma):
	return B*special.erfc((x-mu)/(np.sqrt(2)*sigma))+A
	#https://arxiv.org/pdf/1610.09185.pdf
    

def getArrayIndex(array, low,high):	
	goodValues=[]
	for i in range(len(array)):
		if array[i]>low and array[i]<high:
			goodValues.append(i)
			
	low_i=goodValues[0]
	high_i=goodValues[-1]
	
	return low_i, high_i

#Plot data, fit to Gaussian, calculate energy resolution and draw plot
def enResPlot(settings, integral=0, fitLow=0, fitHigh=np.inf, dataset='run1'):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	savePlots=settings[3]

	#Distinguish between peaks and integral
	if (integral>0):
		datain='_integral'
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName=dataset
	data=f[dsName+datain]
	scaling=f[dsName+'_scaling']
	#yzero=scaling[2]
	#ymult=scaling[3]
	#yoff=scaling[4]
	#print("YOFF")
	#print(yoff)
	#print("YZERO")
	#print(yzero)
	#print("YMULT")
	#print(ymult)
	
	#Create arrays for binning based on scope resolution
	xBinWidth=scaling[3]#YMULT Value
	xMax=np.max(data)
	if integral>0:
		xMin=np.min(data)
	else:
		xMin=0
	if integral>0:
		xBinWidth*=integral
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data
	
	#Create histogram of data
	hist=plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	ydata=hist[0]
	binCenters=hist[1]+xBinWidth/2
	binCenters=binCenters[:-1]
	
	#Set fit range
	low_i,high_i=getArrayIndex(binCenters,fitLow,fitHigh)
	
	#Set up fit with guesses for p01: [Amplitude, Mu, Sigma]
	muGuess=np.mean(data)
	if muGuess>fitHigh or muGuess<fitLow:
		muGuess=binCenters[low_i]+(binCenters[high_i]-binCenters[low_i])/2.
	ampGuess=ydata[low_i:high_i].max()
	p01 = [ampGuess, muGuess, muGuess/2.]
	
	#Fit histogram over desired range
	popt, pcov = curve_fit(Gauss, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000)
	(Amp,Mu,Sigma)=popt
	#range is set with low_i and high_i index values for input arrays
	#bounds keeps all parameters positive
	
	#Calculate energy resolution
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	
	#Display fit on final histogram
	plt.rc('text', usetex=True) #use Latex
	xspace=np.linspace(xMin,binEdges[-1],len(binCenters)*10)
	plt.plot(xspace, Gauss(xspace, *popt), 'r-', label=r'Fit')    
	plt.plot([], [], ' ', label=f"Energy res = {enRes:.2f}%")
	plt.plot([], [], ' ', label=r"Fit parameters:")
	plt.plot([], [], ' ', label=f"\small Amp={Amp:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.legend(loc="best")
	if integral:
		plt.xlabel('Integrated energy [V]')
	else:
		plt.xlabel('Peak Energy [V]')
	plt.ylabel('Counts')
	plt.title(f"Energy Resolution - {title}, pixel {pixel}")
		
	
	#save figure
	saveto=file[:-5]
	saveto=f"{saveto}{datain}EnRes{enRes:.2f}.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()
	
	f.close()

	return popt, enRes, pcov
	
	
#Return lowest recorded data value
def getMin(file, integral=0, dataset='run1'):

	#Distinguish between peaks and integral
	if (integral>0):
		datain='_integral'
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName=dataset
	data=f[dsName+datain]
	
	return np.min(data)
	
	
#Calculate error on fit parameters, save in text file	
def printParams(file, popt, en_res, pcov, savePlots,integral=False):
	#Distinguish between peaks and integral
	if (integral):
		datain='_integral'
	else: #peaks
		datain='_peaks'
	(Amp, Mu, Sigma) = popt
	(Amp_err, Mu_err, Sigma_err) = np.sqrt(np.diag(pcov))
	# Error propagation
	partial_sigma = (2.355*100)/Mu
	partial_mu = (2.355*100*Sigma)/(Mu**2)
	stdev_er = np.sqrt(((partial_sigma**2)*(Sigma_err**2))+((partial_mu**2)*(Mu_err)**2))
	
	if savePlots:
		saveto=file[:-5]
		saveto=f"{saveto}_run1{datain}EnRes.txt"
		k=open(saveto, "a")
		k.write("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		k.write("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		k.write("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		k.write("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")
		k.close()
		#Display contents to terminal
		print(f"FIT FROM {datain}")
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		print("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		print("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		print("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		print("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")


#Use linear fit coefficients to scale data to keV before plot/fit
def enResPlot_linearScale(settings, coef, integral=0, fitLow=0, fitHigh=np.inf, dataset='run1'):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	savePlots=settings[3]

	#Distinguish between peaks and integral
	if (integral>0):
		datain='_integral'
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName=dataset
	data=f[dsName+datain]
	#linear fit
	#data=[((x-coef[1])/coef[0]) for x in data]
	#quadratic fit
	data=[(coef[0]*x*x+coef[1]*x+coef[2]) for x in data]
	
	#Create arrays for binning based on scope resolution
	xBinWidth=1.0 #1.0keV bins
	xMax=np.max(data)
	#give negative space to see full distribution
	if integral>0:
		xMin=np.min(data)-10
	else:
		xMin=-10
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data
	
	#Create histogram of data
	hist=plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	ydata=hist[0]
	binCenters=hist[1]+xBinWidth/2
	binCenters=binCenters[:-1]
	
	#Set fit range
	low_i,high_i=getArrayIndex(binCenters,fitLow,fitHigh)
	
	#Set up fit with guesses for p01: [Amplitude, Mu, Sigma]
	muGuess=np.mean(data)
	if muGuess>fitHigh or muGuess<fitLow:
		muGuess=binCenters[low_i]+(binCenters[high_i]-binCenters[low_i])/2.
	ampGuess=ydata[low_i:high_i].max()
	p01 = [ampGuess, muGuess, muGuess/2.]
	
	#Fit histogram over desired range
	popt, pcov = curve_fit(Gauss, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000)
	(Amp,Mu,Sigma)=popt
	#range is set with low_i and high_i index values for input arrays
	#bounds keeps all parameters positive
	
	#Calculate energy resolution
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	
	#Display fit on final histogram
	plt.rc('text', usetex=True) #use Latex
	xspace=np.linspace(xMin,binEdges[-1],len(binCenters)*10)
	plt.plot(xspace, Gauss(xspace, *popt), 'r-', label=r'Fit')    
	plt.plot([], [], ' ', label=f"Energy res = {enRes:.2f}%")
	plt.plot([], [], ' ', label=r"Fit parameters:")
	plt.plot([], [], ' ', label=f"\small Amp={Amp:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.legend(loc="best")
	plt.xlabel('Calibrated Energy [keV]')
	plt.ylabel('Counts')
	plt.title(f"Energy Resolution - {title}, pixel {pixel}")
		
	
	#save figure
	saveto=file[:-5]
	saveto=f"{saveto}{datain}EnRes{enRes:.2f}.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()
	
	f.close()

	return popt, enRes, pcov
	
	
#Plot data, fit to Gaussian, calculate energy resolution and draw plot
def enResPlot_edge(settings, integral=0, fitLow=0, fitHigh=np.inf, dataset='run1'):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	savePlots=settings[3]

	#Distinguish between peaks and integral
	if (integral>0):
		datain='_integral'
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName=dataset
	data=f[dsName+datain]
	scaling=f[dsName+'_scaling']
	
	#Create arrays for binning based on scope resolution
	xBinWidth=scaling[3]#YMULT Value
	xMax=np.max(data)
	if integral>0:
		xMin=np.min(data)
	else:
		xMin=0
	if integral>0:
		xBinWidth*=integral
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data
	
	#Create histogram of data
	hist=plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	ydata=hist[0]
	binCenters=hist[1]+xBinWidth/2
	binCenters=binCenters[:-1]
	
	#Set fit range
	low_i,high_i=getArrayIndex(binCenters,fitLow,fitHigh)
	
	#Set up fit with guesses for p01: [Amplitude, Mu, Sigma]
	muGuess=np.mean(data)
	if muGuess>fitHigh or muGuess<fitLow:
		muGuess=binCenters[low_i]+(binCenters[high_i]-binCenters[low_i])/2.
	ampGuess=ydata[low_i:high_i].max()
	p01 = [ampGuess, ampGuess, muGuess, muGuess/2.]
	
	#Fit histogram over desired range
	popt, pcov = curve_fit(Edge, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000)
	(Amp1, Amp2, Mu, Sigma)=popt
	#range is set with low_i and high_i index values for input arrays
	#bounds keeps all parameters positive
	
	#Calculate energy resolution
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	
	#Display fit on final histogram
	plt.rc('text', usetex=True) #use Latex
	xspace=np.linspace(xMin,binEdges[-1],len(binCenters)*10)
	plt.plot(xspace, Edge(xspace, *popt), 'r-', label=r'Fit')  
	plt.plot([], [], ' ', label=f"Energy res = {enRes:.2f}%")  
	plt.plot([], [], ' ', label=r"Fit parameters:")
	plt.plot([], [], ' ', label=f"\small Amp1={Amp1:.2f}, Amp2={Amp2:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.legend(loc="best")
	if integral:
		plt.xlabel('Integrated energy [V]')
	else:
		plt.xlabel('Peak Energy [V]')
	plt.ylabel('Counts')
	plt.title(f"Edge Fit - {title}, pixel {pixel}")
		
	
	#save figure
	saveto=file[:-5]
	saveto=f"{saveto}{datain}EdgeFit_EnRes{enRes:.2f}.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()
	
	f.close()

	return popt, enRes, pcov
	
#Calculate error on fit parameters, save in text file	
def printParams_edge(file, popt, en_res, pcov, savePlots,integral=False):
	#Distinguish between peaks and integral
	if (integral):
		datain='_integral'
	else: #peaks
		datain='_peaks'
	(Amp1, Amp2, Mu, Sigma) = popt
	(Amp1_err, Amp2_err, Mu_err, Sigma_err) = np.sqrt(np.diag(pcov))
	# Error propagation
	partial_sigma = (2.355*100)/Mu
	partial_mu = (2.355*100*Sigma)/(Mu**2)
	stdev_er = np.sqrt(((partial_sigma**2)*(Sigma_err**2))+((partial_mu**2)*(Mu_err)**2))
	
	if savePlots:
		saveto=file[:-5]
		saveto=f"{saveto}_run1{datain}EnRes.txt"
		k=open(saveto, "a")
		k.write("Amplitude erfc= %d Amplitude const = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp1, Amp2, Mu, Sigma)+"\n")
		k.write("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		#amanda here
		k.write("Error in erfc amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		k.write("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")
		k.close()
		#Display contents to terminal
		print(f"FIT FROM {datain}")
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		print("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		print("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		print("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		print("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")