import matplotlib.pyplot as plt
import numpy as np
import h5py
import os, glob
import scipy
from scipy.optimize import curve_fit
from scipy import special, interpolate
from scipy.integrate import quad



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

def getSaveto():
	return "/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp1_peaks/linear/"







#Plot data, fit to Gaussian, calculate energy resolution and draw plot
def enResPlot(settings, integral=0, fitLow=0, fitHigh=np.inf, dataset='run1'):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	energy=settings[3]
	savePlots=settings[4]

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
	if xBinWidth<1e-5: #shouldn't be - failsafe
		xBinWidth=0.002
	xMax=np.max(data)
	if integral>0:
		xMin=np.min(data)
	else:
		xMin=0
	if integral>0:
		xBinWidth=integral
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
	sigGuess=muGuess/2.
	p01 = [ampGuess, muGuess, sigGuess]
	
	#Fit histogram over desired range
	try:
		popt, pcov = curve_fit(Gauss, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000)
	except RuntimeError: #fit could not converge
		sigGuess=sum(ydata*(binCenters-muGuess)**2)
		p01 = [ampGuess, muGuess, sigGuess]
		popt, pcov = curve_fit(Gauss, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000)
	(Amp,Mu,Sigma)=popt
	#range is set with low_i and high_i index values for input arrays
	#bounds keeps all parameters positive
	
	#Calculate energy resolution
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	
	#Calculate N (events under fit)
	integ=scipy.integrate.quad(Gauss, -np.inf, np.inf, args=(Amp,Mu,Sigma))
	#returns integral and uncertainty on calculation - only return integral value
	
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
	#saveto=file[:-5]
	saveto=f"{getSaveto()}{title}{datain}_{energy}line.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()
	
	f.close()

	return popt, enRes, pcov, integ[0]
	
	
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
def printParams(settings, integ, popt, en_res, pcov, integral=False):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	energy=settings[3]
	savePlots=settings[4]
	
	#Distinguish between peaks and integral
	if integral:
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
		saveto=f"{getSaveto()}{title}_{energy}line{datain}EnRes.txt"
		k=open(saveto, "w")
		k.write("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		k.write("Events under fit (N) = %0.3f" %(integ)+"\n")
		k.write("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		k.write("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		k.write("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")
		k.close()
		#Display contents to terminal
		print(f"FIT FROM {title}_{energy}line_{datain}")
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		print("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		print("Events under fit (N) = %0.3f" %(integ)+"\n")
		print("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		print("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		print("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")


#Use linear fit coefficients to scale data to keV before plot/fit
def enResPlot_scale(settings, coef, fitLow=0, fitHigh=np.inf, dataset='run1', integral=0):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	energy=settings[3]
	savePlots=settings[4]

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
	#remove noise - neglect peak heights with <20 mV
	if integral==0:
		data=[y for y in data if y > 0.02]
	#linear fit
	data=[(coef[0]*x+coef[1]) for x in data]
	#quadratic fit
	#data=[(coef[0]*x*x+coef[1]*x+coef[2]) for x in data]
	#3rd deg poly fit
	#data=[(coef[0]*x*x*x+coef[1]*x*x+coef[2]*x+coef[3]) for x in data]
	#sqrt fit
	#data=[(coef[0]*np.sqrt(x)+coef[1]) for x in data]
	#spline
	#data=interpolate.splev(data, coef)
	
	#Create arrays for binning based on scope resolution
	xBinWidth=0.5 #1.0keV bins
	xMax=np.max(data)
	#give negative space to see full distribution
	if integral>0:
		xMin=np.min(data)
	else:
		xMin=0
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
	plt.title(f"Calibrated {title}, pixel {pixel}")
		
	
	#save figure
	saveto=f"{getSaveto()}{title}{datain}_{energy}line_calibrated.pdf"
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
	energy=settings[3]
	savePlots=settings[4]

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
	
	#Calculate N (events under fit)
	integ=scipy.integrate.quad(Edge, -2*Sigma, 2*Sigma, args=(Amp1,Amp2,Mu,Sigma))
	
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
	saveto=f"{getSaveto()}{title}{datain}EdgeFit_{energy}edge.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()
	
	f.close()

	return popt, enRes, pcov, integ[0]
	
#Calculate error on fit parameters, save in text file	
def printParams_edge(settings, integ, popt, en_res, pcov, integral=False):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	energy=settings[3]
	savePlots=settings[4]

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
		saveto=f"{getSaveto()}{title}_{energy}edge{datain}EnRes.txt"
		k=open(saveto, "w")
		k.write("Amplitude erfc= %d Amplitude const = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp1, Amp2, Mu, Sigma)+"\n")
		k.write("Events under fit (N) = %0.3f" %(integ)+"\n")
		k.write("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		k.write("Error in erfc amplitude = %0.3f \nError in erfc constant = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp1_err, Amp2_err, Mu_err, Sigma_err)+"\n")
		k.write("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")
		k.close()
		#Display contents to terminal
		print(f"FIT FROM {title}_{energy}edge_{datain}")
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		print("Amplitude erfc= %d Amplitude const = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp1, Amp2, Mu, Sigma)+"\n")
		print("Events under fit (N) = %0.3f" %(integ)+"\n")
		print("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		print("Error in erfc amplitude = %0.3f \nError in erfc constant = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp1_err, Amp2_err, Mu_err, Sigma_err)+"\n")
		print("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")
		
		
		
def getVals_fromTxt(inDir):		
	energyList, muArr1, sigmaArr1, nArr1, enResArr1 = [],[],[],[],[]

	os.chdir(inDir)
	peakFiles = glob.glob('*peaks*.txt')
	for filename in peakFiles:
		energyList.append(float(filename.split('_')[1][:-4]))
		openFile=open(filename,'r')
		lines=openFile.readlines()
		muArr1.append([float(lines[1].split(' = ')[-1])])
		sigmaArr1.append([float(lines[2].split(' = ')[-1])])
		nArr1.append([float(lines[3].split(' = ')[-1])])
		enResArr1.append([float(lines[4].split(' = ')[-1][:-2])])#eliminate % sign at the end
		
	intFiles = glob.glob('*integral*.txt')
	for filename in intFiles:
		energy=float(filename.split('_')[1][:-4])
		energyIndex=energyList.index(energy)
		openFile=open(filename,'r')
		lines=openFile.readlines()
		muArr1[energyIndex].append(float(lines[1].split(' = ')[-1]))
		sigmaArr1[energyIndex].append(float(lines[2].split(' = ')[-1]))
		nArr1[energyIndex].append(float(lines[3].split(' = ')[-1]))
		enResArr1[energyIndex].append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end

	return energyList, muArr1, sigmaArr1, nArr1, enResArr1
	
def getCalibVals_fromTxt(inDir, ele):		
	energyList, muArr1, sigmaArr1, nArr1, enResArr1 = [],[],[],[],[]
	fits=['linear','quad','tri','sqrt','spline1','spline3']

	for fit in fits:
		os.chdir(inDir+fit+'/')
		peakFile = glob.glob('*'+ele+'*peaks*.txt') #returns array with length 1
		energyList.append(float(peakFile[0].split('_')[1][:-4]))
		openFile = open(peakFile[0],'r')
		lines=openFile.readlines()
		muArr1.append(float(lines[1].split(' = ')[-1]))
		sigmaArr1.append(float(lines[2].split(' = ')[-1]))
		nArr1.append(float(lines[3].split(' = ')[-1]))
		enResArr1.append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end
	
		"""
		#when integral-calibrated plots are made
		muArr1.append([float(lines[1].split(' = ')[-1])])
		sigmaArr1.append([float(lines[2].split(' = ')[-1])])
		nArr1.append([float(lines[3].split(' = ')[-1])])
		enResArr1.append([float(lines[4].split(' = ')[-1][:-2])])#eliminate % sign at the end
	

		intFile = glob.glob('*'+ele+'*integral*.txt') #returns array with length 1
		print(intFile)
		energyList.append(float(peakFile[0].split('_')[1][:-4]))
		openFile = open(peakFile[0],'r')
		print(openFile)
		lines=openFile.readlines()
		muArr1[energyIndex].append(float(lines[1].split(' = ')[-1]))
		sigmaArr1[energyIndex].append(float(lines[2].split(' = ')[-1]))
		nArr1[energyIndex].append(float(lines[3].split(' = ')[-1]))
		enResArr1[energyIndex].append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end
		"""
	
	return energyList[0], muArr1, sigmaArr1, nArr1, enResArr1, fits
		
	
def calcError(sigmaArr, nArr):	
	errArr=[]
	#error Arr1ay = sigma/sqrt(N) (for edge, 2sig integral from mu)
	for j in range(len(sigmaArr)):
		err_p = -1 if nArr[j][0]<=0 else sigmaArr[j][0]/np.sqrt(nArr[j][0])
		err_i = -1 if nArr[j][1]<=0 else sigmaArr[j][1]/np.sqrt(nArr[j][1])
		errArr.append([err_p,err_i])
	return errArr