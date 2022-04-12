import matplotlib.pyplot as plt
import numpy as np
import h5py
import os, glob, sys
import scipy
from scipy import special, interpolate
from scipy.optimize import curve_fit
from scipy.integrate import quad


#################
#
# Functions to help with fitting and calculating energy resolution / energy calibration
#
###################

#fit functions
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
	
#Return lowest recorded data value
def getMin(file, integral=False, dataset='run1'):
	scale=1.
	#Distinguish between peaks and integral
	if integral:
		datain='_integral'
		scaling=f[dsName+'_scaling']
		scale=scaling[1] #XINCR value in s (usually around 100 us)
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName=dataset
	data=f[dsName+datain]	
	minn=np.min(data)
	
	#if integral, scale by deltaT
	return minn*scale
	
def getSaveto(savedir=None):
	if savedir:
		return savedir
	else:
	#go to directory where this script (and runOptions) lives
		os.chdir(sys.path[0])
		#pull saveDir variable from runOptions
		try:
			with open("runOptions.txt", 'r') as runOptions:
				lines=runOptions.readlines()
				for line in lines:
					if "saveDir = " in line:
						saveDir = line.split(" = ")[1][:-1]
			runOptions.close()
			return saveDir
		except IOError:
			print("Input file does not exist")
			print("Please provide a save location: ")
			saveDir=input()
			return saveDir


def calc_chisquare(meas, sigma, fit):
	diff = pow(meas-fit, 2.)
	test_statistic = (diff / pow(sigma,2.)).sum()
	return test_statistic
	
def closest(lst, K):
	lst = np.asarray(lst)
	idx = (np.abs(lst - K)).argmin()
	return idx
	
def getGausRange(x,y,errs,popt):
	#y is array of bin content values
	rang=3*popt[2]
	bin_mean=closest(x,popt[1])
	bin_3sig=closest(x,popt[1]+rang)
	bin_3sigN=closest(x,popt[1]-rang)
	
	return x[bin_3sigN:bin_3sig],y[bin_3sigN:bin_3sig],errs[bin_3sigN:bin_3sig]
 

def iterativeFit(fitFn, p01, x, y, low, high, maxIt=25):
	#iterate on fit until max iterations or chi2/ndof<1 (fit within 1 sigma)
	errs=np.sqrt(y)
	#avoid division by zero
	for i,err in enumerate(errs):
		if err==0:
			errs[i]=1e-8
	goodness_last=1e100
	popt_best=[0,0,0]
	pcov_best=0

	i=0
	while i<maxIt:
		print("Fit iteration "+str(i))
		if (high-low<3):
			#too small a range to fit
			print("Too small a range to fit")
			break
		else:
			try:
				popt, pcov = curve_fit(fitFn, xdata=x[low:high], ydata=y[low:high], sigma=errs[low:high], p0=p01, bounds=(0,np.inf), maxfev=8000, absolute_sigma=True)
			except RuntimeError: #fit could not converge
				print("RuntimeError - retry with smaller sigma guess")
				p01[-1]=p01[-1]/2 #make sigma guess smaller
				popt, pcov = curve_fit(fitFn, xdata=x[low:high], ydata=y[low:high], sigma=errs[low:high], p0=p01, bounds=(0,np.inf), maxfev=8000, absolute_sigma=True)
			x1,y1,err1=getGausRange(x,y,errs,popt)
			TS = calc_chisquare(y1, err1, Gauss(x1,*popt))
			NDF = len(y1) - len(popt)
			if NDF==0:
				NDF=1e-8
			goodness = TS/float(NDF)
			#print("chisquare/NDF = {0:.2f} / {1:f} = {2:.2f}".format(TS, NDF, TS / float(NDF)))
		if goodness<goodness_last:
			popt_best=popt
			pcov_best=pcov
			goodness_last=goodness
			p01=popt
			#tighten fit range
			low+=2
			high-=2
			i+=1
		else:
			#fit getting worse
			break

	return popt_best, pcov_best


#Plot data, fit, calculate energy resolution and draw plot
#Defaults for fitting a photopeak from a measured spectrum considering peak heights
#Can fit pulse integral with optional integral input
#Can fit Compton Edge with integrated Gaussian with optional edge input
#Can calibrate measured signal to keV using calibration curve and fit calibrated spectrum with optional edge and fit inputs
def enResPlot(settings, integral=False, edge=False, injection=False, fitLow=0, fitHigh=np.inf, dataset='run1', fit=-1, coef=0, savedir=None, binSize=0):
	#Define inputs
	file=settings[0]
	title=settings[1]
	pixel=settings[2]
	energy=settings[3]
	savePlots=settings[4]
	chip=settings[5]

	#Distinguish between peaks and integral
	if integral:
		datain='_integral'
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName=dataset
	data=f[dsName+datain]
	
	scaling=f[dsName+'_scaling']
	#remove noise - neglect peak heights with <20(30) mV for amp1(amp2)
	if pixel==2 and chip==3:
		noiseCut=0.03
	elif pixel==1 and chip==4:
		noiseCut=0.01
	else:
		noiseCut=0.02
	if not integral:
		data=[y for y in data if y > noiseCut]
	else:
		#scale integral by scope resolution XINCR, noise cut values <0
		deltaT=scaling[1] #XINCR value in s (usually around 100 us)
		data=[y*deltaT*1e6 for y in data if y>0] #[V*ns]
	
	#if calibrating, plug values into calibration function (unless value falls outside of calibration - then discard)
	if fit>=0:
		data=[float(coef(x)) for x in data if float(coef(x))>0]
	
	#Create arrays for binning based on scope resolution
	if binSize>0:
		xBinWidth=binSize
	else:
		xBinWidth=scaling[3]#YMULT Value
		if xBinWidth<1e-5: #shouldn't be - failsafe
			xBinWidth=0.002
		if fit>-1:
			xBinWidth=0.5 #0.5 keV bins for calibrated data by default
		if integral and fit==-1:
			xBinWidth=1 #[V*ns]
	xMax=np.max(data)
	xMin=np.min(data) - np.min(data)*0.05 if injection else 0
	
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data
	
	#Create histogram of data
	hist=plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	ydata=hist[0]	
	errs=np.sqrt(ydata)
	#avoid division by zero
	for i,err in enumerate(errs):
		if err==0:
			errs[i]=1e-8
	binCenters=hist[1]+xBinWidth/2
	binCenters=binCenters[:-1]

	#Set fit range from inputs - coarse
	low_i,high_i=getArrayIndex(binCenters,fitLow,fitHigh)
	
	#Set up fit with guesses for p01: [Amplitude, Mu, Sigma]
	#muGuess=np.mean(data)
	muGuess_i=np.argmax(ydata[low_i:high_i])+low_i
	muGuess=binCenters[muGuess_i]
	if muGuess>fitHigh or muGuess<fitLow:
		muGuess=binCenters[low_i]+(binCenters[high_i]-binCenters[low_i])/2.
	ampGuess=ydata[low_i:high_i].max()
	if injection and not integral:
		#large range requires dynamic parameter estimation - fit VERY dependent on initial sigma value
		sigGuess=(binCenters[high_i]-muGuess)/(muGuess*10)
	else:
		sigGuess=muGuess/10.
	if edge:
		p01 = [ampGuess, ampGuess, muGuess, sigGuess]
	else:
		p01 = [ampGuess, muGuess, sigGuess]
		
	#Fit histogram over desired range
	if edge:
		popt, pcov = curve_fit(Edge, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], sigma=errs[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000, absolute_sigma=True)
		(Amp1, Amp2, Mu, Sigma)=popt
		#Calculate N (events within 2sigma of mean)
		integ=scipy.integrate.quad(Edge, -2*Sigma, 2*Sigma, args=(Amp1,Amp2,Mu,Sigma))
		#returns integral and uncertainty on calculation - only return integral value
	elif injection: #don't need fancy fitting for Gaussian peaks
		popt, pcov = curve_fit(Gauss, xdata=binCenters[low_i:high_i], ydata=ydata[low_i:high_i], sigma=errs[low_i:high_i], p0=p01, bounds=(0,np.inf), maxfev=5000, absolute_sigma=True)
		(Amp, Mu, Sigma)=popt
	else:		
		#Refine fit range
		if fitLow<muGuess-sigGuess:#if too much low data below peak
			fitLow=muGuess-sigGuess
		if ((fitHigh>muGuess+sigGuess) or (high_i==len(binCenters)-1)):#if too much high data above peak or if bound goes to np.inf (end of dataset)
			fitHigh=muGuess+sigGuess
		#Recalculate indices of bounds
		low_i,high_i=getArrayIndex(binCenters,fitLow,fitHigh)
		#prevent artifically small range if peak on the edge of the distribution
		if high_i-low_i<4:
			high_i+=4
			low_i-=4
			#failsafe against going negative or out of bounds
			if low_i<0:
				low_i=0
			if high_i>len(binCenters)-1:
				high_i=len(binCenters)-1
		popt, pcov = iterativeFit(Gauss, p01, binCenters, ydata, low_i, high_i)
		(Amp, Mu, Sigma)=popt
		#Calculate N (events under fit)
		integ=scipy.integrate.quad(Gauss, -np.inf, np.inf, args=(Amp,Mu,Sigma))
		#returns integral and uncertainty on calculation - only return integral value

	
	#Calculate energy resolution
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	
	#Display fit on final histogram
	plt.rc('text', usetex=True) #use Latex
	xspace=np.linspace(xMin,binEdges[-1],len(binCenters)*10)
	if edge:
		plt.plot(xspace, Edge(xspace, *popt), 'r-', label=r'Fit')  
	else:
		plt.plot(xspace, Gauss(xspace, *popt), 'r-', label=r'Fit')    
	plt.plot([], [], ' ', label=f"Energy res = {enRes:.2f}%")
	plt.plot([], [], ' ', label=r"Fit parameters:")
	if edge:
		plt.plot([], [], ' ', label=f"\small Amp1={Amp1:.2f}, Amp2={Amp2:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	else:
		plt.plot([], [], ' ', label=f"\small Amp={Amp:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.legend(loc="best")
	plt.title(f"Energy Resolution - {title}, pixel {pixel}")
	plt.ylabel('Counts')
	if integral and fit==-1:
		plt.xlabel('Integrated energy [V*ns]')
	elif edge:
		plt.title(f"Edge Fit - {title}, pixel {pixel}")
		plt.xlabel('Peak Energy [V]')
	elif not integral and fit==-1:
		plt.xlabel('Peak Energy [V]')
	else:
		plt.xlabel('Calibrated Energy [keV]')
		plt.title(f"Calibrated {title}, pixel {pixel}")

	#plt.yscale('log')
	if savePlots:
		#save figure
		if fit>-1 and integral:
			saveto=f"{getSaveto(savedir)}{title}{datain}_{energy}line_integralCalibrated.pdf"
		elif fit>-1:
			saveto=f"{getSaveto(savedir)}{title}{datain}_{energy}line_peakCalibrated.pdf"
		elif edge:	
			saveto=f"{getSaveto(savedir)}{title}{datain}EdgeFit_{energy}edge.pdf"
		else:
			saveto=f"{getSaveto(savedir)}{title}{datain}_amp{pixel}_{energy}line.pdf"
		plt.savefig(saveto)
	else:
		plt.show()
	plt.clf()
		
	f.close()

	#return popt, enRes, pcov, integ[0]
	return popt, enRes, pcov
		
	
#Calculate error on fit parameters, save in text file	
def printParams(settings, popt, en_res, pcov, integral=False, edge=False, injection=False, savedir=None):
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
		
	if edge:
		(Amp1, Amp2, Mu, Sigma) = popt
		(Amp1_err, Amp2_err, Mu_err, Sigma_err) = np.sqrt(np.diag(pcov))
	else:
		(Amp, Mu, Sigma) = popt
		(Amp_err, Mu_err, Sigma_err) = np.sqrt(np.diag(pcov))
	# Error propagation
	partial_sigma = (2.355*100)/Mu
	partial_mu = (2.355*100*Sigma)/(Mu**2)
	stdev_er = np.sqrt(((partial_sigma**2)*(Sigma_err**2))+((partial_mu**2)*(Mu_err)**2))
	
	if savePlots:
		if edge:		
			saveto=f"{getSaveto(savedir)}{title}_{energy}edge{datain}EnRes.txt"
			k=open(saveto, "w")
			k.write("Amplitude erfc= %d Amplitude const = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp1, Amp2, Mu, Sigma)+"\n")
		elif injection:
			saveto=f"{getSaveto(savedir)}{title}_amp{pixel}{datain}EnRes.txt"
			k=open(saveto, "w")
			k.write("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		else:
			saveto=f"{getSaveto(savedir)}{title}_{energy}line{datain}EnRes.txt"
			k=open(saveto, "w")
			k.write("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma)+"\n")
		k.write("\n")#keep empty line as place holder so code is backwards compatible - this line used to hold output of N calculation
		k.write("Energy resolution = %0.2f" %(abs(en_res))+"%\n")
		if edge:
			k.write("Error in erfc amplitude = %0.3f \nError in erfc constant = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp1_err, Amp2_err, Mu_err, Sigma_err)+"\n")
		else:
			k.write("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err)+"\n")
		k.write("Error in energy resolution = %0.5f"%(stdev_er)+"%\n")
		k.close()
		#Display contents to terminal
		if edge:
			print(f"FIT FROM {title}_{energy}edge{datain}")
		elif injection:
			print(f"FIT FROM {title}{datain}")
		else:
			print(f"FIT FROM {title}_{energy}line{datain}")
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		if edge:
			print("Amplitude erfc= %d Amplitude const = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp1, Amp2, Mu, Sigma))
			print("Energy resolution = %0.2f" %(abs(en_res)))
			print("Error in erfc amplitude = %0.3f \nError in erfc constant = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp1_err, Amp2_err, Mu_err, Sigma_err))
			print("Error in energy resolution = %0.5f"%(stdev_er))
		else:
			print("Amplitude = %d \nMu = %0.4f \nSigma = %0.4f" %(Amp, Mu, Sigma))
			print("Energy resolution = %0.2f" %(abs(en_res)))
			print("Error in amplitude = %0.3f \nError in mu = %0.6f \nError in sigma = %0.6f" %(Amp_err, Mu_err, Sigma_err))
			print("Error in energy resolution = %0.5f"%(stdev_er))	
		
		
def getVals_fromTxt(inDir, integral=False):		
	energyList, muArr1, sigmaArr1,enResArr1, muErrArr1, sigErrArr1, enresErrArr1 = [],[],[],[],[],[],[]

	if integral:
		dsName="integral"
	else:
		dsName="peaks"

	os.chdir(inDir)
	files = glob.glob('*'+dsName+'*.txt')
	for filename in files:
		plus=0
		energyList.append(float(filename.split('_')[1][:-4]))
		if "edge" in filename:
			plus=1
		openFile=open(filename,'r')
		lines=openFile.readlines()
		muArr1.append(float(lines[1].split(' = ')[-1]))
		sigmaArr1.append(float(lines[2].split(' = ')[-1]))
		enResArr1.append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end
		muErrArr1.append(float(lines[6+plus].split(' = ')[-1]))
		sigErrArr1.append(float(lines[7+plus].split(' = ')[-1]))
		enresErrArr1.append(float(lines[8+plus].split(' = ')[-1][:-2]))#eliminate % sign at the end

	return energyList, muArr1, sigmaArr1, enResArr1, muErrArr1, sigErrArr1, enresErrArr1
	
def getCalibVals_fromTxt(inDir, ele, integral=False):		
	energyList, muArr1, sigmaArr1, enResArr1, muErr, sigErr, enresErr = [],[],[],[],[],[],[]
	fits=['linear','quad','tri','sqrt','spline1','spline3','piecewise']

	if integral:
		dsName="integral"
	else:
		dsName="peaks"

	for fit in fits:
		plus=0
		os.chdir(inDir+fit+'/')
		newfile = glob.glob('*'+ele+'*'+dsName+'*.txt') #returns array with length 1
		try:
			energyList.append(float(newfile[0].split('_')[1][:-4]))
			if "edge" in newfile[0]:
				plus=1
			openFile = open(newfile[0],'r')
			lines=openFile.readlines()
			muArr1.append(float(lines[1].split(' = ')[-1]))
			sigmaArr1.append(float(lines[2].split(' = ')[-1]))
			enResArr1.append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end
			muErr.append(float(lines[6+plus].split(' = ')[-1]))
			sigErr.append(float(lines[7+plus].split(' = ')[-1]))
			enresErr.append(float(lines[8+plus].split(' = ')[-1][:-2]))#eliminate % sign at the end
		except IndexError: #file does not exist so newfile[0] returns index error
			#fill with dummy values
			energyList.append(0)
			muArr1.append(0)
			sigmaArr1.append(0)
			enResArr1.append(0)
			muErr.append(0)
			sigErr.append(0)
			enresErr.append(0)
			
	energy= max(energyList)#if missing elements, array populated with zeros. Want to return the real energy, found from populated values
	
	return energy, muArr1, sigmaArr1, enResArr1, fits, muErr, sigErr, enresErr
		
def getCalibVals_fromFit(inDir, fit, integral=False):	
	#do not consider edge files	
	energyList, muArr1, sigmaArr1, enResArr1, muErr, sigErr, enresErr = [],[],[],[],[],[],[]

	if integral:
		dsName="integralEnRes"
	else:
		dsName="peaksEnRes"

	os.chdir(inDir+fit+'/')
	print(f"{inDir}{fit}")
	newfile = glob.glob('*'+dsName+'*.txt') #returns array with length 1
	for found in newfile:
		print(found)
		if "edge" in found:
			continue
		energyList.append(float(found.split('_')[1][:-4]))
		openFile = open(found,'r')
		lines=openFile.readlines()
		muArr1.append(float(lines[1].split(' = ')[-1]))
		sigmaArr1.append(float(lines[2].split(' = ')[-1]))
		enResArr1.append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end
		muErr.append(float(lines[6].split(' = ')[-1]))
		sigErr.append(float(lines[7].split(' = ')[-1]))
		enresErr.append(float(lines[8].split(' = ')[-1][:-2]))#eliminate % sign at the end
				
	return energyList, muArr1, sigmaArr1, enResArr1, muErr, sigErr, enresErr
		
	
def calcError(sigmaArr, nArr):	
	errArr=[]
	#error Arr1ay = sigma/sqrt(N) (for edge, 2sig integral from mu)
	for j in range(len(sigmaArr)):
		err_p = -1 if nArr[j][0]<=0 else sigmaArr[j][0]/np.sqrt(nArr[j][0])
		err_i = -1 if nArr[j][1]<=0 else sigmaArr[j][1]/np.sqrt(nArr[j][1])
		errArr.append([err_p,err_i])
	return errArr