import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
import sys,os



# Define the fit function (a Gaussian)
def Gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))


#Plot data, fit to Gaussian, calculate energy resolution and draw plot
def enResPlot(file, integral=False):
	#Distinguish between peaks and integral
	if (integral):
		datain='_integral'
	else: #peaks
		datain='_peaks'

	#get info from file
	print(file)
	f = h5py.File(file,'r')
	dsName='run1'
	time=f[dsName+'_t']
	data=f[dsName+datain]
	scaling=f[dsName+'_scaling']
	
	#Create arrays for binning based on scope resolution
	xBinWidth=scaling[3]#YMULT Value
	xMax=np.max(data)
	if integral:
		xMin=np.min(data)
	else:
		xMin=0
	if integral:
		xBinWidth*=1000
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data
	binCenters=binEdges+xBinWidth/2
	binCenters=binCenters[:-1]
	
	#Create histogram of data
	hist=plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	
	#Set up fit with guesses for p01: [Amplitude, Mu, Sigma]
	muGuess=np.mean(data)
	ampGuess=hist[0].max()
	p01 = [ampGuess, muGuess, muGuess/2.]
	
	#Fit histogram
	popt, pcov = curve_fit(Gauss, xdata=binCenters, ydata=hist[0], p0=p01, maxfev=5000)
	(Amp,Mu,Sigma)=popt
	
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
	plt.title(f"Energy Resolution - {inj}V injection, pixel {pixel}")
	
	#save figure
	saveto=file[:-5]
	saveto=f"{saveto}{datain}EnRes{enRes:.2f}.pdf"
	print(saveto)
	plt.savefig(saveto)
	plt.clf()
	
	f.close()

	return popt, enRes, pcov
	

#Calculate error on fit parameters, save in text file	
def printParams(file, popt, en_res, pcov, integral=False):
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



def injScanPlot(inj, data, dataName,saveto):
	
	dataNameStr=dataName.replace(" ", "")
	
	amp1_p=[]
	amp2_p=[]
	amp1_i=[]
	amp2_i=[]
	
	for injEnergy in data:
		flatData=np.array(injEnergy).flatten('F')
		amp1_p.append(flatData[0])
		amp2_p.append(flatData[1])
		amp1_i.append(flatData[2])
		amp2_i.append(flatData[3])
		
	plt.plot(inj,amp1_p,label='Amp1',marker="o")
	plt.plot(inj,amp2_p,label='Amp2',marker="o")
	plt.xlabel("Injection voltage [V]")
	plt.ylabel(f"{dataName} (from peak height)")
	plt.legend(loc="best")
	plt.yscale('log')
	#plt.show()
	plt.savefig(f"{saveto}_peaks_{dataNameStr}_log.pdf")
	plt.clf()
	
	plt.plot(inj,amp1_i,label='Amp1',marker="o")
	plt.plot(inj,amp2_i,label='Amp2',marker="o")
	plt.xlabel("Injection voltage [V]")
	plt.ylabel(f"{dataName} (from integral)")
	plt.legend(loc="best")
	plt.yscale('log')
	#plt.show()
	plt.savefig(f"{saveto}_integral_{dataNameStr}_log.pdf")
	plt.clf()
	
	

injection=[i*0.1 for i in range(1,19)]
enResArr=[]
muArr=[]
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp"
enRes_tmp, mu_tmp=[], []
for inj in injection:
	for pixel in [1, 2]:
		file=f"{homeDir}/102221_amp{pixel}/{inj:.1f}Vinj.h5py"
		popt, enRes, pcov = enResPlot(file)
		#printParams(file, popt, enRes, pcov)
		poptI, enResI, pcovI = enResPlot(file, integral=True)
		#printParams(file, poptI, enResI, pcovI, integral=True)
		enRes_tmp.append([enRes, enResI])
		mu_tmp.append([popt[1], poptI[1]])
	enResArr.append(enRes_tmp)
	muArr.append(mu_tmp)
	enRes_tmp,mu_tmp=[],[]
		
		
		
injScanPlot(injection, enResArr, "Energy Resolution [\%]", homeDir+"/102221")	
injScanPlot(injection, muArr, "Fit Mean [V]",homeDir+"/102221")		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		