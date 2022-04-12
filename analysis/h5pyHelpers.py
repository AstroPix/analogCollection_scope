import matplotlib.pyplot as plt
import h5py
import numpy as np
import scipy
from scipy.optimize import curve_fit

#homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/v1/"
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/v2/"



#display histogram to play with
def histDisplay(f_in,ds='run1_peaks', xBinWidth=0.002):
	f=h5py.File(homeDir+f_in, 'r')

	data=f[ds]
	#data=[x*1e-7 for x in data]
	
	#scale=f['run1_scaling']
	#print(np.array(scale))
	xMax=np.max(data)
	xMin=np.min(data)
	
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data

	#Create histogram of data
	plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	#plt.yscale('log')
	plt.show()
	
#investigate histograms of DC baseline for AC runs
#fit functions
def Gauss(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.0*sigma**2))
def histDisplay_baseline(f_in,ds='run1_peaks', xBinWidth=0.002):
	f=h5py.File(homeDir+f_in, 'r')

	data=f[ds]
	peaks=f['run1_peaks']
	
	scale=f['run1_scaling']
	#print(np.array(scale))
	xMax=np.max(data)
	xMin=np.min(data)
	
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data

	xBinWidth2=0.002 #for peaks hist
	binEdges2=np.arange(np.min(peaks),np.max(peaks)+xBinWidth2,xBinWidth2)

	#Create histogram of data
	plt.hist(data,bins=binEdges,color='blue')
	#plt.yscale('log')
	plt.title("Baseline values")
	plt.xlabel("Average DC baseline (first 2000 trace points) [V]")
	plt.show()
	
	plt.clf()
	hist=plt.hist(peaks,bins=binEdges2,label=r'Data', color='blue')
	binCenters=hist[1]+xBinWidth2/2
	binCenters=binCenters[:-1]
	p01=[50,0.23,0.01]
	popt, pcov = curve_fit(Gauss, xdata=binCenters[-15:], ydata=hist[0][-15:], bounds=(0,np.inf), p0=p01, maxfev=8000, absolute_sigma=True)
	(Amp,Mu,Sigma) = popt
	(Amp_err, Mu_err, Sigma_err) = np.sqrt(np.diag(pcov))
	print(Amp_err, Mu_err, Sigma_err)
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	plt.plot([], [], ' ', label=f"Amp={Amp:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.plot([], [], ' ', label=f"Energy Resolution= {enRes:.3f}%")
	xspace=np.linspace(np.min(peaks),binEdges2[-1],len(binCenters)*10)
	plt.plot(xspace, Gauss(xspace, *popt), 'r-', label=r'Fit')    
	#plt.yscale('log')
	plt.title("Peak height")
	plt.xlabel("Peak height [V]")
	plt.legend(loc="best")
	plt.show()
	
	data_lim,peaks_lim=[],[]
	for i,item in enumerate(data):
		if item>0.005:
			data_lim.append(item)
			peaks_lim.append(peaks[i])
	
	plt.clf()
	plt.scatter(peaks_lim,data_lim,label="data")
	plt.axvline(x=np.mean(peaks),color="red",label="Ave. peak height")
	plt.xlabel("peak height [V]")
	plt.ylabel("baseline value [V]")
	plt.legend(loc="best")
	plt.title("Baseline values > 5mV")
	plt.show()

	plt.clf()
	y=np.linspace(0.07,0.25,100)
	lim_sum=[a + b for a, b in zip(data_lim, peaks_lim)]
	plt.scatter(peaks_lim, lim_sum)
	plt.plot(y,y,'r',label="x=y")
	plt.xlabel("peak height [V]")
	plt.ylabel("baseline+peak height [V]")
	plt.legend(loc="best")
	plt.show()
	
	plt.clf()
	tot_sum=[a + b for a, b in zip(data, peaks)]
	hist=plt.hist(tot_sum,bins=binEdges2,label=r'Data', color='blue')
	popt, pcov = curve_fit(Gauss, xdata=binCenters[-15:], ydata=hist[0][-15:], bounds=(0,np.inf), p0=p01, maxfev=8000, absolute_sigma=True)
	(Amp,Mu,Sigma) = popt
	(Amp_err, Mu_err, Sigma_err) = np.sqrt(np.diag(pcov))
	print(Amp_err, Mu_err, Sigma_err)
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	plt.plot([], [], ' ', label=f"Amp={Amp:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.plot([], [], ' ', label=f"Energy Resolution= {enRes:.3f}%")
	plt.plot(xspace, Gauss(xspace, *popt), 'r-', label=r'Fit')    
	#plt.yscale('log')
	plt.title("Peak height + DC baseline")
	plt.xlabel("Peak height+baseline [V]")
	plt.legend(loc="best")
	plt.show()

#combine separate runs into one new file
def combineFiles(files,outFile):

	ds="run1"
	#dsNames=["_t", "_peaks", "_integral", "_baseline", "_scaling"]
	dsNames=["_peaks", "_integral", "_baseline", "_scaling"]


	dataArrays=[]
	for fin in files:
		f=h5py.File(homeDir+fin,'r')
		print(homeDir+fin)
		print(list(f.keys()))
		i=0
		for dsn in dsNames:
			dataArrays.append(np.array([]))
			data_tmp=np.array(f[ds+dsn])
			dataArrays[i]=np.concatenate([dataArrays[i],data_tmp])
			i+=1
		f.close()

	h=h5py.File(homeDir+outFile, 'w')
	j=0
	for dsn in dsNames:
		h.create_dataset(ds+dsn, data=dataArrays[j])
		j+=1
	print(list(h.keys()))
	h.close()
	


#copy scaling dataset into main dataset
def copyScalingDS(f_in_scale, f_in):

	f=h5py.File(homeDir+f_in_scale, 'r')

	my_array=f['run1_scaling']

	g = h5py.File(homeDir+f_in, 'a')
	g.create_dataset('run1_scaling', data=my_array)
	print(list(g.keys()))
	f.close()
	g.close()

#fit Gaussian to a distribution, display both on plot, save if desired
def fitGaussian(f_in,traceInteg,savePlots):
	#binsize=0.0025
	binsize=0.001
	f = h5py.File(homeDir+f_in,'r')
	data=f["run1_integral"] if traceInteg else f["run1_peaks"]
	bins=np.arange(min(data),max(data)+binsize,binsize)
	hist=plt.hist(data,bins=bins,label=r'Data', color='blue')
	ydata=hist[0]
	errs=np.sqrt(ydata)
	#avoid division by zero
	for i,err in enumerate(errs):
		if err==0:
			errs[i]=1e-8
	binCenters=hist[1]+(binsize/2)
	binCenters=binCenters[:-1]
	
	muGuess=np.mean(data)
	sigGuess=muGuess/100.
	ampGuess=ydata.max()
	p01 = [ampGuess, muGuess, sigGuess]
	print(p01)

	popt, pcov = curve_fit(Gauss, xdata=binCenters, ydata=ydata, sigma=errs, p0=p01, bounds=(0,np.inf), maxfev=5000, absolute_sigma=True)
	(Amp, Mu, Sigma)=popt
	enRes = (2.355*abs(Sigma)*100)/Mu # Calculates energy resolution, 2.355 converts sigma to FWHM
	print(Sigma)

	#Display fit on final histogram
	plt.rc('text', usetex=True) #use Latex
	xspace=np.linspace(bins[0],bins[-1],len(binCenters)*10)
	plt.plot(xspace, Gauss(xspace, *popt), 'r-', label=r'Fit')    
	plt.plot([], [], ' ', label=f"Energy res = {enRes:.2f}%")
	plt.plot([], [], ' ', label=r"Fit parameters:")
	plt.plot([], [], ' ', label=f"\small Amp={Amp:.2f}, $\mu$={Mu:.3f}, $\sigma$={Sigma:.3f}")
	plt.legend(loc="best")
	plt.ylabel('Counts')
	if traceInteg:
		plt.xlabel('Integrated energy [V*ns]')
	else:
		plt.xlabel('Peak Energy [V]')
		
	if savePlots:
		plt.savefig(homeDir+"gaussianFit.png")
	else:
		plt.show()



#################################################################
# main
#################################################################

if __name__ == "__main__":
	
	savePlots=False
	traceInteg=False
	f="110421_amp1/Americium_480min_combined.h5py"
	filesIn=["102921_amp1/americium241_90min.h5py", "110421_amp1/Americium_120min.h5py","110421_amp1/test_Americium_30min.h5py","110421_amp1/_Americium_240min.h5py"]
	outFile="110421_amp1/Americium_480min_combined.h5py"

	filesInCd = ["120721_amp1/highPeak_chip004_AC_cobalt57_960min.h5py","121521_amp1/150mV_chip004_cobalt57_1020min.h5py"]	
	#todel="011022_amp2/chip004_cobalt57_combined_2220min.h5py"
	filesInCo = ["121421_amp2/150mV_chip004_cobalt57_1200min.h5py","011022_amp2/chip004_cobalt57_1020min.h5py"]
	#todel="121521_amp1/150mV_chip004_cobalt57_combined_2040min.h5py"
	filesInv2Co = ["030422_amp1/chip2_75mV_cobalt57_120min.h5py","030422_amp1/chip2_75mV_cobalt57_120min.h5py","030422_amp1/chip2_150mV_cobalt57_150min.h5py"]
	#todel="030422_amp1/chip2_varTrig_cobalt57_300min_combined.h5py"
	filesInv2Cd=["030322_amp1/TEST_chip2_200mV_cadmium109_150min.h5py","030822_amp1/chip2_200mV_cadmium109_180min.h5py"]
	todel="030822_amp1/chip2_200mV_cadmium109_330min_combined.h5py"
	
	#histDisplay("121521_amp1/150mV_chip004_cobalt57_combined_2040min.h5py",xBinWidth=0.002)
	#histDisplay_baseline("120921_amp1/90mV_chip004_AC_Cadmium_1200min.h5py",ds='run1_baseline', xBinWidth=0.00005)
	#histDisplay("120221_amp2/calib_cadmium190_1080min.h5py",ds='run1_baseline', xBinWidth=0.00005)
	#combineFiles(filesInv2Cd, todel)
	#copyScalingDS("110821_amp1/test_barium133_30min.h5py",todel)
	#fitGaussian("030822_amp1/chip2_64mV_background_120min.h5py",traceInteg,savePlots)
		
	histDisplay("030822_amp1/chip2_200mV_cadmium109_330min_combined.h5py")
		
	fileName_eb2=["032922_amp1/chip2_baseline_1.0Vinj_2min.h5py","032922_amp1/chip2_EBt_0.5in_1.0Vinj_2min.h5py","032922_amp1/chip2_EBt_1.7in_1.0Vinj_2min.h5py","032922_amp1/chip2_EBt_3.2in_1.0Vinj_2min.h5py"]
	fileName_eb3=["032922_amp1/chip3_baseline_1.0Vinj_2min.h5py","032922_amp1/chip3_EBt_0.5in_1.0Vinj_2min.h5py","032922_amp1/chip3_EBt_1.7in_1.0Vinj_2min.h5py","032922_amp1/chip3_EBt_3.2in_1.0Vinj_2min.h5py"]
	histDisplay("032922_amp1/chip3_baseline_1.0Vinj_2min.h5py",xBinWidth=0.002)
	fitGaussian("032922_amp1/chip3_baseline_1.0Vinj_2min.h5py",traceInteg,savePlots)		

	for entry in fileName_eb2:
		fitGaussian(entry,traceInteg,savePlots)		
	for entry in fileName_eb3:
		fitGaussian(entry,traceInteg,savePlots)				
		
		


