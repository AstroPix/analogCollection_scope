import matplotlib.pyplot as plt
import h5py
import numpy as np
import scipy
from scipy.optimize import curve_fit


homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/"



#display histogram to play with
def histDisplay(f_in,ds='run1_peaks', xBinWidth=0.002):
	f=h5py.File(homeDir+f_in, 'r')

	data=f[ds]
	#data=[x*1e-7 for x in data]
	
	scale=f['run1_scaling']
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
	#dsNames=["_t", "_peaks", "_integral", "_baseline", "_trigTime"]
	dsNames=["_t", "_peaks", "_integral", "_baseline"]

	dataArrays=[]
	print(dataArrays)
	#datap=np.array([])
	#datai=np.array([])
	#datas=np.array([])
	for fin in files:
		f=h5py.File(homeDir+fin,'r')
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



#################################################################
# main
#################################################################

if __name__ == "__main__":
	
	f="110421_amp1/Americium_480min_combined.h5py"
	filesIn=["102921_amp1/americium241_90min.h5py", "110421_amp1/Americium_120min.h5py","110421_amp1/test_Americium_30min.h5py","110421_amp1/_Americium_240min.h5py"]
	outFile="110421_amp1/Americium_480min_combined.h5py"
	
	todel="110821_amp1/barium133_combined_65min.h5py"
	filesInBa = ["110821_amp1/test_barium133_30min.h5py","110821_amp1/test_barium133_35min.h5py"]
	
	#histDisplay("120721_amp1/lowPeak_chip004_AC_cobalt57_60min.h5py",xBinWidth=0.002)
	histDisplay_baseline("120921_amp1/90mV_chip004_AC_Cadmium_1200min.h5py",ds='run1_baseline', xBinWidth=0.00005)
	#histDisplay("120221_amp2/calib_cadmium190_1080min.h5py",ds='run1_baseline', xBinWidth=0.00005)
	#combineFiles(filesInBa, todel)
	#copyScalingDS("110821_amp1/test_barium133_30min.h5py",todel)

		
		
		
		
		


