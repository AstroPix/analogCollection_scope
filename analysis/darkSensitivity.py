import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit


def get_time(filename, dataset):	

	f = h5py.File(filename, 'r')
	time = np.array(f[dataset+"_t"])
	#reduce number of points for smoother curve
	time=time[0::50]
	
	return time
	
def get_average_trace( filename, dataset ):

	f = h5py.File(filename, 'r')
	traces = f[dataset] #baseline subtracted already

	mean=np.mean(traces, axis = 0)
	#reduce number of points for smoother curve
	mean=mean[0::50]
	
	return mean

def get_baseline_plt(filename, dataset, xrange=0.004):

	f = h5py.File(filename, 'r')
	baseline = f[dataset]
	baseline=baseline[:,0:2000]
	rmsArr=[np.std(trace) for trace in baseline] 
	
	
	xspace=np.linspace(0,xrange,50)
	hist,binEdges=np.histogram(rmsArr, bins=xspace)
	binWidth=binEdges[1]-binEdges[0]
	binCenters=binEdges+binWidth
	binCenters=binCenters[:-1]

	return binCenters, hist


###############################
#Runnable choices
def plotTraces():
	fileName=["","_dark"]
	homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp"
	for name in fileName:
		for pixel in [1, 2]:
			file=f"{homeDir}/102221_amp{pixel}/0.05Vinj{name}.h5py"
			print(file)
			time=get_time(file,'run1')
			gain1_trace=get_average_trace(file, 'run1')
			plt.plot(time*1e6, gain1_trace, label=f"Amp{pixel}, gain 1")
			time2=get_time(file,'run2')
			gain2_trace=get_average_trace(file, 'run2')
			plt.plot(time*1e6, gain2_trace, label=f"Amp{pixel}, gain 2")
		plt.legend(loc="best")
		plt.xlabel( "time ($\mu$s)" )
		plt.ylabel("trace - baseline [V]")
		plt.title(f"0.05V Injection {name}")
		plt.savefig(f"{homeDir}/102221_0.05Vinj{name}_traces.pdf")
		plt.clf()
		

		
#Compare two traces by plotting together, plotting ratio, and calculating STD of noise from each and plotting
def plotTraces_compare(fileName, dataset,labels, xrange=0.004, ratioBool=True):
	for pixel in [1, 2]:
		traces=[]
		ratio=None
		noiseRMS=[]
		if ratioBool:
			fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [2,1]})
		else:
			#cheat so that ratio plot is not displayed in what is technically a second subfigure
			fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [100,1]})
		fig.suptitle(f"0.05V Injection, pixel {pixel}")
		for name in fileName:
			for set in dataset:
				file=f"{homeDir}/102221_amp{pixel}/0.05Vinj{name}.h5py"
				print(file)
				time=get_time(file,'run1')
				trace=get_average_trace(file, set)
				binCenters, hist= get_baseline_plt(file, set, xrange)
				noiseRMS.append(binCenters)
				noiseRMS.append(hist)
				axs[0].plot(time*1e6, trace, label=f"Amp{pixel} {name}")
				traces.append(trace)
		ratio=traces[0]/traces[1]
		axs[0].legend(loc="best")
		plt.xlabel( "time ($\mu$s)" )
		plt.setp(axs[0], ylabel="trace - baseline [V]")
		if ratioBool:
			axs[1].plot(time*1e6,ratio)
			plt.setp(axs[1],ylabel=f"Ratio {labels[1]}/{labels[2]}")
		#plt.show()
		plt.savefig(f"{homeDir}/102221_amp{pixel}_0.05Vinj_traces_compare{labels[0]}.pdf")
		plt.clf()
	
		labelIndex=1
		for i in range(len(noiseRMS)): 
			if i%2!=0: #Find end of plotting pair
				plt.plot(noiseRMS[i-1],noiseRMS[i],label=labels[labelIndex])
				labelIndex+=1
		#plt.plot(noiseRMS[2],noiseRMS[3],label=labels[2])
		plt.title(f"Noise RMS pixel{pixel}")
		plt.xlabel("RMS [V]")
		plt.ylabel("counts")
		plt.legend(loc="best")
		#plt.show()
		plt.savefig(f"{homeDir}/102221_amp{pixel}_0.05Vinj_noiseRMS_{labels[0]}.pdf")
		plt.clf()
			
		
		
		
		
		
		
		
if __name__ == "__main__":
	fileName=["","_dark"]
	trigScanVal=[5,10,20,50,100,200,500,1000]
	fileName_scan=["_trigScan_"+str(i)+"mV" for i in trigScanVal]
	dataset=["run1","run2"]
	dark=["","Light","Dark"]
	gain=["Gain","gain1","gain2"]
	trigScan=[str(i)+"mV" for i in trigScanVal]
	trigScan.insert(0,"scaleScan")
	homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp"
	
	plotTraces()		
	plotTraces_compare([""],dataset,gain)	
	plotTraces_compare(fileName,["run1"],dark)	
	plotTraces_compare(fileName_scan,["run1"],trigScan,xrange=0.03,ratioBool=False)	
		
		
		
		
		
		
		
		
		
		
		
		