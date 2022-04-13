import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit

###############################
#Global variables
vers=1
homeDir = f"/Users/asteinhe/AstroPixData/astropixOut_tmp/v{vers}/"
saveDir= f"/Users/asteinhe/AstroPixData/astropixOut_tmp/noise_gain_thresholdScan/v{vers}/"

###############################
#helper functions
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
	
def first_nonzero(lst):
	for i, value in enumerate(lst):
		if value > 0.02:
			return i
	return -1

def last_nonzero(lst):
	for i, value in enumerate(reversed(lst)):
		if value > 0.02:
			return len(lst)-i-1
	return -1

	
def get_height_duration(filename, dataset):

	f = h5py.File(filename, 'r')
	traces = np.array(f[dataset]) #baseline subtracted already
	time = np.array(f[dataset+"_t"])
	
	toHist=[]
	for trace in traces:
		height=np.max(trace)
		duration_start=first_nonzero(trace)
		duration_end=last_nonzero(trace)
		duration=time[duration_end]-time[duration_start]
		if duration==0:
			duration=1e-8
		toHist.append(height/duration)

	return toHist
	
def saveFromInput():
	inp= 'a'
	while inp!='y' and inp!='n':
		print("Save plot? y/n")
		inp = input()
	
	result=True if inp=='y' else False
	return result

###############################
#Runnable choices
def plotTraces(files,labels,fileOut,ds=["run1"]):
	i=0
	for f in files:
		file=homeDir+f
		j=0
		while j<len(ds):
			time=get_time(file,ds[j])
			trace=get_average_trace(file, ds[j])
			i+=1
			plt.plot(time*1e6, trace, label=labels[i])
			j+=1
	plt.legend(loc="best")
	plt.xlabel( "time ($\mu$s)" )
	plt.ylabel("trace - baseline [V]")
	plt.title(labels[0])
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	savePlt=saveFromInput()
	if savePlt: #save plt stored with plt.gcf()
		saveFile=saveDir+fileOut+".pdf"
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)	
	plt.clf()

		
#Compare two traces by plotting together, plotting ratio, and calculating STD of noise from each and plotting
def plotTraces_compare(fileName,labels, fileOut, ds=["run1"], xrange=0.004, ratioBool=True):
	traces=[]
	ratio=None
	noiseRMS=[]
	if ratioBool:
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [2,1]})
	else:
		#cheat so that ratio plot is not displayed in what is technically a second subfigure
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [100,1]})
	fig.suptitle(labels[0])
	i=0
	for f in fileName:
		file=homeDir+f
		for set in ds:
			time=get_time(file,set)
			trace=get_average_trace(file, set)
			binCenters, hist= get_baseline_plt(file, set, xrange)
			noiseRMS.append(binCenters)
			noiseRMS.append(hist)
			i+=1
			axs[0].plot(time*1e6, trace, label=labels[i])
			traces.append(trace)
	ratio=traces[0]/traces[1]
	axs[0].legend(loc="best")
	plt.xlabel( "time ($\mu$s)" )
	plt.setp(axs[0], ylabel="trace - baseline [V]")
	if ratioBool:
		axs[1].plot(time*1e6,ratio)
		plt.setp(axs[1],ylabel=f"Ratio {labels[1]}/{labels[2]}")
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	savePlt=saveFromInput()
	if savePlt: #save plt stored with plt.gcf()
		saveFile=saveDir+fileOut+".pdf"
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)	
	plt.clf()

	labelIndex=1
	for i in range(len(noiseRMS)): 
		if i%2!=0: #Find end of plotting pair
			plt.plot(noiseRMS[i-1],noiseRMS[i],label=labels[labelIndex])
			labelIndex+=1
	#plt.plot(noiseRMS[2],noiseRMS[3],label=labels[2])
	plt.title(f"Noise RMS - {labels[0]}")
	plt.xlabel("RMS [V]")
	plt.ylabel("counts")
	plt.legend(loc="best")
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	savePlt=saveFromInput()
	if savePlt: #save plt stored with plt.gcf()
		saveFile=saveDir+fileOut+"_noiseRMS.pdf"
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)	
	plt.clf()
			
#Compare two traces by plotting together, plotting ratio, and calculating correlation of height and duration
def plotTraces_compareVersions(fileName,labels, xrange=0.004, ratioBool=True):
	traces=[]
	ratioHist=[]
	ratio=None
	noiseRMS=[]
	dataset="run1"
	if ratioBool:
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [2,1]})
	else:
		#cheat so that ratio plot is not displayed in what is technically a second subfigure
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [100,1]})
	fig.suptitle(f"1.0V Injection, v1 Chip003 pixel 1 vs v2")
	for i,name in enumerate(fileName):
		file=homeDir+name
		time=get_time(file,dataset)
		trace=get_average_trace(file, dataset)
		binCenters, hist= get_baseline_plt(file, dataset, xrange)
		ratioHist_tmp = get_height_duration(file,dataset)
		noiseRMS.append(binCenters)
		noiseRMS.append(hist)
		axs[0].plot(time*1e6, trace, label=f"{labels[i]}")
		traces.append(trace)
		ratioHist.append(ratioHist_tmp)
	ratio=traces[0]/traces[1]
	axs[0].legend(loc="best")
	plt.xlabel( "time ($\mu$s)" )
	plt.setp(axs[0], ylabel="trace - baseline [V]")
	if ratioBool:
		axs[1].plot(time*1e6,ratio)
		plt.setp(axs[1],ylabel=f"Ratio {labels[0]}/{labels[1]}")
	plt.show()
	#plt.savefig(f"{homeDir}/noise_gain_threshold/102221_amp{pixel}_0.05Vinj_traces_compare{labels[0]}.pdf")
	plt.clf()

	#plot RMS
	labelIndex=0
	for i in range(len(noiseRMS)): 
		if i%2!=0: #Find end of plotting pair
			plt.plot(noiseRMS[i-1],noiseRMS[i],label=labels[labelIndex])
			labelIndex+=1
	#plt.plot(noiseRMS[2],noiseRMS[3],label=labels[2])
	plt.title(f"Noise RMS")
	plt.xlabel("RMS [V]")
	plt.ylabel("counts")
	plt.legend(loc="best")
	plt.show()
	#plt.savefig(f"{homeDir}/noise_gain_threshold/102221_amp{pixel}_0.05Vinj_noiseRMS_{labels[0]}.pdf")
	plt.clf()	
		
	#plot ratio of height to duration	
	h1,bins=np.histogram(ratioHist[0])
	binWidth=(bins[1]-bins[0])/2
	bins1=bins+binWidth
	plt.plot(bins1[:-1],h1,label="v1")
	h2,bins=np.histogram(ratioHist[1])
	binWidth=(bins[1]-bins[0])/2
	bins2=bins+binWidth
	plt.plot(bins2[:-1],h2,label="v2")
	plt.xscale('log')
	plt.xlabel("height/duration")
	plt.ylabel("counts")
	plt.legend(loc="best")
	plt.show()
		
		
###############################
#Main
if __name__ == "__main__":

	##Plot multiple average traces on one plot
	##Required arguments: input files, legend labels, output file name
	##Optional argument: array of datasets to compare (must be in same file: ex. run1 and and run2 from same input file). Default = only "run1"
	#filesIn=["102221_amp1/0.05Vinj.h5py","102221_amp2/0.05Vinj.h5py"]
	#labels=["0.05V Injection","Amp1, gain1", "Amp1, gain2", "Amp2, gain1", "Amp2, gain2"]
	#plotTraces(filesIn,labels,"102221_0.05Vinj_traces",ds=['run1','run2'])
	#filesIn=["102221_amp1/0.05Vinj_dark.h5py","102221_amp2/0.05Vinj_dark.h5py"]
	#labels=["0.05V Injection dark","Amp1, gain1", "Amp1, gain2", "Amp2, gain1", "Amp2, gain2"]
	#plotTraces(filesIn,labels,"102221_0.05Vinj_dark_traces",ds=['run1','run2'])

	##Compare two traces by plotting together, plotting ratio (entry 0 / entry 1), and calculating STD of noise from each and plotting
	##Required arguments: input files, labels, outFile
	##Optional arguments: array of datasets [default = 'run1'], xrange max [default = 0.004], bool to display ratio plot [default=True]
	pix=[1,2]
	trigScanVal=[5,10,20,50,100,200,500,1000]
	trigScan=[str(i)+"mV" for i in trigScanVal]
	for p in pix:
		filesIn=[f"102221_amp{p}/0.05Vinj.h5py",f"102221_amp{p}/0.05Vinj_dark.h5py"]
		dark=[f"0.05V Injection Amp{p}","Light","Dark"]
		#plotTraces_compare(filesIn, dark, f"102221_amp{p}_0.05Vinj_traces_compare")	
		
		gain=[f"0.05V Injection (Light) Amp{p}","gain1","gain2"]
		#plotTraces_compare([f"102221_amp{p}/0.05Vinj.h5py"], gain,f"102221_amp{p}_0.05Vinj_traces_compareGain" ,ds=['run1','run2'])	

		fileName_scan=[f"102221_amp{p}/0.05Vinj_trigScan_"+str(i)+"mV.h5py" for i in trigScanVal]
		trigScan.insert(0,f"Trigger Scan, Amp{p}")
		plotTraces_compare(fileName_scan, trigScan, f"102221_amp{p}_0.05Vinj_scaleScan", xrange=0.03,ratioBool=False)	

	
	versions=["v1","v2"]
	fileName_versions=["/v1/102221_amp1/1.0Vinj.h5py","/v2/030122_amp1/scan_1.0Vinj_2min.h5py"]	
	#plotTraces_compareVersions(fileName_versions,versions,xrange=0.01)	
		
		
		
		
		
		
		
		
		
		
		