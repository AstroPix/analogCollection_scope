import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit

###############################
#Global variables
###############################

chip=3 #choose 2 or 3
homeDir = f"/Users/asteinhe/AstroPixData/astropixOut_tmp/v2/"
saveDir = f"/Users/asteinhe/AstroPixData/astropixOut_tmp/heavyIonBeam/chip{chip}"

###############################
#helper functions
###############################

def get_time(filename, dataset):	

	f = h5py.File(filename, 'r')
	time = np.array(f[dataset+"_t"])
	#reduce number of points for smoother curve
	time=time[0::50]
	
	return time
	
def get_average_trace( filename, dataset, bins ):

	if bins:
		binArr=[0.09,0.2]#bin high edges, +max

	f = h5py.File(filename, 'r')
	traces = f[dataset] #baseline subtracted already
	mean=np.mean(traces, axis = 0)

	#return array of average traces
	mean=[mean[0::50]] #reduce number of points for smoother curve
	
	#if breaking into bins, calculate average for each bin
	if bins:
		events=len(traces)
		peaks=np.array(f[f"{dataset}_peaks"])[:events] #pull first peak values that are associated with the recorded traces
		traces_low=[]
		traces_mid=[]
		traces_high=[]
		
		test_i=0
		for trace,peak in zip(traces,peaks):
			if peak<binArr[0]:
				traces_low.append(trace)
			elif peak<binArr[1]:
				traces_mid.append(trace)
			else:
				traces_high.append(trace)
				print(f"High trace: element {test_i}")
			test_i+=1
				
		print(f"Contents of bins \n {len(traces)} total recorded traces. \n Low bin: {len(traces_low)} \n Mid bin: {len(traces_mid)} \n High bin: {len(traces_high)}")
	
		mean_low=np.mean(traces_low,axis=0)
		mean_mid=np.mean(traces_mid,axis=0)
		mean_high=np.mean(traces_high,axis=0)
		
		mean.append(mean_low[0::50])
		mean.append(mean_mid[0::50])
		mean.append(mean_high[0::50])
		
		mean.append(traces[59][0::50])
	
	return mean

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
###############################

def plotTraces(files,labels,fileOut,ds=["run1"],bins=False):
	i=0
	
	#update labels for legend
	if bins:
		binArr=["0-0.09", "0.09-0.2", "0.2+"]
		for index,label in enumerate(labels):
			if index==0:
				continue
				#plot title - not related to legend
			labelArr=[label]
			for binstr in binArr:
				labelArr.append(labelArr[0]+"- Bin "+binstr)
			labels[index]=labelArr
			
			labels[1].append("trigger before tall one")
			
	#for each file, calculate average trace (also calculate in bins, if bool is true)
	for f in files:
		file=homeDir+f
		j=0
		while j<len(ds):
			time=get_time(file,ds[j])
			trace=get_average_trace(file, ds[j], bins)
			i+=1
			for k,aveTrace in enumerate(trace):
				plt.plot(time*1e6, aveTrace, label=labels[i][k])
			j+=1
	#plot average traces together
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
def plotTraces_compare(fileName,labels, fileOut, ds=["run1"], xrange=0.004, ratioBool=True, versionBool=False, versionArray=[]):
	traces=[]
	ratio=None
	noiseRMS=[]
	if versionBool:
		global saveDir
		saveDir="/Users/asteinhe/AstroPixData/astropixOut_tmp/noise_gain_thresholdScan/"
	if ratioBool:
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [2,1]})
	else:
		#cheat so that ratio plot is not displayed in what is technically a second subfigure
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [100,1]})
	fig.suptitle(labels[0])
	i=0
	if versionBool and len(versionArray)==0:
		print("ERROR - Must input array of AstroPix version associated with input files if versionBool=True")
		return
	for v,f in enumerate(fileName):
		if versionBool:
			global homeDir
			homeDir=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/v{versionArray[v]}/"
		file=homeDir+f
		print(file)
		for set in ds:
			time=get_time(file,set)
			trace=get_average_trace(file, set)
			binCenters, hist= get_baseline_plt(file, set, xrange)
			noiseRMS.append(binCenters)
			noiseRMS.append(hist)
			i+=1
			axs[0].plot(time*1e6, trace, label=labels[i])
			traces.append(trace)
	ratio=traces[0]/traces[-1]
	axs[0].legend(loc="best")
	plt.xlabel( "time ($\mu$s)" )
	plt.setp(axs[0], ylabel="trace - baseline [V]")
	if ratioBool:
		axs[1].plot(time*1e6,ratio)
		plt.setp(axs[1],ylabel=f"Ratio {labels[1]}/{labels[-1]}")
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
		if i%2!=0: #Find end of plotting pair [binCenters,hist]
			plt.plot(noiseRMS[i-1],noiseRMS[i],label=labels[labelIndex])
			labelIndex+=1
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

		
		
###############################
#Main
###############################

if __name__ == "__main__":

	##Plot multiple average traces on one plot
	##Required arguments: input files, legend labels [first entry = plot title], output file name
	##Optional argument: array of datasets to compare (must be in same file: ex. run1 and and run2 from same input file). Default = only "run1"
	
	filesIn=["040722_amp1/LBNL_inCave_withBeam_ion1_chip3_EBt_1.8Vinj_2min.h5py"]
	labels=["1.8V Injection in HI beam","Chip003, amp1"]
	plotTraces(filesIn,labels,"v2_1.8Vinj_chip3_HIbeam_ion1",bins=True)
	

	##Compare two traces by plotting together, plotting ratio (entry 0 / entry 1), and calculating STD of noise from each and plotting
	##Required arguments: input files, labels, outFile
	##Optional arguments: array of datasets [default = 'run1'], xrange max [default = 0.004] - for RMS, ratioBool bool to display ratio plot [default=True] between first and last input files, versionBool bool to compare versions [default=False] - must set desired version in global vars
	"""
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
		#plotTraces_compare(fileName_scan, trigScan, f"102221_amp{p}_0.05Vinj_scaleScan", xrange=0.03,ratioBool=False)	
	"""

#individual traces that feed into average bins
#one plot - average trace binned in peak height [0-0.09,0.09-0.2,0.2+] (trigger at 0.09)
#individual traces on one plot with peak height > 0.2
#do any traces look like expected? (one single peak)		
#histogram of height duration on triggered pulse - full collection and also binned
#compare pulse duration to good injected run

#baseline values?
#redo average plot without subtracting baseline (adding it back in?)- do any traces even have reliable baseline baseline estimates?
#distribution of peak times. peak time vs height, peak time vs duration
		
		
		
		
		
		
		