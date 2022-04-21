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
	
def get_binned_traces (binArr, filename, dataset):
	f = h5py.File(filename, 'r')
	traces = f[dataset]
	events=len(traces)
	peaks=np.array(f[f"{dataset}_peaks"])[:events] #pull first peak values that are associated with the recorded traces

	#return 2d array containing indices of entries for each bin
	indices=[[] for j in range(len(binArr)+1)] 
	for i,peak in enumerate(peaks):
		if peak<binArr[0]:
			indices[0].append(i)
		elif peak<binArr[1]:
			indices[1].append(i)
		else:
			indices[2].append(i)	
	return indices
	
def get_average_trace( filename, dataset, bins):

	if bins:
		binArr=[0.09,0.2]#bin high edges, +max

	f = h5py.File(filename, 'r')
	traces = f[dataset] #baseline subtracted already
	mean=np.mean(traces, axis = 0)

	#return array of average traces
	mean=[mean[0::50]] #reduce number of points for smoother curve
	
	#if breaking into bins, calculate average for each bin
	if bins:
		traces_low=[]
		traces_mid=[]
		traces_high=[]
		
		binIndices=get_binned_traces(binArr,filename,dataset)
		
		for i,trace in enumerate(traces):
			if i in binIndices[0]:
				traces_low.append(trace)
			elif i in binIndices[1]:
				traces_mid.append(trace)
			else:
				traces_high.append(trace)
		
		print(f"Contents of bins \n {len(traces)} total recorded traces. \n Low bin: {len(traces_low)} \n Mid bin: {len(traces_mid)} \n High bin: {len(traces_high)}")
	
		mean_low=np.mean(traces_low,axis=0)
		mean_mid=np.mean(traces_mid,axis=0)
		mean_high=np.mean(traces_high,axis=0)
		
		mean.append(mean_low[0::50])
		mean.append(mean_mid[0::50])
		mean.append(mean_high[0::50])
	
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

def plotBinnedTraces(file, fileOut, ds, bins, strBins, aveTraces):

	f = h5py.File(file, 'r')
	time=get_time(file, ds)
	allTraces = f[ds] #baseline subtracted already
	binIndices=get_binned_traces(bins,file,ds)
	colors=['r','b','g']
	
	binnedTraces=[[] for j in range(len(strBins))]
	for i,trace in enumerate(allTraces):
		if i in binIndices[0]:
			binnedTraces[0].append(trace)
		elif i in binIndices[1]:
			binnedTraces[1].append(trace)
		else:
			binnedTraces[2].append(trace)

	#skip first average trace (contains average of full file)
	for j,trace in enumerate(aveTraces[1:]):
		plt.clf()
		plt.plot(time*1e6, trace, linewidth=3, color="black")
		for eachTrace in binnedTraces[j]:
			plt.plot(time*1e6, eachTrace[0::50], alpha=0.1, color=colors[j]) 
		plt.title(f"Bin {strBins[j]}")
		plt.xlabel( "time ($\mu$s)" )
		plt.ylabel("trace - baseline [V]")
		plot=plt.gcf() #get current figure - saves fig in case savePlt==True
		plt.show() #creates new figure for display
		savePlt=saveFromInput()
		if savePlt: #save plt stored with plt.gcf()
			saveFile=saveDir+fileOut+"_bin"+strBins[j]+".pdf"
			print(f"Saving {saveFile}")
			plot.savefig(saveFile)	
		

###############################
#Runnable choices
###############################

def plotTraces(files,labels,fileOut,ds=["run1"],bins=False):

	i=0
	
	#update labels for legend
	if bins:
		strBinArr=["0-0.09", "0.09-0.2", "0.2+"]
		for index,label in enumerate(labels):
			if index==0:
				continue
				#plot title - not related to legend
			labelArr=[label]
			for binstr in strBinArr:
				labelArr.append(labelArr[0]+"- Bin "+binstr)
			labels[index]=labelArr
			
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

	#plot individual traces along with average
	if bins:
		for f in files:
			file=homeDir+f
			for dataset in ds:
				trace=get_average_trace(file, dataset, bins)
				plotBinnedTraces(file, fileOut, dataset, [0.09,0.2], strBinArr, trace) #bin array is high end of each bin - do not include infinity

		
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
	

#individual traces that feed into average bins
#one plot - average trace binned in peak height [0-0.09,0.09-0.2,0.2+] (trigger at 0.09)
#individual traces on one plot with peak height > 0.2
#do any traces look like expected? (one single peak)		
#histogram of height duration on triggered pulse - full collection and also binned
#compare pulse duration to good injected run

#baseline values?
#redo average plot without subtracting baseline (adding it back in?)- do any traces even have reliable baseline baseline estimates?
#distribution of peak times. peak time vs height, peak time vs duration
		
		
		
		
		
		
		