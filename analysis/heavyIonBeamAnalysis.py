import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from datetime import datetime

###############################
#Global variables
###############################

chip=3 #choose 2 or 3
homeDir = f"/Users/asteinhe/AstroPixData/astropixOut_tmp/v2/"
saveDir = f"/Users/asteinhe/AstroPixData/astropixOut_tmp/heavyIonBeam/chip{chip}/"

###############################
#helper functions
###############################

def get_time(filename, dataset):	

	f = h5py.File(filename, 'r')
	time = np.array(f[dataset+"_t"])
	#reduce number of points for smoother curve
	time=time[0::50]
	
	return time

###############################
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

###############################	
def get_average_trace( filename, dataset, bins, addBaseline):

	if bins:
		binArr=[0.09,0.2]#bin high edges, +max

	f = h5py.File(filename, 'r')
	traces = f[dataset] #baseline subtracted already
	events=len(traces)

	if addBaseline:
		baseline = f[dataset+"_baseline"][:events]
	else:
		baseline = [0 for i in range(events)]
		
	scaledTraces=[]
	for i,trace in enumerate(traces):
		scaledTraces.append(trace+baseline[i])
	
	mean=np.mean(scaledTraces, axis = 0)

	#return array of average scaledTraces
	mean=[mean[0::50]] #reduce number of points for smoother curve
	
	#if breaking into bins, calculate average for each bin
	if bins:
		scaledTraces_low=[]
		scaledTraces_mid=[]
		scaledTraces_high=[]
		
		binIndices=get_binned_traces(binArr,filename,dataset)
		
		for i,trace in enumerate(scaledTraces):
			if i in binIndices[0]:
				scaledTraces_low.append(trace)
			elif i in binIndices[1]:
				scaledTraces_mid.append(trace)
			else:
				scaledTraces_high.append(trace)
		
		print(f"Contents of bins \n {len(scaledTraces)} total recorded scaledTraces. \n Low bin: {len(scaledTraces_low)} \n Mid bin: {len(scaledTraces_mid)} \n High bin: {len(scaledTraces_high)}")
	
		try:
			mean_low=np.mean(scaledTraces_low,axis=0)
			mean.append(mean_low[0::50])
		except IndexError: #if array is empty
			mean.append(np.zeros(200))
		try:
			mean_mid=np.mean(scaledTraces_mid,axis=0)
			mean.append(mean_mid[0::50])
		except IndexError: #if array is empty
			mean.append(np.zeros(200))
		try:
			mean_high=np.mean(scaledTraces_high,axis=0)
			mean.append(mean_high[0::50])
		except IndexError: #if array is empty
			mean.append(np.zeros(200))
	
	return mean

###############################
def first_nonzero(lst):
	for i, value in enumerate(lst):
		if value > 0.01:
			return i
	return -1

###############################
def last_nonzero(lst):
	for i, value in enumerate(reversed(lst)):
		if value > 0.01:
			return len(lst)-i-1
	return -1
	
###############################
def next_zero(lst, startIndex):	
	
	for i,value in enumerate(lst[startIndex:]):
		if value<0.01:
			return startIndex+i
	return -1

###############################
def get_height_duration(filename, dataset):

	f = h5py.File(filename, 'r')
	traces = np.array(f[dataset]) #baseline subtracted already
	time = np.array(f[dataset+"_t"])
	peakTime=np.array(f[dataset+"_peakTime"])
	
	toHist=[]
	durationArr=[]
	for trace in traces:
		height=np.max(trace)
		duration_start=np.argmax(trace)
		duration_end=next_zero(trace,duration_start)
		duration=time[duration_end]-time[duration_start]
		if duration==0:
			duration=1e-8
		#toHist.append(height/duration)
		durationArr.append(duration)

	#return toHist
	return durationArr

###############################	
def saveFromInput():
	inp= 'a'
	while inp!='y' and inp!='n':
		print("Save plot? y/n")
		inp = input()
	
	result=True if inp=='y' else False
	return result

###############################
def plotBinnedTraces(file, fileOut, ds, bins, strBins, aveTraces, addBaseline):

	f = h5py.File(file, 'r')
	time=get_time(file, ds)
	allTraces = f[ds] #baseline subtracted already	
	if addBaseline:
		baseline = f[ds+"_baseline"][:len(allTraces)]
	else:
		baseline = [0 for i in range(len(allTraces))]
	binIndices=get_binned_traces(bins,file,ds)
	colors=['r','b','g']
	
	#add baseline value back to traces
	#if addBaseline==False, then add 0 to trace	
	scaledTraces=[]
	for i,trace in enumerate(allTraces):
		scaledTraces.append(trace+baseline[i])
	
	binnedTraces=[[] for j in range(len(strBins))]
	for i,trace in enumerate(scaledTraces):
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
		if addBaseline:
			plt.ylabel("raw trace [V]")
		else:
			plt.ylabel("trace - baseline [V]")
		plot=plt.gcf() #get current figure - saves fig in case savePlt==True
		plt.show() #creates new figure for display
		savePlt=saveFromInput()
		if savePlt: #save plt stored with plt.gcf()
			if addBaseline:
				saveFile=saveDir+fileOut+"_bin"+strBins[j]+"_addBaseline.pdf"
			else:
				saveFile=saveDir+fileOut+"_bin"+strBins[j]+".pdf"
			print(f"Saving {saveFile}")
			plot.savefig(saveFile)	

###############################		
def plotBaseline(file, ds, fileOut):

	f = h5py.File(file, 'r')
	traces = f[ds] #baseline subtracted already
	baseline = f[ds+"_baseline"][:len(traces)]

	plt.hist(baseline, bins=50)
	plt.xlabel("Baseline [V] (average of points 0-2k)")
	plt.ylabel("Counts")
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	savePlt=saveFromInput()
	if savePlt: #save plt stored with plt.gcf()
		saveFile=saveDir+fileOut+"_baselineHist.pdf"
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)	


###############################
#Runnable choices
###############################

def plotTraces(files,labels,fileOut,ds=["run1"],bins=False, addBaseline=False):

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
	else:
		labels=[[x] for x in labels]
		#make label entries arrays of length 1 so full string is included in legend
			
	#for each file, calculate average trace (also calculate in bins, if bool is true)
	for f in files:
		file=homeDir+f
		j=0
		while j<len(ds):
			time=get_time(file,ds[j])
			trace=get_average_trace(file, ds[j], bins, addBaseline)
			i+=1
			for k,aveTrace in enumerate(trace):
				plt.plot(time*1e6, aveTrace, label=labels[i][k])
			j+=1	
	#plot average traces together
	plt.legend(loc="best")
	plt.xlabel( "time ($\mu$s)" )
	if addBaseline:
		plt.ylabel("raw trace [V]")
	else:
		plt.ylabel("trace - baseline [V]")
	plt.title(labels[0])
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	savePlt=saveFromInput()
	if savePlt: #save plt stored with plt.gcf()
		if addBaseline:
			saveFile=saveDir+fileOut+"_addBaseline.pdf"
		else:
			saveFile=saveDir+fileOut+".pdf"
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)	
	plt.clf()

	#plot individual traces along with average and/or baseline histograms
	for f in files:
		file=homeDir+f
		for dataset in ds:
			trace=get_average_trace(file, dataset, bins, addBaseline)
			if bins:
				plotBinnedTraces(file, fileOut, dataset, [0.09,0.2], strBinArr, trace, addBaseline) #bin array is high end of each bin - do not include infinity
			if addBaseline:
				plotBaseline(file,dataset,fileOut)	
		
#############################################################
def responseHistograms(filesIn, fileOut, ds="run1"):
	for f in filesIn:
		file=homeDir+f
		f = h5py.File(file, 'r')
		baselines=f[ds+"_baseline"]
		peaks=f[ds+"_peaks"]
		integrals=f[ds+"_integral"]
		peakTimes=np.array(f[ds+"_peakTime"])*1e6
		trigTimes=np.array(f[ds+"_trigTime"])
		deltaTrigs=[0]#placeholder for first element
		for i in range(len(trigTimes)-1):
			deltaTrigs.append(trigTimes[i+1]-trigTimes[i])
		durations=get_height_duration(file, ds)
		durations=[x*1e6 for x in durations]
		
		data=[peaks,integrals,peakTimes,durations,deltaTrigs]
		strdata=["peaks","integral","peakTimes", "durations","deltaTrigs"]
		xlabels=["Peak height [V]", "Integral [V]", "Peak Time (from trigger) [us]", "Pulse duration [us]","Time between triggers [s]"]
		
		for i,var in enumerate(data):
			plt.hist(var,bins=60)
			plt.xlabel(xlabels[i])
			plt.ylabel("counts")
			plot=plt.gcf() #get current figure - saves fig in case savePlt==True
			plt.show() #creates new figure for display
			savePlt=saveFromInput()
			if savePlt: #save plt stored with plt.gcf()
				saveFile=saveDir+fileOut+"_"+strdata[i]+"_full.pdf"
				print(f"Saving {saveFile}")
				plot.savefig(saveFile)	
			plt.clf()	
		
		#post processing to remove flat lines with baseline =-0.34 and have peak max within 50us after trigger (injections expect 10us)
		peaks_pp=[]
		integrals_pp=[]
		peakTimes_pp=[]
		durations_pp=[]
		deltaTrigs_pp=[]
		for i,bl in enumerate(baselines):
			if bl>0 and peakTimes[i]>0 and peakTimes[i]<50:
				peaks_pp.append(peaks[i])
				integrals_pp.append(integrals[i])
				peakTimes_pp.append(peakTimes[i])
				if i<len(durations):
					durations_pp.append(durations[i])
				deltaTrigs_pp.append(deltaTrigs[i])
				
		print(f"Before PP - {len(peaks)} events")
		print(f"After PP - {len(peaks_pp)} events (loose {(len(peaks)-len(peaks_pp))/len(peaks)*100}%)")

		data_pp=[peaks_pp,integrals_pp,peakTimes_pp,durations_pp,deltaTrigs_pp]

		for i,data in enumerate(data_pp):
			plt.hist(data, bins=60)
			plt.xlabel(xlabels[i])
			plt.ylabel("counts")
			plt.title("After Postprocessing - baseline>0 and 0<TrigTime<50us")
			plot=plt.gcf() #get current figure - saves fig in case savePlt==True
			plt.show() #creates new figure for display
			savePlt=saveFromInput()
			if savePlt: #save plt stored with plt.gcf()
				saveFile=saveDir+fileOut+"_"+strdata[i]+"_pp.pdf"
				print(f"Saving {saveFile}")
				plot.savefig(saveFile)	
			plt.clf()	
			
		#SCATTERPLOTS
		#peak time vs height
		plt.scatter(peakTimes,peaks, s=3)
		plt.xlabel("Peak time (from trigger) [us]")
		plt.ylabel("Pulse height [V]")
		plot=plt.gcf() #get current figure - saves fig in case savePlt==True
		plt.show() #creates new figure for display
		savePlt=saveFromInput()
		if savePlt: #save plt stored with plt.gcf()
			saveFile=saveDir+fileOut+"_peakTimeVsHeight_full.pdf"
			print(f"Saving {saveFile}")
			plot.savefig(saveFile)	
		plt.clf()
		
		#peak time vs duration
		plt.scatter(peakTimes[:100],durations,s=3)#only have 100 durations because they come from saved traces
		plt.xlabel("Peak time (from trigger) [us]")
		plt.ylabel("Pulse duration [us]")
		plot=plt.gcf() #get current figure - saves fig in case savePlt==True
		plt.show() #creates new figure for display
		savePlt=saveFromInput()
		if savePlt: #save plt stored with plt.gcf()
			saveFile=saveDir+fileOut+"_peakTimeVsDuration_full.pdf"
			print(f"Saving {saveFile}")
			plot.savefig(saveFile)	
		plt.clf()
		
		#peak height vs duration
		plt.scatter(peaks[:100],durations,s=3)#only have 100 durations because they come from saved traces
		plt.xlabel("Pulse height [V]")
		plt.ylabel("Pulse duration [us]")
		plot=plt.gcf() #get current figure - saves fig in case savePlt==True
		plt.show() #creates new figure for display
		savePlt=saveFromInput()
		if savePlt: #save plt stored with plt.gcf()
			saveFile=saveDir+fileOut+"_peakHeightVsDuration_full.pdf"
			print(f"Saving {saveFile}")
			plot.savefig(saveFile)	
		plt.clf()

#############################################################
def tracesInSeries(filesIn,fileOut,nmbTraces,ds="run1",firstTrace=0):

#plot traces as fn of time in series - see if can reconstruct peaks
#do not baseline correct - add back in
	for f in filesIn:
		file=homeDir+f
		f = h5py.File(file, 'r')
		traces=f[ds]
		baselines=f[ds+"_baseline"]
		traces=[x+y for x,y in zip(traces,baselines)]#add baselines back into traces
		peaks=f[ds+"_peaks"]
		integrals=f[ds+"_integral"]
		trigTimes=f[ds+"_trigTime"] #timestamps - in computer time, units=seconds
		timeAxis=f[ds+"_t"] #units=seconds

		#use all traces if nmbTraces==-1
		if nmbTraces<0:
			nmbTraces=len(traces)-1
		
		deltaTrigs=[0]#placeholder 0 for first element
		for i in range(len(trigTimes)-1):
			deltaTrigs.append(trigTimes[i+1]-trigTimes[i])
		

		for i in range(nmbTraces):
			plt.plot([x+sum(deltaTrigs[:i+1]) for x in timeAxis],traces[i])
			#advance plot by all time between all triggers that came before the trigger in question
		#plt.xscale('log')
		plt.xlabel('time [s]')
		plt.ylabel('recorded signal [V]')
		plot=plt.gcf() #get current figure - saves fig in case savePlt==True
		plt.show() #creates new figure for display
		savePlt=saveFromInput()
		if savePlt: #save plt stored with plt.gcf()
			saveFile=saveDir+fileOut+"_"+str(nmbTraces)+"tracesInSeries.pdf"
			print(f"Saving {saveFile}")
			plot.savefig(saveFile)	
		plt.clf()
		
		
		trigTimes=[datetime.fromtimestamp(x) for x in trigTimes]
		print(trigTimes[0])
		print(trigTimes[1])
	
###############################
#Main
###############################

if __name__ == "__main__":

	##Plot multiple average traces on one plot
	##Required arguments: input files, legend labels [first entry = plot title], output file name
	##Optional argument: array of datasets to compare (must be in same file: ex. run1 and and run2 from same input file). Default = only "run1"

	#filesIn=["040722_amp1/LBNL_inCave_withBeam_ion1_chip3_EBt_1.8Vinj_2min.h5py","040722_amp1/LBNL_inCave_withBeam_ion2_chip3_EBt_1.8Vinj_2min.h5py","040722_amp1/LBNL_inCave_withBeam_ion3_chip3_EBt_1.8Vinj_2min.h5py","040722_amp1/LBNL_inCave_withBeam_ion4_chip3_EBt_1.8Vinj_2min.h5py","040722_amp1/LBNL_inCave_withBeam_ion5_chip3_EBt_1.8Vinj_4min.h5py"]	
	#labels=["1.8V Injection in HI beam, Chip003, amp1", "ion1", "ion2", "ion3", "ion4", "ion5"]
	#fileOut="1.8Vinj_chip3_HIbeam_all"
	filesIn=["040722_amp1/LBNL_inCave_withBeam_ion2_chip3_EBt_1.8Vinj_2min.h5py"]	
	labels=["1.8V Injection in HI beam","Chip003, amp1"]
	fileOut="1.8Vinj_chip3_HIbeam_ion2"
	#filesIn=["040722_amp1/testLBNL_inCave_chip3_EBt_1.0Vinj_0.5min.h5py"]	
	#labels=["1.0V Injection in HI beam","Chip003, amp1"]
	#fileOut="1.0Vinj_chip3_HIbeam_inCave"
	plotTraces(filesIn,labels,fileOut,bins=True, addBaseline=True)
	#responseHistograms(filesIn,fileOut)
	#tracesInSeries(filesIn,fileOut,20) #not useful - triggers too far apart
	#tracesInSeries(filesIn,fileOut,2) 
	

#do any traces look like expected? (one single peak)		
#histogram of height duration on triggered pulse - full collection and also binned
#compare pulse duration to good injected run

#distribution of peak time vs height, peak time vs duration
		
		
		
		
		
		
		