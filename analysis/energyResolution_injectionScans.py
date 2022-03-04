import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
import sys,os,glob
sys.path.insert(1, 'energyCalibration/')
import enResFitting


savePlots=True
#Fit peak heights or integral value (scaled by scope resolution)
#If traceInteg is true, use integral. If false, use peak height
traceInteg = True
v1=False #if v1=False, generate plots for v2

def injScanPlot(inj, data, err, dataName, saveto, chip, fit=False):
	
	dataNameStr=dataName.replace(" ", "")
	
	amp1_p=data[0]
	err1=err[0]
	amp2_p=data[1]
	err2=err[1]

	if v1:
		plt.title(f"Chip00{chip}")
		labels=['Amp1','Amp2']
	else:
		plt.title("v2")
		labels=['Chip1','Chip2']
	plt.errorbar(inj,amp1_p,yerr=err1,label=labels[0],marker="o",linestyle='None',color='r')
	plt.errorbar(inj,amp2_p,yerr=err2,label=labels[1],marker="o",linestyle='None')

	plt.xlabel("Injection voltage [V]")
	plt.legend(loc="best")
	if traceInteg:
		plt.ylabel(f"{dataName} (from integral)")
		plt.savefig(f"{saveto}integ_chip{chip}_{dataNameStr}.pdf") if savePlots else plt.show()
	else:
		plt.ylabel(f"{dataName} (from peak height)")
		plt.savefig(f"{saveto}peaks_chip{chip}_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()

	if fit:
		x = np.linspace(inj[1],inj[-3], 180)
		coef1i = np.polyfit(inj[1:-3],amp1_p[1:-3],1)
		poly1d_fn1i = np.poly1d(coef1i) 
		plt.plot(x, poly1d_fn1i(x), '--b',label=f"y={coef1i[0]:.3f}x+{coef1i[1]:.3f}")
		coef2i = np.polyfit(inj[1:-3],amp2_p[1:-3],1)
		poly1d_fn2i = np.poly1d(coef2i) 
		plt.plot(x, poly1d_fn2i(x), '--',color='orange',label=f"y={coef2i[0]:.3f}x+{coef2i[1]:.3f}")


	
###########################################################################
#Lowest injection values (0.1 - 0.3V) have scope scalings that are too large relative to the distribution
injection=[i*0.1 for i in range(5,14)] #1,19
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp"
if v1:
	saveDir = ["/Users/asteinhe/AstroPixData/astropixOut_tmp/injectionScans/v1/chip003/","/Users/asteinhe/AstroPixData/astropixOut_tmp/injectionScans/v1/chip003/"]
	runList = [[102221,1],[102221,2]] #day of run, pixel number #v1
else:
	saveDir = ["/Users/asteinhe/AstroPixData/astropixOut_tmp/injectionScans/v2/chip1/","/Users/asteinhe/AstroPixData/astropixOut_tmp/injectionScans/v2/chip2/"]
	runList=[['030122',1],['030322',1]]
binsize=0.002 
if traceInteg:
	binsize=0.01
elif v1:
	binsize=0 #default - based off scope resolution

enResArr, muArr, sigmaArr, muErrArr, sigErrArr = [],[],[],[],[]

#2 pixels of chip3 for v1
#compare 2 chips for v2
chip=3
for i,element in enumerate(runList):
	enResArr_tmp, muArr_tmp, sigmaArr_tmp, nArr_tmp, muErrArr_tmp, sigErrArr_tmp = [],[],[],[],[],[]
	for inj in injection:
		if v1:
			file=f"{homeDir}/v1/{element[0]}_amp{element[1]}/{inj:.1f}Vinj.h5py"
		elif i==0: #v2 chip1
			file=f"{homeDir}/v2/{element[0]}_amp{element[1]}/scan_{inj:.1f}Vinj_2min.h5py"
		else: #v2 chip2
			file=f"{homeDir}/v2/{element[0]}_amp{element[1]}/scan_chip2_{inj:.1f}Vinj_2min.h5py"
		settings=[file, f"{inj:.1f}V-injection", element[1], inj, savePlots, chip]
		popt, enRes, pcov = enResFitting.enResPlot(settings, savedir=saveDir[i], integral=traceInteg, injection=True, binSize=binsize)
		enResFitting.printParams(settings, popt, enRes, pcov, savedir=saveDir[i], integral=traceInteg, injection=True)
		enResArr_tmp.append(enRes)
		muArr_tmp.append(popt[1])
		sigmaArr_tmp.append(popt[2])
		muErrArr_tmp.append(np.sqrt(np.diag(pcov))[1])
		sigErrArr_tmp.append(np.sqrt(np.diag(pcov))[2])
	enResArr.append(enResArr_tmp)
	muArr.append(muArr_tmp)
	sigmaArr.append(sigmaArr_tmp)
	muErrArr.append(muErrArr_tmp)
	sigErrArr.append(sigErrArr_tmp)
#plot two pixels in each chip against each other				
injScanPlot(injection, enResArr, np.zeros(len(enResArr)), "Energy Resolution [\%]", saveDir[0], chip)#no easily returned value for enResErr
injScanPlot(injection, muArr, muErrArr, "Fit Mean [V]",saveDir[0],chip,fit=True)		
		
if v1:	
	#compare 2 pixels of chip4
	saveDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/injectionScans/v1/chip004/"
	enResArr, muArr, sigmaArr, muErrArr, sigErrArr = [],[],[],[],[]
	runList = [[120721,1],[120721,2]] #day of run, pixel number
	#changed the name of the data files
	chip=4
	for i,element in enumerate(runList):
		enResArr_tmp, muArr_tmp, sigmaArr_tmp, nArr_tmp, muErrArr_tmp, sigErrArr_tmp = [],[],[],[],[],[]
		for inj in injection:
			file=f"{homeDir}/{element[0]}_amp{element[1]}/chip004_AC_{inj:.1f}Vinj_2min.h5py"
			settings=[file, f"{inj:.1f}V injection", element[1], inj, savePlots, chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings, savedir=saveDir, integral=traceInteg, injection=True)
			enResFitting.printParams(settings, popt, enRes, pcov, savedir=saveDir, integral=traceInteg, injection=True)
			enResArr_tmp.append(enRes)
			muArr_tmp.append(popt[1])
			sigmaArr_tmp.append(popt[2])
			muErrArr_tmp.append(np.sqrt(np.diag(pcov))[1])
			sigErrArr_tmp.append(np.sqrt(np.diag(pcov))[2])
		enResArr.append(enResArr_tmp)
		muArr.append(muArr_tmp)
		sigmaArr.append(sigmaArr_tmp)
		muErrArr.append(muErrArr_tmp)
		sigErrArr.append(sigErrArr_tmp)
	#plot two pixels in each chip against each other				
	injScanPlot(injection, enResArr, np.zeros(len(enResArr)), "Energy Resolution [\%]", saveDir, chip)#no easily returned value for enResErr
	injScanPlot(injection, muArr, muErrArr, "Fit Mean [V]",saveDir,chip,fit=True)			
	
		