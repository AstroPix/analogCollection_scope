import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
import sys,os
import energyCalibration/enResFitting


savePlots=True



def injScanPlot(inj, data, dataName, saveto, fit=False):
	
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

	plt.scatter(inj,amp1_p,label='Amp1',marker="o",linestyle='None')
	plt.scatter(inj,amp2_p,label='Amp2',marker="o",linestyle='None')
	plt.xlabel("Injection voltage [V]")
	plt.ylabel(f"{dataName} (from peak height)")
	plt.legend(loc="best")
	plt.savefig(f"{saveto}_peaks_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()

	if fit:
		x = np.linspace(inj[1],inj[-3], 180)
		coef1i = np.polyfit(inj[1:-3],amp1_i[1:-3],1)
		poly1d_fn1i = np.poly1d(coef1i) 
		coef2i = np.polyfit(inj[1:-3],amp2_i[1:-3],1)
		poly1d_fn2i = np.poly1d(coef2i) 
		plt.plot(x, poly1d_fn1i(x), '--b',label=f"y={coef1i[0]:.3f}x+{coef1i[1]:.3f}")
		plt.plot(x, poly1d_fn2i(x), '--',color='orange',label=f"y={coef2i[0]:.3f}x+{coef2i[1]:.3f}")

	plt.plot(inj,amp1_i,label='Amp1',marker="o",linestyle='None')
	plt.plot(inj,amp2_i,label='Amp2',marker="o",linestyle='None')
	plt.xlabel("Injection voltage [V]")
	plt.ylabel(f"{dataName} (from integral)")
	plt.legend(loc="best")
	plt.savefig(f"{saveto}_integral_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
	
###########################################################################
injection=[i*0.1 for i in range(1,19)]
enResArr=[]
muArr=[]
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp"
enRes_tmp, mu_tmp=[], []
for inj in injection:
	for pixel in [1, 2]:
		file=f"{homeDir}/102221_amp{pixel}/{inj:.1f}Vinj.h5py"
		settings=[file, f"{inj:.1f}V injection", pixel, savePlots]
		popt, enRes, pcov, integ = enResFitting.enResPlot(settings)
		enResFitting.printParams(file, popt, enRes, pcov, savePlots)
		#integral argument = integral bin size (in V)
		poptI, enResI, pcovI, integ = enResFitting.enResPlot(settings, integral=5*(inj+0.1))
		enResFitting.printParams(file, poptI, enResI, pcovI, savePlots, integral=True)
		enRes_tmp.append([enRes, enResI])
		mu_tmp.append([popt[1], poptI[1]])
	enResArr.append(enRes_tmp)
	muArr.append(mu_tmp)
	enRes_tmp,mu_tmp=[],[]
		
		
		
injScanPlot(injection, enResArr, "Energy Resolution [\%]", homeDir+"/102221")	
injScanPlot(injection, muArr, "Fit Mean [V]",homeDir+"/102221", fit=True)		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		