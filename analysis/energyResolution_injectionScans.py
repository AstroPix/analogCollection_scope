import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
import sys,os
import enResFitting


savePlots=False



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

	coef1p = np.polyfit(trueEn,amp1_p,2)
	poly2d_fn1p = np.poly2d(coef1p) 
	coef2p = np.polyfit(trueEn,amp2_p,2)
	poly2d_fn2p = np.poly2d(coef2p) 
	plt.scatter(inj,amp1_p,label='Amp1',marker="o")
	plt.scatter(inj,amp2_p,label='Amp2',marker="o")
	plt.plot(trueEn, poly2d_fn1p(trueEn), '--k',label=f"y={coef1p[0]:.3f}x$^2$+{coef1p[1]:.3f}x+{coef1p[2]:.3f}")
	plt.plot(trueEn, poly2d_fn2p(trueEn), '--k',label=f"y={coef2p[0]:.3f}x$^2$+{coef2p[1]:.3f}x+{coef2p[2]:.3f}")
	plt.xlabel("Injection voltage [V]")
	plt.ylabel(f"{dataName} (from peak height)")
	plt.legend(loc="best")
	plt.yscale('log')
	#plt.show()
	#plt.savefig(f"{saveto}_peaks_{dataNameStr}_log.pdf")
	plt.savefig(f"{saveto}_peaks_{dataNameStr}_log.pdf") if savePlots else plt.show()
	plt.clf()

	coef1i = np.polyfit(trueEn,amp1_p,2)
	poly2d_fn1i = np.poly2d(coef1i) 
	coef2i = np.polyfit(trueEn,amp2_p,2)
	poly2d_fn2i = np.poly2d(coef2i) 	
	plt.plot(inj,amp1_i,label='Amp1',marker="o")
	plt.plot(inj,amp2_i,label='Amp2',marker="o")
	plt.plot(trueEn, poly2d_fn1i(trueEn), '--k',label=f"y={coef1i[0]:.3f}x$^2$+{coef1i[1]:.3f}x+{coef1i[2]:.3f}")
	plt.plot(trueEn, poly2d_fn2i(trueEn), '--k',label=f"y={coef2i[0]:.3f}x$^2$+{coef2i[1]:.3f}x+{coef2i[2]:.3f}")
	plt.xlabel("Injection voltage [V]")
	plt.ylabel(f"{dataName} (from integral)")
	plt.legend(loc="best")
	plt.yscale('log')
	#plt.show()
	#plt.savefig(f"{saveto}_integral_{dataNameStr}_log.pdf")
	plt.savefig(f"{saveto}_integral_{dataNameStr}_log.pdf") if savePlots else plt.show()
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
		popt, enRes, pcov = enResFitting.enResPlot(settings)
		enResFitting.printParams(file, popt, enRes, pcov, savePlots)
		poptI, enResI, pcovI = enResFitting.enResPlot(settings, integral=1000)
		enResFitting.printParams(file, poptI, enResI, pcovI, savePlots, integral=1000)
		enRes_tmp.append([enRes, enResI])
		mu_tmp.append([popt[1], poptI[1]])
	enResArr.append(enRes_tmp)
	muArr.append(mu_tmp)
	enRes_tmp,mu_tmp=[],[]
		
		
		
injScanPlot(injection, enResArr, "Energy Resolution [\%]", homeDir+"/102221")	
injScanPlot(injection, muArr, "Fit Mean [V]",homeDir+"/102221")		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		