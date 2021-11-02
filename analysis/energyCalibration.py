import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from scipy import stats
import sys,os
import enResFitting


savePlots=False



def energyCalibFit(trueEn, data, dataName,saveto):
	
	dataNameStr=dataName.replace(" ", "")
	
	amp1_p=[]
	amp1_i=[]
	
	for photopeak in data:
		flatData=np.array(photopeak).flatten('F')
		amp1_p.append(flatData[0])
		amp1_i.append(flatData[1])
		
	coef = np.polyfit(amp1_p,trueEn,1)
	poly1d_fn = np.poly1d(coef) 
	plt.scatter(amp1_p,trueEn,label='Amp1',marker="o")
	plt.plot(amp1_p, poly1d_fn(amp1_p), '--k',label=f"y={coef[0]:.3f}x+{coef[1]:.3f}")
	print(f"Peak linear fit: y={coef[0]:.3f}x+{coef[1]:.3f}")
	plt.xlabel(f"{dataName} (from peak height)")
	plt.ylabel("True Energy [keV]")
	plt.legend(loc="best")
	#plt.yscale('log')
	#plt.show()
	#plt.savefig(f"{saveto}_peaks_{dataNameStr}_log.pdf")
	plt.savefig(f"{saveto}_peaks_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
	coef2 = np.polyfit(amp1_i,trueEn,1)
	poly1d_fn2 = np.poly1d(coef2) 	
	plt.scatter(amp1_i,trueEn,label='Amp1',marker="o")
	plt.plot(amp1_i, poly1d_fn2(amp1_i), '--k', label=f"y={coef2[0]:.3f}x+{coef2[1]:.3f}")
	print(f"Integral linear fit: y={coef2[0]:.3f}x+{coef2[1]:.3f}")
	plt.xlabel(f"{dataName} (from integral)")
	plt.ylabel("True Energy [keV]")
	plt.legend(loc="best")
	#plt.yscale('log')
	plt.legend(loc="best")
	#plt.show()
	#plt.savefig(f"{saveto}_integral_{dataNameStr}_log.pdf")
	plt.savefig(f"{saveto}_integral_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
	return coef,coef2
	
###########################################################################
#injection=[i*0.1 for i in range(1,19)]injection=[i*0.1 for i in range(1,3)]


enResArr=[]
muArr=[]
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/"

fileList=["102021_amp1/cadmium109_45min.h5py", "102821_amp1/cadmium109_16h.h5py", "102021_amp1/cobalt57_14h.h5py"]
energyList=[22.16, 88.03, 122.06, 39.46]
nameList=["Cadmium109", "Cadmium109", "Cobalt57"]

fitLow_p=[0.05,0.25,0.29]
fitLow_i=[150,725,1700]


"""
i=0
enRes_tmp, mu_tmp=[], []
for file in fileList:
	for pixel in [1]:
		settings=[homeDir+file, nameList[i], pixel, savePlots]
		popt, enRes, pcov = enResFitting.enResPlot(settings,fitLow=fitLow_p[i])
		enResFitting.printParams(file, popt, enRes, pcov, savePlots)
		poptI, enResI, pcovI = enResFitting.enResPlot(settings, fitLow=fitLow_i[i], integral=10000)
		enResFitting.printParams(file, poptI, enResI, pcovI, savePlots, integral=10000)
		enRes_tmp.append([enRes, enResI])
		mu_tmp.append([popt[1], poptI[1]])
	enResArr.append(enRes_tmp)
	muArr.append(mu_tmp)
	enRes_tmp,mu_tmp=[],[]
	i+=1
	
#Compton edge of Cobalt
file="102021_amp1/cobalt57_14h.h5py"
pixel=1
settings=[homeDir+file, "Cobalt57", pixel, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_edge(settings,fitLow=0.13, fitHigh=0.17)
popt, enRes, pcov = enResFitting.enResPlot_edge(settings,fitLow=500, fitHigh=1000,integral=10000)

"""







#trigScan fit: y=0.001x-0.6 (y=measured, x=trigThreshold)
#for Am, trig = 550 mV -> -0.05 measured
#zero point from baseline pulled from scope
#muArr=[[-0.05,-0.05],[0.0705,189.6605],[0.2597,742.7046],[0.3049,1812.8481]]
#energyList.insert(0,0)#zero energy point


muArr=[[0.0705,189.6605],[0.2597,742.7046],[0.3049,1812.8481],[0.151,641.77]]

coef_p,coef_i=energyCalibFit(energyList, muArr, "Fit Mean [V]",homeDir)	

file="102921_amp1/americium241_90min.h5py"
settings=[homeDir+file, "Americium241", 1, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_linearScale(settings,coef_p,fitLow=50)
enResFitting.printParams(file, popt, enRes, pcov, savePlots)
	

file="102021_amp1/cadmium109_45min.h5py"
settings=[homeDir+file, "Cadmium109", 1, savePlots]
popt, pcov = enResFitting.enResPlot_linearScale(settings,coef_p)
enResFitting.printParams(file, popt, enRes, pcov, savePlots)

		
		

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		