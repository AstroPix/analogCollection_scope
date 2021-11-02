import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from scipy import stats
import sys,os
import enResFitting


savePlots=True



	
###########################################################################
# fit for peaks and compare location

enResArr=[]
muArr=[]
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/"

#trigThreshold=[550]
trigThreshold=[550,555,560,565,570,580,600,620,640,650]
dsList=[f"run1_{x}trig" for x in trigThreshold]
print(dsList)
file="102821_amp1/cadmium109_2min.h5py"

i=0
pixel=1
for ds in dsList:
	settings=[homeDir+file, f"Cd109 {trigThreshold[i]}mV trig", pixel, savePlots]
	popt, enRes, pcov = enResFitting.enResPlot(settings, dataset=ds)
	enResFitting.printParams(homeDir+file, popt, enRes, pcov, savePlots)
	poptI, enResI, pcovI = enResFitting.enResPlot(settings, integral=10000, dataset=ds)
	enResFitting.printParams(homeDir+file, poptI, enResI, pcovI, savePlots, integral=10000)
	enResArr.append([enRes, enResI])
	muArr.append([popt[1], poptI[1]])
	i+=1



# look at lowest recorded datapoint and linear fit to correlate threshold with measured energy
minPeakArr=[]
for ds in dsList:
	minPeakArr.append(enResFitting.getMin(homeDir+file,dataset=ds))
		
		
coef = np.polyfit(trigThreshold,minPeakArr,1)
poly1d_fn = np.poly1d(coef) 
plt.scatter(trigThreshold,minPeakArr,marker="o")
plt.plot(trigThreshold, poly1d_fn(trigThreshold), '--k',label=f"y={coef[0]:.3f}x{coef[1]:.3f}")
print(f"Full Linear fit: y={coef[0]:.3f}x{coef[1]:.3f}")
coef2 = np.polyfit(trigThreshold[:6],minPeakArr[:6],1)
poly1d_fn2 = np.poly1d(coef2) 
plt.plot(trigThreshold[:6], poly1d_fn2(trigThreshold[:6]), '--r',label=f"y={coef2[0]:.3f}x{coef2[1]:.3f}")
print(f"Partial Linear fit: y={coef2[0]:.3f}x{coef2[1]:.3f}")
plt.xlabel("Trigger threshold [mV]")
plt.ylabel("Minimum measured energy [V]")
plt.legend(loc="best")
plt.savefig(homeDir+file[:-5]+"_thresholdScan_fit.png") if savePlots else plt.show()



		
		
		
		
		
		
		
		
		