import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from scipy import stats
from scipy.odr import ODR, Model, Data, RealData
import sys,os
import enResFitting


savePlots=False

def sqrtFit(x,A,B,mu):
	#return A*np.sqrt(x-mu)+B
	return A*np.sqrt(x)+B
	
def linFit(x,m,b):
	return m*x + b
	
def quadFit(x,a,b,c):
	return a*x*x + b*x + c
	
def lin_odr(p,x):
	m,b=p
	return m*x + b
def quad_odr(p, x):
	a, b, c = p
	return a * x *x + b*x + c
def sqrt_odr(p,x):
	A,B=p
	return A*np.sqrt(x)+B
	
def energyCalibFit(trueEn, data, err, dataName, saveto):
	plt.clf()
	
	#get name for plotting
	dataNameStr=dataName.replace(" ", "")
	
	#fill arrays with data passed in
	amp1_p=[]
	amp1_i=[]
	for photopeak in data:
		flatData=np.array(photopeak).flatten('F')
		amp1_p.append(flatData[0])
		amp1_i.append(flatData[1])
	err_p=[]
	err_i=[]
	for noise in err:
		flatErr=np.array(noise).flatten('F')
		err_p.append(flatErr[0])
		err_i.append(flatErr[1])
		

	#Fit different functions to data
	coef, coef_pcov = curve_fit(linFit,trueEn,amp1_p,sigma=err_p,absolute_sigma=True) #linear 
	coef2, coef2_pcov = curve_fit(quadFit,trueEn,amp1_p,sigma=err_p,absolute_sigma=True) #quadratic
	popt, pcov = curve_fit(sqrtFit, trueEn, amp1_p, sigma=err_p,absolute_sigma=True) #square root
	
	
	#Plot data and fit functions
	x = np.linspace(np.min(trueEn), np.max(trueEn), 100)
	plt.errorbar(trueEn, amp1_p, yerr=err_p, fmt='o', label="data")
	linPlot=linFit(x,*coef)
	plt.plot(x, linPlot, '--k',label=f"y={coef[0]:.3f}x+{coef[1]:.3f}")
	quadPlot=quadFit(x,*coef2)
	plt.plot(x, quadPlot, '--r',label=f"y={coef2[0]:.4f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
	sqrt_fn = sqrtFit(x, *popt)
	plt.plot(x,sqrt_fn,'--g',label=f"y={popt[0]:.3f}*sqrt(x)$+{popt[1]:.3f}")
	plt.xlabel("True Energy [keV]")
	plt.ylabel(f"{dataName} (from peak height)")
	plt.legend(loc="best")
	plt.grid()
	plt.xlim([-10,1.1*np.max(trueEn)])
	plt.ylim([-0.05,1.1*np.max(amp1_p)])
	#plt.yscale('log')
	plt.savefig(f"{saveto}_peaks_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
	#Print fit equations to terminal
	print(f"Peak linear fit: y={coef[0]:.3f}x+{coef[1]:.3f}")
	print(f"Peak quadratic fit: y={coef2[0]:.4f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
	print(f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")

	#fit opposite orientation so that scaling is easier
	#need different regression technique for X-error bars
	#coef_fit, coef_pcov_fit = curve_fit(linFit,amp1_p,trueEn) #linear
	#coef2_fit, coef2_pcov_fit = curve_fit(quadFit,amp1_p,trueEn) #quadratic
	#popt_fit,pcov_fit = curve_fit(sqrtFit, amp1_p,trueEn)#sqrt
	
	
	datain = RealData(amp1_p, trueEn, sx=err_p)
	lin_model = Model(lin_odr)
	quad_model = Model(quad_odr)
	sqrt_model = Model(sqrt_odr)
	odr_lin = ODR(datain,lin_model,[1,1e-8])
	odr_quad = ODR(datain, quad_model,[1e-8,1e-8,1e-8])
	odr_sqrt = ODR(datain, sqrt_model,[1e-8,1e-8])
	out_lin=odr_lin.run()
	out_quad = odr_quad.run()
	out_sqrt=odr_sqrt.run()
	coef_fit = out_lin.beta
	coef2_fit=out_quad.beta
	coef3_fit=out_quad.beta
	

	
	#AMANDA - goodness of fit value
	#return ideal fit only - quad
	
	return coef2_fit
	
###########################################################################
#injection=[i*0.1 for i in range(1,19)]injection=[i*0.1 for i in range(1,3)]


enResArr=[]
muArr=[]
homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/"

fileList=["102021_amp1/cadmium109_45min.h5py", "102821_amp1/cadmium109_16h.h5py", "102021_amp1/cobalt57_14h.h5py","102921_amp1/americium241_90min.h5py"]
energyList=[22.16, 88.03, 122.06, 59.54]
nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Americium241"]

fitLow_p=[0.05,0.25,0.29,0.17]
fitLow_i=[150,725,1700,300,200]


"""
i=0
enRes_tmp, mu_tmp=[], []
for file in fileList:
	for pixel in [1]:
		settings=[homeDir+file, nameList[i], pixel, savePlots]
		popt, enRes, pcov = enResFitting.enResPlot(settings,fitLow=fitLow_p[i])
		enResFitting.printParams(homeDir+file, popt, enRes, pcov, savePlots)
		poptI, enResI, pcovI = enResFitting.enResPlot(settings, fitLow=fitLow_i[i], integral=10000)
		enResFitting.printParams(homeDir+file, poptI, enResI, pcovI, savePlots, integral=10000)
		enRes_tmp.append([enRes, enResI])
		mu_tmp.append([popt[1], poptI[1]])
	enResArr.append(enRes_tmp)
	muArr.append(mu_tmp)
	enRes_tmp,mu_tmp=[],[]
	i+=1
	


#Compton edge of Cobalt
energyList.append(39.46)
file="102021_amp1/cobalt57_14h.h5py"
pixel=1
settings=[homeDir+file, "Cobalt57", pixel, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_edge(settings,fitLow=0.13, fitHigh=0.17)
enResFitting.printParams_edge(homeDir+file, popt, enRes, pcov, savePlots)
poptI, enResI, pcovI = enResFitting.enResPlot_edge(settings,fitLow=500, fitHigh=1000,integral=10000)
enResFitting.printParams_edge(homeDir+file, poptI, enResI, pcovI, savePlots, integral=10000)
enResArr.append([enRes,enResI])
muArr.append([popt[2],poptI[2]])

"""







energyList.append(39.46)

muArr=[[0.0705,189.6605],[0.2597,742.7046],[0.3049,1812.8481],[0.2042,471.45],[0.151,641.77]]
muArrErr=[[0.000105,0.518],[0.0006, 1.657],[0.0004,2.687],[0.0007,5.805],[0.0007,8.437]]
sigmaArr=[[0.012,53.817],[0.0054,28.206],[0.0076,56.683],[0.0061,33.808],[0.0033,65.1507]]
nArr=[[6.335,67289],[0.678,2.164],[0.809,2.794],[0.428,1.417],[2.234,67744.383]]
errArr=[]

#error array = sigma/sqrt(N) (for edge, 2sig integral from mu)
for j in range(len(sigmaArr)):
	err_p=sigmaArr[j][0]/np.sqrt(nArr[j][0])
	err_i=sigmaArr[j][1]/np.sqrt(nArr[j][1])
	errArr.append([err_p,err_i])


#zero energy point
muArr.insert(0,[1e-8,1e-8])
energyList.insert(0,1e-8)
errArr.insert(0,[1e-8,1e-8])

coef_p =energyCalibFit(energyList, muArr, errArr, "Fit Mean [V]",homeDir)	

file="102921_amp1/americium241_90min.h5py"
settings=[homeDir+file, "Americium241", 1, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_scale(settings,coef_p,fitLow=50)
enResFitting.printParams(file, popt, enRes, pcov, savePlots)
	

file="102021_amp1/cadmium109_45min.h5py"
settings=[homeDir+file, "Cadmium109", 1, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_scale(settings,coef_p)
enResFitting.printParams(file, popt, enRes, pcov, savePlots)



		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		