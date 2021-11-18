import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from scipy import stats, interpolate
from scipy.odr import ODR, Model, Data, RealData
import sys,os,glob
import enResFitting


savePlots=True
fitSpectra=False #False if spectra have already been fit
pix=1 #1 or 2 - which amp to consider

#Fit functions for curve_fit
def linFit(x,m,b):
	return m*x + b
def quadFit(x,a,b,c):
	return a*x*x + b*x + c
def triFit(x,a,b,c,d):
	return a*x*x*x + b*x*x + c*x + d	
def sqrtFit(x,A,B,mu):
	#return A*np.sqrt(x-mu)+B
	return A*np.sqrt(x)+B
	
	
#fit functions for ODR	
def lin_odr(p,x):
	m,b=p
	return m*x + b
def quad_odr(p, x):
	a, b, c = p
	return a * x *x + b*x + c
def tri_odr(p,x):
	a,b,c,d = p
	return a*x*x*x + b*x*x + c*x + d	
def sqrt_odr(p,x):
	A,B=p
	return A*np.sqrt(x)+B
	
def odr_polyfit(fitdata, deg):
	if deg==1:
		mod =  Model(lin_odr)
	elif deg==2:
		mod = Model(quad_odr)
	elif deg==3:
		mod = Model(tri_odr)
	else:
		printf("Not possible - choose degree 1, 2, or 3")
		#AMANDA - better break
	
	
	odrfit = ODR(fitdata,mod,[1e-8 for x in range(deg+1)])		
	out=odrfit.run()
	coef=out.beta
	
	return coef

def getFiles(amp):
	if amp==1:
		fileList=["102021_amp1/cobalt57_14h.h5py","102021_amp1/cadmium109_45min.h5py", "102821_amp1/cadmium109_16h.h5py", "102021_amp1/cobalt57_14h.h5py","110421_amp1/Americium_480min_combined.h5py","110821_amp1/barium133_combined_65min.h5py"]
		energyList=[14.41, 22.16, 88.03, 122.06, 59.54, 30.97]
		nameList=["Cobalt57", "Cadmium109", "Cadmium109", "Cobalt57", "Americium241", "Barium133"]
		#fine-tune the range so that individual peaks are picked out
		fitLow_p=[0.02, 0.05,0.25,0.29,0.17,0.09]
		fitLow_i=[0,150,725,1700,300,200,100]
		return fileList,energyList,nameList, fitLow_p, fitLow_i
	elif amp==2:
		fileList=["111221_amp2/weekend_Cobalt57_4020min.h5py","111221_amp2/weekend_Cobalt57_4020min.h5py","111521_amp2/overnight_Americium241_960min.h5py", "111621_amp2/day_Cadmium109_300min.h5py"]
		energyList=[14.41, 122.06, 59.54, 22.16]
		nameList=["Cobalt57","Cobalt57", "Americium241", "Cadmium109"]
		#fine-tune the range so that individual peaks are picked out
		fitLow_p=[0.03,0.29,0.09,0.05]
		fitLow_i=[0,1700,300,100,150]
		return fileList,energyList,nameList, fitLow_p, fitLow_i
	else:
		print("Choose amp1 or amp2")
		return 0
		
def getEdgeFiles(amp):
	if amp==1:
		fileList=["102021_amp1/cobalt57_14h.h5py"]
		energyList=[39.46]
		nameList=["Cobalt57"]
		fitLow_p=[0.13]
		fitHigh_p=[0.17]
		fitLow_i=[0]
		fitHigh_i=[1000]
		return fileList,energyList,nameList, fitLow_p, fitLow_i, fitHigh_p, fitHigh_i
	elif amp==2:
		print("No files")
		return 0
	else:
		print("Choose amp1 or amp2")
		return 0
	
def energyCalibFit(trueEn, data, err, dataName, saveto):
	plt.clf()
	
	#get name for plotting
	dataNameStr=dataName.replace(" ", "")
	
	#fill arrays with data passed in
	amp_p=[]
	amp_i=[]
	for photopeak in data:
		flatData=np.array(photopeak).flatten('F')
		amp_p.append(flatData[0])
		amp_i.append(flatData[1])
	err_p=[]
	err_i=[]
	for noise in err:
		flatErr=np.array(noise).flatten('F')
		err_p.append(flatErr[0])
		err_i.append(flatErr[1])
		
	#organize arrays for interpolation
	trueEn_sorted= sorted(trueEn)
	amp_p_sorted = [x for _, x in sorted(zip(trueEn, amp_p), key=lambda pair: pair[0])]
	err_p_sorted = [x for _, x in sorted(zip(trueEn, err_p), key=lambda pair: pair[0])]
	amp_i_sorted = [x for _, x in sorted(zip(trueEn, amp_i), key=lambda pair: pair[0])]
	err_i_sorted = [x for _, x in sorted(zip(trueEn, err_i), key=lambda pair: pair[0])]


	#Fit different functions to data
	coef, coef_pcov = curve_fit(linFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #linear 
	coef2, coef2_pcov = curve_fit(quadFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #quadratic
	coef3, coef3_pcov = curve_fit(triFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #3rd deg poly
	popt, pcov = curve_fit(sqrtFit, trueEn, amp_p, sigma=err_p,absolute_sigma=True) #square root
	tck = interpolate.splrep(trueEn_sorted, amp_p_sorted) #cubic spline
	tck1 = interpolate.splrep(trueEn_sorted, amp_p_sorted,k=1) #linear spline

	
	#Plot data and fit functions
	x = np.linspace(np.min(trueEn), np.max(trueEn), 100)
	plt.errorbar(trueEn, amp_p, yerr=err_p, fmt='o', label="data")
	linPlot=linFit(x,*coef)
	#plt.plot(x, linPlot, '--k',label=f"Linear fit")
	quadPlot=quadFit(x,*coef2)
	#plt.plot(x, quadPlot, '--r',label=f"y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
	triPlot=triFit(x,*coef3)
	#plt.plot(x, triPlot, '--b',label=f"3rd deg. polynomial")
	sqrt_fn = sqrtFit(x, *popt)
	#plt.plot(x,sqrt_fn,'--g',label=f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
	splinePlot = interpolate.splev(x, tck)
	#plt.plot(x, splinePlot, '--g', label="Cubic Spline")
	spline1Plot = interpolate.splev(x, tck)
	plt.plot(x, spline1Plot, '--g', label="linear spline")

	plt.xlabel("True Energy [keV]")
	plt.ylabel(f"{dataName} (from peak height)")	
	#plt.ylabel(f"{dataName} (from integral)")
	plt.legend(loc="best")
	plt.grid()
	plt.xlim([-10,1.1*np.max(trueEn)])
	plt.ylim([-0.05,1.1*np.max(amp_p)])
	#plt.yscale('log')
	plt.savefig(f"{saveto}peaks_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
	#Print fit equations to terminal
	print(f"Peak linear fit: y={coef[0]:.3f}x+{coef[1]:.3f}")
	print(f"Peak quadratic fit: y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
	print(f"Peak 3rd deg poly fit: y={coef3[0]:.3f}x$^3$+{coef3[1]:.3f}x$^2$+{coef3[2]:.3f}x + {coef3[3]:.3f}")
	print(f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")

	#fit opposite orientation so that scaling is easier
	#need different regression technique for X-error bars
	datain = RealData(amp_p, trueEn, sx=err_p)
	coef_fit = odr_polyfit(datain,1)
	coef2_fit = odr_polyfit(datain,2)
	coef3_fit = odr_polyfit(datain,3)
	sqrt_model = Model(sqrt_odr)
	odr_sqrt = ODR(datain, sqrt_model,[1e-8,1e-8])
	out_sqrt=odr_sqrt.run()
	coef4_fit=out_sqrt.beta
	#spline
	splineFn = interpolate.splrep(amp_p_sorted, trueEn_sorted, w=err_p_sorted)
	spline1Fn = interpolate.splrep(amp_p_sorted, trueEn_sorted, w=err_p_sorted, k=1)

	
	#print(f"Peak 3rd deg poly fit, x [keV]: y={coef3[0]}x$^3$+{coef3[1]}x$^2$+{coef3[2]}x + {coef3[3]}")
	#print(f"Peak 3rd deg poly fit, x [V]: y={coef3_fit[0]}x$^3$+{coef3_fit[1]}x$^2$+{coef3_fit[2]}x + {coef3_fit[3]}")

		
	
	#AMANDA - goodness of fit value
	#return ideal fit only - 3rd deg polynomial
	
	#return coef3_fit
	return spline1Fn
	
	
	
########################################################################################
###########################################################################
##############################################################

homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/"
saveDir = enResFitting.getSaveto()
dataDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp1_peaks/fitSpectra"


#files to be used for energy calibration curve
#amp1
enResArr1, muArr1, sigmaArr1, nArr1 = [],[],[], []

#amp2
enResArr2, muArr2, sigmaArr2, nArr2 = [],[],[], []
#fileList,energyList,nameList = getFiles(2)


if fitSpectra:
	#amp1
	fileList,energyList,nameList,fitLow_p,fitLow_i = getFiles(pix)

	
	#loop through all files, fit with Gaussian, return mean/sigma/energy resolution and store in arrays - separately for each amp
	i=0
	for file in fileList:
		settings=[homeDir+file, nameList[i], pix, energyList[i], savePlots]
		popt, enRes, pcov, integ = enResFitting.enResPlot(settings,fitLow=fitLow_p[i])
		enResFitting.printParams(settings, integ, popt, enRes, pcov)
		#integral argument = integral bin size (in V)
		poptI, enResI, pcovI, integI = enResFitting.enResPlot(settings, fitLow=fitLow_i[i], integral=20)
		enResFitting.printParams(settings, integI, poptI, enResI, pcovI, integral=True)
		enResArr1.append([enRes,enResI])
		muArr1.append([popt[1],poptI[1]])
		sigmaArr1.append([popt[2],poptI[2]])
		nArr1.append([integ,integI])
		i+=1
	
	#loop through additional files, fit with integrated Gaussian (for Compton edge), return mean/sigma/energy resolution and store in arrays
	edgeFileList,edgeEnergyList,edgeNameList,edgeFitLow_p,edgeFitLow_i,edgeFitHigh_p,edgeFitHigh_i = getEdgeFiles(pix)
	energyList.extend(edgeEnergyList)
	
	i=0
	for file in edgeFileList:
		settings=[homeDir+file, edgeNameList[i], pix, edgeEnergyList[i],savePlots]
		popt, enRes, pcov, integ= enResFitting.enResPlot_edge(settings,fitLow=edgeFitLow_p[i], fitHigh=edgeFitHigh_p[i])
		enResFitting.printParams_edge(settings, integ, popt, enRes, pcov)
		poptI, enResI, pcovI, integI = enResFitting.enResPlot_edge(settings,fitLow=edgeFitLow_i[i], fitHigh=edgeFitHigh_i[i],integral=10000)
		enResFitting.printParams_edge(settings, integI, poptI, enResI, pcovI, integral=True)
		enResArr1.append([enRes,enResI])
		muArr1.append([popt[2],poptI[2]])
		sigmaArr1.append([popt[3],poptI[3]])
		nArr1.append([integ,integI])
		i+=1
		
else: #if spectra have been fit before, pull out values from txt files
	print("spectra are already fit")
	energyList=[]

	os.chdir(dataDir)
	peakFiles = glob.glob('*peaks*.txt')
	for filename in peakFiles:
		energyList.append(float(filename.split('_')[1][:-4]))
		openFile=open(filename,'r')
		lines=openFile.readlines()
		muArr1.append([float(lines[1].split(' = ')[-1])])
		sigmaArr1.append([float(lines[2].split(' = ')[-1])])
		nArr1.append([float(lines[3].split(' = ')[-1])])
		enResArr1.append([float(lines[4].split(' = ')[-1][:-2])])#eliminate % sign at the end
		
	intFiles = glob.glob('*integral*.txt')
	for filename in intFiles:
		energy=float(filename.split('_')[1][:-4])
		energyIndex=energyList.index(energy)
		openFile=open(filename,'r')
		lines=openFile.readlines()
		muArr1[energyIndex].append(float(lines[1].split(' = ')[-1]))
		sigmaArr1[energyIndex].append(float(lines[2].split(' = ')[-1]))
		nArr1[energyIndex].append(float(lines[3].split(' = ')[-1]))
		enResArr1[energyIndex].append(float(lines[4].split(' = ')[-1][:-2]))#eliminate % sign at the end



"""
#to debug fitting - hardcode values
energyList.append(39.46)
muArr1=[[0.0705,190.149],[0.2597,742.7046],[0.3049,1812.8481],[0.2042,463.842],[0.116,153.671],[0.151,641.77]]
sigmaArr1=[[0.012,53.393],[0.0054,28.206],[0.0076,56.683],[0.0061,41.639],[0.012,51.320],[0.0033,65.1507]]
nArr1=[[6.335,166564.19],[0.678,2.164],[0.809,2.794],[0.428,21722.894],[1.312, 14189.19],[2.234,67744.383]]
#AMANDA - N calculation is failing for integral
"""


#calculate error
errArr1=[]

#error Arr1ay = sigma/sqrt(N) (for edge, 2sig integral from mu)
for j in range(len(sigmaArr1)):
	err_p=sigmaArr1[j][0]/np.sqrt(nArr1[j][0])
	err_i=sigmaArr1[j][1]/np.sqrt(nArr1[j][1])
	errArr1.append([err_p,err_i])


#use fit mean of measured peaks and associated error to create calibration curve
coef_p = energyCalibFit(energyList, muArr1, errArr1, "Fit Mean [V]",saveDir)


#use calibration curve to calibrate a spectrum
file="110421_amp1/Americium_480min_combined.h5py"
settings=[homeDir+file, "Americium241-calib", 1, 59.54, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_scale(settings,coef_p,fitLow=50)
enResFitting.printParams(settings, -1, popt, enRes, pcov, savePlots)
	

file="102021_amp1/cadmium109_45min.h5py"
settings=[homeDir+file,  "Cadmium109-calib", 1, 22.16, savePlots]
popt, enRes, pcov = enResFitting.enResPlot_scale(settings,coef_p)
enResFitting.printParams(settings, -1, popt, enRes, pcov, savePlots)



		
#AMANDA
#18 NOV
#add 14 keV Co point before creating all fit plots
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		