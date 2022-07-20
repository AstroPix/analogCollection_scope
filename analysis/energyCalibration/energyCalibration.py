import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import lstsq
import h5py
import scipy
from scipy import stats, interpolate, optimize
from scipy.optimize import curve_fit
from scipy.odr import ODR, Model, Data, RealData
from scipy.interpolate import UnivariateSpline
import sys,os,glob
import enResFitting


#Read in settings from runOptions.txt
inputDict={}
with open(sys.argv[1], 'r') as runOptions:
	lines=runOptions.readlines()
	for line in lines:
		if "=" in line:
			readLine = line.split(' = ')
			inputDict[readLine[0]] = readLine[1][:-1]
#convert dictionary values into variables
for key,val in inputDict.items():
        exec(key + '=val')
#typeset if necessary
savePlots=bool(int(savePlots))
fitSpectra=bool(int(fitSpectra))
traceInteg=bool(int(traceInteg))
vers=int(version)
pix=int(pix)
fit=int(fit)
chip=int(chip)

#AMANDA UPDATE
#if using v2, use same noise cut as v1 chip3 pix 1 (20 mV)
#eliminate v2 and amp2 option

########################################################################################
###########################################################################
##############################################################

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
def piecewise_linear(x, x0, y0, k1, k2):
	return np.piecewise(x, [x < x0, x>=x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
	
	
def getSumSq(trueEn, voltage, err, fn):
	sum_square=0
	for i in range(len(trueEn)):
		fitVal=fn(voltage[i])
		num = (trueEn[i]-fitVal)**2
		fitErr=fn(err[i])
		denom = fitErr**2
		sum_square += (num/denom)
		i+=1
	return sum_square

def fixFileList(fileList):
	#complete path if not already done in input file, eliminate empty elements
	todel=[]
	for i, inp in enumerate(fileList):
		if inp=='\n':
			todel.append(i)
		elif homeDir not in inp:
			fileList[i]=homeDir+inp	
		fileList[i]=fileList[i][:-1]#remove new line characters
	for i in sorted(todel, reverse=True):
		del fileList[i]	
	return fileList

def getFiles(amp):
	os.chdir(sys.path[0])
	if amp==1 and vers==1:
		with open(input1, 'r') as files1:
			fileList=files1.readlines()
		if chip==3:
			energyList=[22.16, 88.03, 14.41, 122.06, 59.54, 30.97]
			binsize=np.zeros(len(energyList))
			nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Cobalt57", "Americium241", "Barium133"]
			#fine-tune the range so that individual peaks are picked out
			if not traceInteg:
				fitLowArr=[0.05,0.244,0.02,0.3,0.19,0.09]
				fitHighArr=[0.11,1,0.07,0.32,1,0.15]
			else:
				fitLowArr=[0,70,0,70,40,10]
				fitHighArr=[30,80,100,80,60,100]
		elif chip==4:
			energyList=[22.16, 88.03, 14.41, 122.06]
			binsize=np.zeros(len(energyList))
			nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Cobalt57"]
			if not traceInteg:
				fitLowArr=[0.0,0.244,0.0,0.29]
				fitHighArr=[0.11,1,0.04,0.3]
			else:
				fitLowArr=[0,50,0,80]
				fitHighArr=[30,80,100,120]
	elif amp==1 and vers==2:
		with open(input1, 'r') as files1:
			fileList=files1.readlines()
		if chip==1:
			energyList=[22.16, 88.03,122.06,14.41,30.97,59.54]
			nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Cobalt57","Barium133", "Americium241"]
			#fine-tune the range so that individual peaks are picked out
			if not traceInteg:
				fitLowArr=[0.0,0.36,0.3,0.06,0.15,0.32]
				fitHighArr=[0.25,0.39,0.42,0.1,0.25,0.4]
				binsize=[0.01,0.002,0.01,0.007,0.005,0.004] 
			else:
				fitLowArr=[2,12,12,0,0,9]
				fitHighArr=[np.inf,np.inf,np.inf,2,np.inf,np.inf]
				#binsize=np.zeros(len(energyList))
				binsize=np.ones(len(energyList))*0.25
		if chip==2 :
			energyList=[22.16, 88.03,122.06,14.41,30.97]
			nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Cobalt57","Barium133"]
			#fine-tune the range so that individual peaks are picked out
			if not traceInteg:
				fitLowArr=[0.0,0.44,0.48,0,0]
				fitHighArr=[0.25,0.6,0.55,0.15,0.3]
				binsize=[0.01,0.007,0.005,0.005,0.005]
			else:
				fitLowArr=[0,70,80,0,0]
				fitHighArr=[30,80,120,120,100]
				binsize=np.zeros(len(energyList))
	elif amp==2:
		with open(input2, 'r') as files2:
			fileList=files2.readlines()
		if chip==3:
			energyList=[14.41, 122.06, 59.54, 22.16, 88.03, 30.97]
			binsize=np.zeros(len(energyList))
			nameList=["Cobalt57","Cobalt57", "Americium241", "Cadmium109", "Cadmium109","Barium133"]
			#fine-tune the range so that individual peaks are picked out
			if not traceInteg:
				fitLowArr=[0.03,0.25,0.18,0.06,0.25,0.1]
				fitHighArr=[0.08,0.35,0.24,0.12,0.28,0.15]
			else:
				fitLowArr=[0,45,20,2,33,4]
				fitHighArr=[20,60,30,15,43,20]
		elif chip==4:
			energyList=[22.16, 88.03, 14.41, 122.06]
			binsize=np.zeros(len(energyList))
			nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Cobalt57"]
			if not traceInteg:	
				fitLowArr=[0.02,0.18,0.0,0.24]
				fitHighArr=[0.05,1,0.07,0.26]
			else:
				fitLowArr=[0,0,0,0]
				fitHighArr=[30,80,10,90]
	else:
		print("Choose amp1 or amp2")
		fileList,energyList,nameList,fitLowArr, fitHighArr = [],[],[],[],[]

	fileList=fixFileList(fileList)
	return fileList, energyList, nameList, fitLowArr, fitHighArr, binsize

def getEdgeFiles(amp):
	os.chdir(sys.path[0])
	if amp==1:
		with open(input1_edge, 'r') as files1:
			fileList=files1.readlines()
		energyList=[39.46]
		nameList=["Cobalt57"]
		if not traceInteg:
			fitLowArr=[0.13]
			fitHighArr=[0.17]
		else:
			fitLowArr=[20]
			fitHighArr=[35]
	elif amp==2:
		with open(input2_edge, 'r') as files2:
			fileList=files2.readlines()
		energyList=[39.46]
		nameList=["Cobalt57"]
		if not traceInteg:
			fitLowArr=[0.15]
			fitHighArr=[0.19]
		else:
			fitLowArr=[7]
			fitHighArr=[25]
	else:
		print("Choose amp1 or amp2")
		fileList,energyList,nameList,fitLowArr,fitHighArr = [],[],[],[],[]
	
	if len(fileList)==0: #no edge files used
		return [],[],[],[],[]
	
	fileList=fixFileList(fileList)
	return fileList, energyList, nameList, fitLowArr, fitHighArr

	
def energyCalibFit(trueEn, data, err, dataName, saveto, integral):
	plt.clf()
	
	#get name for plotting
	dataNameStr=dataName.replace(" ", "")	
		
	#organize arrays for interpolation
	trueEn_sorted = sorted(trueEn)
	#sort by energy
	data_sorted = [x for _, x in sorted(zip(trueEn, data), key=lambda pair: pair[0])]
	err_sorted = [x for _, x in sorted(zip(trueEn, err), key=lambda pair: pair[0])]

	
	#Plot data and fit functions
	x = np.linspace(np.min(trueEn), np.max(trueEn), 1000)
	xx = np.linspace(np.min(data), np.max(data), 1000)
	plt.errorbar(trueEn, data, yerr=err, fmt='o', label="data")
	
	if fit==0:
		coef, coef_pcov = curve_fit(linFit,trueEn,data,sigma=err,absolute_sigma=True) #linear 
		linPlot=linFit(x,*coef)
		plt.plot(x, linPlot, '--k',label=f"Linear fit")
		print(f"Peak linear fit: y={coef[0]:.3f}x+{coef[1]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		#invertedFn=interpolate.interp1d(linPlot,x,fill_value="extrapolate")
		invertedFn=interpolate.interp1d(linPlot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-1)#loose 1 dof 
		res_var = sum_square/ndof
		fn="linear"
	elif fit==1:
		coef2, coef2_pcov = curve_fit(quadFit,trueEn,data,sigma=err,absolute_sigma=True) #quadratic
		quadPlot=quadFit(x,*coef2)
		plt.plot(x, quadPlot, '--r',label=f"y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
		print(f"Peak quadratic fit: y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(quadPlot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-2)#loose 2 dof 
		res_var = sum_square/ndof
		fn="quadratic"
	elif fit==2:
		coef3, coef3_pcov = curve_fit(triFit,trueEn,data,sigma=err,absolute_sigma=True) #3rd deg poly
		triPlot=triFit(x,*coef3)
		plt.plot(x, triPlot, '--b',label=f"3rd deg. polynomial")
		print(f"Peak 3rd deg poly fit: y={coef3[0]:.10f}x$^3$+{coef3[1]:.5f}x$^2$+{coef3[2]:.5f}x + {coef3[3]:.5f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(triPlot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-3)#loose 3 dof 
		res_var = sum_square/ndof
		fn="3rd degree polynomial"
	elif fit==3:
		popt, pcov = curve_fit(sqrtFit, trueEn, data, sigma=err,absolute_sigma=True) #square root
		sqrt_fn = sqrtFit(x, *popt)
		plt.plot(x,sqrt_fn,'--g',label=f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
		print(f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(sqrt_fn,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-2)#loose 2 dof 
		res_var = sum_square/ndof
		fn="sqrt"
	elif fit==4:
		tck1 = interpolate.splrep(trueEn_sorted, data_sorted,k=1) #linear spline
		spline1Plot = interpolate.splev(x, tck1)
		plt.plot(x, spline1Plot, '--g', label="linear spline")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(spline1Plot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-3)#loose 3 dof for continuinty requirements
		res_var = sum_square/ndof
		fn="spline,k=1"
	elif fit==5:	
		tck = interpolate.splrep(trueEn_sorted, data_sorted) #cubic spline
		splinePlot = interpolate.splev(x, tck)
		plt.plot(x, splinePlot, '--g', label="Cubic Spline")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(splinePlot,x,fill_value="extrapolate")		
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-3)#loose 3 dof for continuinty requirements
		res_var = sum_square/ndof
		fn="spline,k=3"
	elif fit==6: 
		#piecewise linear, floating breakpoint
		coef, coef_pcov = curve_fit(piecewise_linear, trueEn, data,sigma=err,absolute_sigma=True,p0=[50,0.17,0,0])
		piecew=piecewise_linear(x,*coef)
		print(coef)
		plt.plot(x, piecew, label="Piecewise linear")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(piecew,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, data_sorted, err_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-4)#loose 4 dof 
		res_var = sum_square/ndof
		fn="piecewise"


	plt.xlabel("True Energy [keV]")
	if integral:
		plt.ylabel(f"{dataName} (from integral)")
	else:
		plt.ylabel(f"{dataName} (from peak height)")	
	plt.legend(loc="best")
	plt.grid()
	plt.xlim([-10,1.1*np.max(trueEn)])
	plt.ylim([-0.05,1.1*np.max(data)])
	if not savePlots:
		plt.show()
	elif savePlots and integral:
		plt.savefig(f"{saveto}integral_{dataNameStr}.pdf")
	elif savePlots and not integral:
		plt.savefig(f"{saveto}peaks_{dataNameStr}.pdf") 
	plt.clf()
	
				
	#print/save params
	if savePlots:
		if integral:			
			saveto=f"{saveto}integral_fit{fit}_params.txt"
		else:
			saveto=f"{saveto}peaks_fit{fit}_params.txt"
		k=open(saveto, "w")
		k.write("ODR fit (x=measured energy [V], y=true energy [keV])\n")
		k.write("chi2/ndf = %0.3f/%d = %0.3f" %(sum_square,ndof,res_var) +"\n")
		k.write(fn+"\n")
		k.close()
		#Display contents to terminal
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		print("Fit (x=measured energy [V], y=true energy [keV])")
		print("chi2/ndf = %0.3f/%d = %0.3f" %(sum_square,int(sum_square/res_var)+1,res_var))
		print(fn)
		
	return invertedFn
	
	
########################################################################################
###########################################################################
##############################################################

#files to be used for energy calibration curve
enResArr, muArr, sigmaArr, nArr, muErrArr, sigErrArr = [],[],[],[],[],[]


if fitSpectra:
	fileList,energyList,nameList,fitLowArr,fitHighArr,binsize = getFiles(pix)

	#loop through all files, fit with Gaussian, return mean/sigma/energy resolution and store in arrays - separately for each amp
	i=0
	for file in fileList:
		settings=[file, nameList[i], pix, energyList[i], savePlots, chip]
		popt, enRes, pcov = enResFitting.enResPlot(settings,fitLow=fitLowArr[i], fitHigh=fitHighArr[i], savedir=dataDir, integral=traceInteg, binSize=binsize[i])
		enResFitting.printParams(settings, popt, enRes, pcov, savedir=dataDir, integral=traceInteg)
		enResArr.append(enRes)
		muArr.append(popt[1])
		sigmaArr.append(popt[2])
		muErrArr.append(np.sqrt(np.diag(pcov))[1])
		sigErrArr.append(np.sqrt(np.diag(pcov))[2])
		i+=1

	#loop through additional files, fit with integrated Gaussian (for Compton edge), return mean/sigma/energy resolution and store in arrays
	edgeFileList,edgeEnergyList,edgeNameList,edgeFitLowArr,edgeFitHighArr = getEdgeFiles(pix)
	energyList.extend(edgeEnergyList)

	i=0
	for file in edgeFileList:
		settings=[file, edgeNameList[i], pix, edgeEnergyList[i],savePlots,chip]
		popt, enRes, pcov, integ= enResFitting.enResPlot(settings,edge=True,fitLow=edgeFitLowArr[i], fitHigh=edgeFitHighArr[i], savedir=dataDir)
		enResFitting.printParams(settings, popt, enRes, pcov, edge=True, savedir=dataDir)
		enResArr.append(enRes)
		muArr.append(popt[2])
		sigmaArr.append(popt[3])
		muErrArr.append(np.sqrt(np.diag(pcov))[2])
		sigErrArr.append(np.sqrt(np.diag(pcov))[3])
		i+=1
		
else: #if spectra have been fit before, pull out values from txt files
	print("spectra are already fit")	
	energyList, muArr, sigmaArr, enResArr, muErrArr, sigErrArr, enresErrArr1= enResFitting.getVals_fromTxt(dataDir, traceInteg)


#use fit mean of measured peaks and associated error to create calibration curve
fitFn = energyCalibFit(energyList, muArr, muErrArr, "Fit Mean [V]", saveDir, traceInteg)


#use calibration curve to calibrate a spectrum
if pix==1:
	if chip==1:
	#chip001 v2#not fitting edge
		files=["050422_amp1/35mV_cadmium109_5min.h5py","050422_amp1/200mV_cadmium109_960min.h5py","050422_amp1/100mV_cobalt57_180min.h5py","050422_amp1/35mV_cobalt57_60min.h5py","050522_amp1/test_100mV_barium133_330min.h5py","051822_amp1/125mV_digitalPaired_americium241_1020min.h5py"]
		name=["Cadmium109-calib","Cadmium109-calib","Cobalt57-calib","Cobalt57-calib","Barium133-calib", "Americium241-calib"]
		trueEn=[22.16,88.03,122.06,14.41,30.97,59.54]
		binSizeArr=[0.5,1.2,2.5,0.5,1.5,1.0]#default 0
		if traceInteg:
			fitLowArr=[18,60,100,0,0,50]
			fitHighArr=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
		else:
			fitLowArr=[10,70,100,0,0,55] #default 0
			fitHighArr=[40,np.inf,140,np.inf,np.inf,np.inf] #default np.inf
			
		for i,f in enumerate(files):
			settings=[homeDir+f,name[i],pix,trueEn[i],savePlots,chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings,fit=fit,coef=fitFn,fitLow=fitLowArr[i], fitHigh=fitHighArr[i], binSize=binSizeArr[i], integral=traceInteg)
			enResFitting.printParams(settings, popt, enRes, pcov, integral=traceInteg)

	if chip==2:
	#chip002 v2#not fitting edge
		files=["030322_amp1/TEST_chip2_cadmium109_10min.h5py","030822_amp1/chip2_200mV_cadmium109_330min_combined.h5py","030422_amp1/chip2_150mV_cobalt57_150min.h5py","030422_amp1/chip2_75mV_cobalt57_120min.h5py","030722_amp1/chip2_100mV_barium133_960min.h5py"]
		name=["Cadmium109-calib","Cadmium109-calib","Cobalt57-calib","Cobalt57-calib","Barium133-calib"]
		trueEn=[22.16,88.03,122.06,14.41,30.97]
		binSizeArr=[0.5,1.2,2.5,0.5,1]#default 0
		if traceInteg:
			fitLowArr=[0,0,100,0,0]
			fitHighArr=[np.inf,np.inf,np.inf,np.inf,np.inf]
		else:
			fitLowArr=[0,84,100,0,0] #default 0
			fitHighArr=[40,np.inf,np.inf,np.inf,np.inf] #default np.inf
			
		for i,f in enumerate(files):
			settings=[homeDir+f,name[i],pix,trueEn[i],savePlots,chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings,fit=fit,coef=fitFn,fitLow=fitLowArr[i], fitHigh=fitHighArr[i], binSize=binSizeArr[i], integral=traceInteg)
			enResFitting.printParams(settings, popt, enRes, pcov, integral=traceInteg)
	
	if chip==3:
	#chip003 v1
		#not fitting edge
		files=["110421_amp1/Americium_480min_combined.h5py","102021_amp1/cadmium109_45min.h5py","102821_amp1/cadmium109_16h.h5py","110821_amp1/barium133_combined_65min.h5py","102021_amp1/cobalt57_14h.h5py","102021_amp1/cobalt57_14h.h5py"]
		name=["Americium241-calib","Cadmium109-calib","Cadmium109-calib","Barium133-calib","Cobalt57-calib","Cobalt57-calib"]
		trueEn=[59.54,22.16,88.03,30.97,122.06,14.41]
		binSizeArr=[0,0,1,1,1,0,0]#default 0
		if traceInteg:
			fitLowArr=[50,0,0,60,100,0]
			fitHighArr=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]
		else:
			fitLowArr=[50,0,60,0,98,0,35] #default 0
			fitHighArr=[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,60] #default np.inf
			
		for i,f in enumerate(files):
			settings=[homeDir+f,name[i],pix,trueEn[i],savePlots,chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings,fit=fit,coef=fitFn,fitLow=fitLowArr[i], fitHigh=fitHighArr[i], binSize=binSizeArr[i], integral=traceInteg)
			enResFitting.printParams(settings, popt, enRes, pcov, integral=traceInteg)
	
	elif chip==4:		
	#chip004 v1
		files=["120921_amp1/15mV_chip004_AC_Cadmium_120min.h5py","120921_amp1/90mV_chip004_AC_Cadmium_1200min.h5py","120721_amp1/lowPeak_chip004_AC_cobalt57_60min.h5py","120721_amp1/highPeak_chip004_AC_cobalt57_960min.h5py"]
		name=["Cadmium109-calib","Cadmium109-calib","Cobalt57-calib","Cobalt57-calib"]
		trueEn=[22.16,88.03,14.41,122.06]
		binSizeArr=[0,1,1,2]
		if traceInteg:
			fitLowArr=[0,80,0,110]
			fitHighArr=[np.inf,np.inf,np.inf,np.inf]
		else:
			fitLowArr=[0,80,0,110] #default 0
			fitHighArr=[np.inf,np.inf,np.inf,np.inf]
			
		for i,f in enumerate(files):
			settings=[homeDir+f,name[i],pix,trueEn[i],savePlots,chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings,fit=fit,coef=fitFn,fitLow=fitLowArr[i], binSize=binSizeArr[i], integral=traceInteg)
			enResFitting.printParams(settings, popt, enRes, pcov, integral=traceInteg)
		
					
elif pix==2:
	if chip==3:
	#chip003 v1
		files=["111521_amp2/overnight_Americium241_960min.h5py","111621_amp2/day_Cadmium109_300min.h5py","111221_amp2/weekend_Cobalt57_4020min.h5py","111221_amp2/weekend_Cobalt57_4020min.h5py","120221_amp2/calib_barium133_180min.h5py","120221_amp2/calib_cadmium190_1080min.h5py"]
		name=["Americium241-calib","Cadmium109-calib","Cobalt57-calib","Cobalt57-calib","Barium133-calib", "Cadmium109-calib"]
		trueEn=[59.54,22.16,122.06,14.41,30.97, 88.03]
		binSizeArr=[0,0,1,0,0,1] #default 0
		if traceInteg:
			fitLowArr=[10,15,110,0,0,0]
			fitHighArr=[14,np.inf,np.inf,25,np.inf,np.inf]
		else:
			fitLowArr=[48,17,100,0,20,60] #default 0
			fitHighArr=[np.inf,30,140,np.inf,45,60] #default np.inf
			
		for i,f in enumerate(files):
			settings=[homeDir+f,name[i],pix,trueEn[i],savePlots,chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings,fit=fit,coef=fitFn,fitLow=fitLowArr[i], fitHigh=fitHighArr[i], binSize=binSizeArr[i],integral=traceInteg)
			enResFitting.printParams(settings, popt, enRes, pcov, integral=traceInteg)
	

	elif chip==4:
	#chip004 v1
		files=["121321_amp2/15mV_chip004_AC_cadmium109_120min.h5py","121321_amp2/90mV_chip004_AC_cadmium109_1200min.h5py","121521_amp2/15mV_chip004_cobalt57_60min.h5py", "011022_amp2/chip004_cobalt57_combined_2220min.h5py"]
		name=["Cadmium109-calib","Cadmium109-calib","Cobalt57-calib","Cobalt57-calib"]
		trueEn=[22.16,88.03,14.41, 122.06]
		binSizeArr=[0,1,1,2] #default 0
		if traceInteg:
			fitLowArr=[0,60,0,100]
			fitHighArr=[np.inf,np.inf,np.inf,np.inf]
		else:
			#fitLowArr=[17,82,0,116] #default 0
			fitLowArr=[17,50,0,80] #default 0
			fitHighArr=[30,95,20,np.inf]
			
		for i,f in enumerate(files):
			settings=[homeDir+f,name[i],pix,trueEn[i],savePlots,chip]
			popt, enRes, pcov = enResFitting.enResPlot(settings,fit=fit,coef=fitFn,fitLow=fitLowArr[i], fitHigh=fitHighArr[i],binSize=binSizeArr[i],integral=traceInteg)
			enResFitting.printParams(settings, popt, enRes, pcov, integral=traceInteg)
			