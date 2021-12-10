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
pix=int(pix)
fit=int(fit)


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
	
"""
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
	res_var=out.res_var
	sum_square=out.sum_square
	return coef, res_var, sum_square
"""
	
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
	if amp==1:
		with open(input1, 'r') as files1:
			fileList=files1.readlines()
		energyList=[22.16, 88.03, 14.41, 122.06, 59.54, 30.97]
		nameList=["Cadmium109", "Cadmium109", "Cobalt57", "Cobalt57", "Americium241", "Barium133"]
		#fine-tune the range so that individual peaks are picked out
		fitLow_p=[0.05,0.244,0.02,0.3,0.19,0.09]
		#fitLow_i=[150,700,0,1700,350,100]
		fitLow_i=[0,70,0,70,40,10]
		fitHigh_p=[0.11,1,0.07,0.32,1,0.15]
		#fitHigh_i=[400,1000,250,2000,600,300]
		fitHigh_i=[30,80,100,80,60,100]
	elif amp==2:
		with open(input2, 'r') as files2:
			fileList=files2.readlines()
		energyList=[14.41, 122.06, 59.54, 22.16, 88.03, 30.97]
		nameList=["Cobalt57","Cobalt57", "Americium241", "Cadmium109", "Cadmium109","Barium133"]
		#fine-tune the range so that individual peaks are picked out
		fitLow_p=[0.03,0.3,0.18,0.06,0.25,0.1]
		#fitLow_i=[0,450,200,0,75]
		fitLow_i=[0,45,20,2,33,4]
		fitHigh_p=[0.08,0.33,0.24,0.12,0.28,0.15]
		#fitHigh_i=[100,1000,300,200,175]
		fitHigh_i=[100,60,30,15,43,20]
	else:
		print("Choose amp1 or amp2")
		fileList,energyList,nameList,fitLow_p,fitLow_i, fitHigh_p, fitHigh_i = [],[],[],[],[],[],[]

	fileList=fixFileList(fileList)
	return fileList, energyList, nameList, fitLow_p, fitLow_i, fitHigh_p, fitHigh_i

def getEdgeFiles(amp):
	os.chdir(sys.path[0])
	if amp==1:
		with open(input1_edge, 'r') as files1:
			fileList=files1.readlines()
		energyList=[39.46]
		nameList=["Cobalt57"]
		fitLow_p=[0.13]
		fitHigh_p=[0.17]
		fitLow_i=[20]
		fitHigh_i=[35]
	elif amp==2:
		with open(input2_edge, 'r') as files2:
			fileList=files2.readlines()
		energyList=[39.46]
		nameList=["Cobalt57"]
		fitLow_p=[0.15]
		fitHigh_p=[0.19]
		fitLow_i=[7]
		fitHigh_i=[25]
	else:
		print("Choose amp1 or amp2")
		fileList,energyList,nameList,fitLow_p,fitLow_i, fitHigh_p, fitHigh_i = [],[],[],[],[],[],[]
	
	fileList=fixFileList(fileList)
	return fileList, energyList, nameList, fitLow_p, fitLow_i, fitHigh_p, fitHigh_i

	
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
	trueEn_sorted = sorted(trueEn)
	#sort by energy
	amp_p_sorted = [x for _, x in sorted(zip(trueEn, amp_p), key=lambda pair: pair[0])]
	err_p_sorted = [x for _, x in sorted(zip(trueEn, err_p), key=lambda pair: pair[0])]
	amp_i_sorted = [x for _, x in sorted(zip(trueEn, amp_i), key=lambda pair: pair[0])]
	err_i_sorted = [x for _, x in sorted(zip(trueEn, err_i), key=lambda pair: pair[0])]

	
	#Plot data and fit functions
	x = np.linspace(np.min(trueEn), np.max(trueEn), 1000)
	xx = np.linspace(np.min(amp_p), np.max(amp_p), 1000)
	plt.errorbar(trueEn, amp_p, yerr=err_p, fmt='o', label="data")
	
	if fit==0:
		coef, coef_pcov = curve_fit(linFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #linear 
		linPlot=linFit(x,*coef)
		plt.plot(x, linPlot, '--k',label=f"Linear fit")
		print(f"Peak linear fit: y={coef[0]:.3f}x+{coef[1]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		#invertedFn=interpolate.interp1d(linPlot,x,fill_value="extrapolate")
		invertedFn=interpolate.interp1d(linPlot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-1)#loose 1 dof 
		res_var = sum_square/ndof
		fn="linear"
	elif fit==1:
		coef2, coef2_pcov = curve_fit(quadFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #quadratic
		quadPlot=quadFit(x,*coef2)
		plt.plot(x, quadPlot, '--r',label=f"y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
		print(f"Peak quadratic fit: y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(quadPlot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-2)#loose 2 dof 
		res_var = sum_square/ndof
		fn="quadratic"
	elif fit==2:
		coef3, coef3_pcov = curve_fit(triFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #3rd deg poly
		triPlot=triFit(x,*coef3)
		plt.plot(x, triPlot, '--b',label=f"3rd deg. polynomial")
		print(f"Peak 3rd deg poly fit: y={coef3[0]:.3f}x$^3$+{coef3[1]:.3f}x$^2$+{coef3[2]:.3f}x + {coef3[3]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(triPlot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-3)#loose 3 dof 
		res_var = sum_square/ndof
		fn="3rd degree polynomial"
	elif fit==3:
		popt, pcov = curve_fit(sqrtFit, trueEn, amp_p, sigma=err_p,absolute_sigma=True) #square root
		sqrt_fn = sqrtFit(x, *popt)
		plt.plot(x,sqrt_fn,'--g',label=f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
		print(f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(sqrt_fn,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-2)#loose 2 dof 
		res_var = sum_square/ndof
		fn="sqrt"
	elif fit==4:
		tck1 = interpolate.splrep(trueEn_sorted, amp_p_sorted,k=1) #linear spline
		spline1Plot = interpolate.splev(x, tck1)
		plt.plot(x, spline1Plot, '--g', label="linear spline")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(spline1Plot,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-3)#loose 3 dof for continuinty requirements
		res_var = sum_square/ndof
		fn="spline,k=1"
	elif fit==5:	
		tck = interpolate.splrep(trueEn_sorted, amp_p_sorted) #cubic spline
		splinePlot = interpolate.splev(x, tck)
		plt.plot(x, splinePlot, '--g', label="Cubic Spline")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(splinePlot,x,fill_value="extrapolate")		
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-3)#loose 3 dof for continuinty requirements
		res_var = sum_square/ndof
		fn="spline,k=3"
	elif fit==6: 
		#piecewise linear, floating breakpoint
		coef, coef_pcov = curve_fit(piecewise_linear, trueEn, amp_p,sigma=err_p,absolute_sigma=True,p0=[50,0.17,0,0])
		piecew=piecewise_linear(x,*coef)
		print(coef)
		plt.plot(x, piecew, label="Piecewise linear")
		#invert - if value outside of calibrated range, assign it a value of -1 (underflow)
		invertedFn=interpolate.interp1d(piecew,x,fill_value="extrapolate")
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, invertedFn)
		ndof = np.float64(len(trueEn)-4)#loose 4 dof 
		res_var = sum_square/ndof
		fn="piecewise"


	plt.xlabel("True Energy [keV]")
	plt.ylabel(f"{dataName} (from peak height)")	
	#plt.ylabel(f"{dataName} (from integral)")
	plt.legend(loc="best")
	plt.grid()
	plt.xlim([-10,1.1*np.max(trueEn)])
	plt.ylim([-0.05,1.1*np.max(amp_p)])
	plt.savefig(f"{saveto}peaks_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
				
	#print/save params
	if savePlots:
		saveto=f"{saveto}fit{fit}_params.txt"
		k=open(saveto, "w")
		k.write("ODR fit (x=measured energy [V], y=true energy [keV])\n")
		k.write("chi2/ndf = %0.3f/%d = %0.3f" %(sum_square,ndof,res_var) +"\n")
		k.write(fn+"\n")
		#k.write("Coefficients: \n")
		#k.write(str(coef_fit))
		k.close()
		#Display contents to terminal
		m = open(saveto, "r")
		text = m.read()
		print(text)
		m.close()
	else:
		print("ODR fit (x=measured energy [V], y=true energy [keV])")
		print("chi2/ndf = %0.3f/%d = %0.3f" %(sum_square,int(sum_square/res_var)+1,res_var))
		print(fn)
		#print("Coefficients:")
		#print(coef_fit)
		
	#return coef_fit
	return invertedFn
	
	
########################################################################################
###########################################################################
##############################################################

#files to be used for energy calibration curve
enResArr1, muArr1, sigmaArr1, nArr1, muErrArr1 = [],[],[],[],[]


if fitSpectra:
	fileList,energyList,nameList,fitLow_p,fitLow_i,fitHigh_p, fitHigh_i = getFiles(pix)

	#loop through all files, fit with Gaussian, return mean/sigma/energy resolution and store in arrays - separately for each amp
	i=0
	for file in fileList:
		settings=[file, nameList[i], pix, energyList[i], savePlots]
		popt, enRes, pcov, integ = enResFitting.enResPlot(settings,fitLow=fitLow_p[i], fitHigh=fitHigh_p[i], savedir=dataDir)
		enResFitting.printParams(settings, integ, popt, enRes, pcov, savedir=dataDir)
		#integral argument = integral bin size (in V)
		#AMANDA - N calculation is failing for integral
		poptI, enResI, pcovI, integI = enResFitting.enResPlot(settings, fitLow=fitLow_i[i], fitHigh=fitHigh_i[i], integral=True, savedir=dataDir)
		enResFitting.printParams(settings, integI, poptI, enResI, pcovI, integral=True, savedir=dataDir)
		enResArr1.append([enRes,enResI])
		muArr1.append([popt[1],poptI[1]])
		sigmaArr1.append([popt[2],poptI[2]])
		nArr1.append([integ,integI])
		muErrArr1.append([np.sqrt(np.diag(pcov))[1],np.sqrt(np.diag(pcovI))[1]])
		i+=1

	#loop through additional files, fit with integrated Gaussian (for Compton edge), return mean/sigma/energy resolution and store in arrays
	edgeFileList,edgeEnergyList,edgeNameList,edgeFitLow_p,edgeFitLow_i,edgeFitHigh_p,edgeFitHigh_i = getEdgeFiles(pix)
	energyList.extend(edgeEnergyList)

	i=0
	for file in edgeFileList:
		settings=[file, edgeNameList[i], pix, edgeEnergyList[i],savePlots]
		popt, enRes, pcov, integ= enResFitting.enResPlot(settings,edge=True,fitLow=edgeFitLow_p[i], fitHigh=edgeFitHigh_p[i], savedir=dataDir)
		enResFitting.printParams(settings, integ, popt, enRes, pcov, edge=True, savedir=dataDir)
		poptI, enResI, pcovI, integI = enResFitting.enResPlot(settings,edge=True,fitLow=edgeFitLow_i[i], fitHigh=edgeFitHigh_i[i],integral=20, savedir=dataDir)
		enResFitting.printParams(settings, integI, poptI, enResI, pcovI, integral=True, edge=True, savedir=dataDir)
		enResArr1.append([enRes,enResI])
		muArr1.append([popt[2],poptI[2]])
		sigmaArr1.append([popt[3],poptI[3]])
		nArr1.append([integ,integI])
		muErrArr1.append([np.sqrt(np.diag(pcov))[1],np.sqrt(np.diag(pcovI))[1]])
		i+=1
		
else: #if spectra have been fit before, pull out values from txt files
	print("spectra are already fit")	
	energyList, muArr1, sigmaArr1, nArr1, enResArr1, muErrArr1 = enResFitting.getVals_fromTxt(dataDir)



#calculate error
#errArr1 = enResFitting.calcError(sigmaArr1, nArr1)
errArr1 = muErrArr1

#use fit mean of measured peaks and associated error to create calibration curve
coef_p = energyCalibFit(energyList, muArr1, errArr1, "Fit Mean [V]", saveDir)


#use calibration curve to calibrate a spectrum


if pix==1:
	file="110421_amp1/Americium_480min_combined.h5py"
	settings=[homeDir+file, "Americium241-calib", 1, 59.54, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,fit=fit,coef=coef_p,fitLow=50)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)

	file="102021_amp1/cadmium109_45min.h5py"
	settings=[homeDir+file,  "Cadmium109-calib", 1, 22.16, savePlots]
	#popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitLow=15, fitHigh=30)
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)


	file="102821_amp1/cadmium109_16h.h5py"
	settings=[homeDir+file,  "Cadmium109-calib", 1, 88.03, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitLow=60)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)		


	file="110821_amp1/barium133_combined_65min.h5py"
	settings=[homeDir+file,  "Barium133-calib", 1, 30.97, savePlots]
	#popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitLow=20, fitHigh=45)
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)		
	

	file="102021_amp1/cobalt57_14h.h5py"
	settings=[homeDir+file,  "Cobalt57-calib", 1, 122.06, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit,fitLow=98,binSize=1)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)	
	
	file="102021_amp1/cobalt57_14h.h5py"
	settings=[homeDir+file,  "Cobalt57-calib", 1, 14.41, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)	


	file="102021_amp1/cobalt57_14h.h5py"
	settings=[homeDir+file,  "Cobalt57-edge-calib", 1, 39.46, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit,edge=True, fitLow=35, fitHigh=60)
	enResFitting.printParams(settings, -1, popt, enRes, pcov, edge=True)
			
elif pix==2:
	file="111521_amp2/overnight_Americium241_960min.h5py"
	settings=[homeDir+file, "Americium241-calib", 2, 59.54, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,fit=fit,coef=coef_p,fitLow=48)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)
	
	file="111621_amp2/day_Cadmium109_300min.h5py"
	settings=[homeDir+file,  "Cadmium109-calib", 2, 22.16, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitLow=17,fitHigh=30)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)
	
	file="111221_amp2/weekend_Cobalt57_4020min.h5py"
	settings=[homeDir+file,  "Cobalt57-calib", 2, 122.06, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit,fitLow=100, fitHigh=140)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)
		
	file="111221_amp2/weekend_Cobalt57_4020min.h5py"
	settings=[homeDir+file,  "Cobalt57-calib", 2, 14.41, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitHigh=30)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)
		
	file="120221_amp2/calib_barium133_180min.h5py"
	settings=[homeDir+file,  "Barium133-calib", 2, 30.97, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitLow=20, fitHigh=45)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)			
		
	file="120221_amp2/calib_cadmium190_1080min.h5py"
	settings=[homeDir+file,  "Cadmium109-calib", 2, 88.03, savePlots]
	popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit, fitLow=50)
	enResFitting.printParams(settings, -1, popt, enRes, pcov)			


		
		
		
		
		
		
		