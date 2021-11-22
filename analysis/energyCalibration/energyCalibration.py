import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
from scipy import stats, interpolate
from scipy.odr import ODR, Model, Data, RealData
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

def getSumSq(trueEn, voltage, err, coef):
	sum_square=0
	for i in range(len(trueEn)):
		fitVal=interpolate.splev(voltage[i], coef)
		num = (trueEn[i]-fitVal)**2
		fitErr=interpolate.splev(err[i],coef)
		denom = fitErr**2
		sum_square += (num/denom)
		i+=1
	print(f"final sum_square {sum_square}")
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
		fitLow_p=[0.05,0.25,0.02,0.29,0.17,0.09]
		fitLow_i=[150,725,0,1700,300,200,100]
	elif amp==2:
		with open(input2, 'r') as files2:
			fileList=files2.readlines()
		energyList=[14.41, 122.06, 59.54, 22.16]
		nameList=["Cobalt57","Cobalt57", "Americium241", "Cadmium109"]
		#fine-tune the range so that individual peaks are picked out
		fitLow_p=[0.03,0.3,0.19,0.06]
		fitLow_i=[0,900,200,0]
	else:
		print("Choose amp1 or amp2")
		fileList,energyList,nameList,fitLow_p,fitLow_i = [],[],[],[],[]

	fileList=fixFileList(fileList)
	return fileList, energyList, nameList, fitLow_p, fitLow_i

		
def getEdgeFiles(amp):
	os.chdir(sys.path[0])
	if amp==1:
		with open(input1_edge, 'r') as files1:
			fileList=files1.readlines()
		energyList=[39.46]
		nameList=["Cobalt57"]
		fitLow_p=[0.13]
		fitHigh_p=[0.17]
		fitLow_i=[0]
		fitHigh_i=[1000]
	elif amp==2:
		with open(input2_edge, 'r') as files2:
			fileList=files2.readlines()
		energyList=[39.46]
		nameList=["Cobalt57"]
		fitLow_p=[0.15]
		fitHigh_p=[0.19]
		fitLow_i=[100]
		fitHigh_i=[200]
	else:
		print("Choose amp1 or amp2")
		fileList,energyList,nameList,fitLow_p,fitLow_i, fithith_p, fitHigh_i = [],[],[],[],[],[],[]
	
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
	trueEn_sorted= sorted(trueEn)
	amp_p_sorted = [x for _, x in sorted(zip(trueEn, amp_p), key=lambda pair: pair[0])]
	err_p_sorted = [x for _, x in sorted(zip(trueEn, err_p), key=lambda pair: pair[0])]
	amp_i_sorted = [x for _, x in sorted(zip(trueEn, amp_i), key=lambda pair: pair[0])]
	err_i_sorted = [x for _, x in sorted(zip(trueEn, err_i), key=lambda pair: pair[0])]

	
	#Plot data and fit functions
	x = np.linspace(np.min(trueEn), np.max(trueEn), 100)
	plt.errorbar(trueEn, amp_p, yerr=err_p, fmt='o', label="data")
	
	if fit==0:
		coef, coef_pcov = curve_fit(linFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #linear 
		linPlot=linFit(x,*coef)
		plt.plot(x, linPlot, '--k',label=f"Linear fit")
		print(f"Peak linear fit: y={coef[0]:.3f}x+{coef[1]:.3f}")
	elif fit==1:
		coef2, coef2_pcov = curve_fit(quadFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #quadratic
		quadPlot=quadFit(x,*coef2)
		plt.plot(x, quadPlot, '--r',label=f"y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
		print(f"Peak quadratic fit: y={coef2[0]:.5f}x$^2$+{coef2[1]:.3f}x+{coef2[2]:.3f}")
	elif fit==2:
		coef3, coef3_pcov = curve_fit(triFit,trueEn,amp_p,sigma=err_p,absolute_sigma=True) #3rd deg poly
		triPlot=triFit(x,*coef3)
		plt.plot(x, triPlot, '--b',label=f"3rd deg. polynomial")
		print(f"Peak 3rd deg poly fit: y={coef3[0]:.3f}x$^3$+{coef3[1]:.3f}x$^2$+{coef3[2]:.3f}x + {coef3[3]:.3f}")
	elif fit==3:
		popt, pcov = curve_fit(sqrtFit, trueEn, amp_p, sigma=err_p,absolute_sigma=True) #square root
		sqrt_fn = sqrtFit(x, *popt)
		plt.plot(x,sqrt_fn,'--g',label=f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
		print(f"y={popt[0]:.3f}*sqrt(x)+{popt[1]:.3f}")
	elif fit==4:
		tck1 = interpolate.splrep(trueEn_sorted, amp_p_sorted,k=1) #linear spline
		spline1Plot = interpolate.splev(x, tck1)
		plt.plot(x, spline1Plot, '--g', label="linear spline")
	elif fit==5:	
		tck = interpolate.splrep(trueEn_sorted, amp_p_sorted) #cubic spline
		splinePlot = interpolate.splev(x, tck)
		plt.plot(x, splinePlot, '--g', label="Cubic Spline")

	plt.xlabel("True Energy [keV]")
	plt.ylabel(f"{dataName} (from peak height)")	
	#plt.ylabel(f"{dataName} (from integral)")
	plt.legend(loc="best")
	plt.grid()
	plt.xlim([-10,1.1*np.max(trueEn)])
	plt.ylim([-0.05,1.1*np.max(amp_p)])
	plt.savefig(f"{saveto}peaks_{dataNameStr}.pdf") if savePlots else plt.show()
	plt.clf()
	
	
	#fit opposite orientation so that scaling is easier
	#need different regression technique for X-error bars
	if fit==0:
		datain = RealData(amp_p, trueEn, sx=err_p)
		coef_fit, res_var, sum_square = odr_polyfit(datain,1)
		fn="m*x + b"
	elif fit==1:
		datain = RealData(amp_p, trueEn, sx=err_p)
		coef_fit, res_var, sum_square = odr_polyfit(datain,2)
		fn="a * x *x + b*x + c"
	elif fit==2:
		datain = RealData(amp_p, trueEn, sx=err_p)
		coef_fit, res_var, sum_square = odr_polyfit(datain,3)
		fn="a*x*x*x + b*x*x + c*x + d"
	elif fit==3:
		datain = RealData(amp_p, trueEn, sx=err_p)
		sqrt_model = Model(sqrt_odr)
		odr_sqrt = ODR(datain, sqrt_model,[1e-8,1e-8])
		out_sqrt=odr_sqrt.run()
		coef_fit=out_sqrt.beta
		sum_square=out_sqrt.sum_square
		res_var=out_sqrt.res_var
		fn="A*np.sqrt(x)+B"
	elif fit==4:
		coef_fit = interpolate.splrep(amp_p_sorted, trueEn_sorted, w=err_p_sorted, k=1)
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, coef_fit)
		ndof = float(len(trueEn)-3-1)#loose 3 dof for continuinty requirements - need extra 1 for poor counting!?
		print(f"len(trueEn) {len(trueEn)}")
		res_var = sum_square/ndof
		fn="spline, k=1"
	elif fit==5:
		coef_fit = interpolate.splrep(amp_p_sorted, trueEn_sorted, w=err_p_sorted)
		sum_square=getSumSq(trueEn_sorted, amp_p_sorted, err_p_sorted, coef_fit)
		ndof = float(len(trueEn)-3-1)#loose 3 dof for continuinty requirements - need extra 1 for poor counting!?
		res_var = sum_square/float(ndof)
		fn="spline,k=3"
	
	#print/save params
	if savePlots:
		saveto=f"{saveto}fit{fit}_params.txt"
		k=open(saveto, "w")
		k.write("ODR fit (x=measured energy [V], y=true energy [keV])\n")
		k.write("chi2/ndf = %0.3f/%d = %0.3f" %(sum_square,int(sum_square/res_var)+1,res_var) +"\n")
		k.write(fn+"\n")
		k.write("Coefficients: \n")
		k.write(str(coef_fit))
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
		print("Coefficients:")
		print(coef_fit)
		
	return coef_fit
	
	
	
########################################################################################
###########################################################################
##############################################################

#files to be used for energy calibration curve
#amp1
enResArr1, muArr1, sigmaArr1, nArr1 = [],[],[], []

#amp2
enResArr2, muArr2, sigmaArr2, nArr2 = [],[],[], []


if fitSpectra:
	fileList,energyList,nameList,fitLow_p,fitLow_i = getFiles(pix)

	#loop through all files, fit with Gaussian, return mean/sigma/energy resolution and store in arrays - separately for each amp
	i=0
	for file in fileList:
		settings=[file, nameList[i], pix, energyList[i], savePlots]
		popt, enRes, pcov, integ = enResFitting.enResPlot(settings,fitLow=fitLow_p[i])
		enResFitting.printParams(settings, integ, popt, enRes, pcov)
		#integral argument = integral bin size (in V)
		#AMANDA - N calculation is failing for integral
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
		settings=[file, edgeNameList[i], pix, edgeEnergyList[i],savePlots]
		popt, enRes, pcov, integ= enResFitting.enResPlot(settings,edge=True,fitLow=edgeFitLow_p[i], fitHigh=edgeFitHigh_p[i])
		enResFitting.printParams_edge(settings, integ, popt, enRes, pcov)
		poptI, enResI, pcovI, integI = enResFitting.enResPlot(settings,edge=True,fitLow=edgeFitLow_i[i], fitHigh=edgeFitHigh_i[i],integral=20)
		enResFitting.printParams_edge(settings, integI, poptI, enResI, pcovI, integral=True)
		enResArr1.append([enRes,enResI])
		muArr1.append([popt[2],poptI[2]])
		sigmaArr1.append([popt[3],poptI[3]])
		nArr1.append([integ,integI])
		i+=1
		
else: #if spectra have been fit before, pull out values from txt files
	print("spectra are already fit")	
	energyList, muArr1, sigmaArr1, nArr1, enResArr1 = enResFitting.getVals_fromTxt(dataDir)



#calculate error
errArr1 = enResFitting.calcError(sigmaArr1, nArr1)


"""
#error Arr1ay = sigma/sqrt(N) (for edge, 2sig integral from mu)
for j in range(len(sigmaArr1)):
	err_p=sigmaArr1[j][0]/np.sqrt(nArr1[j][0])
	err_i=sigmaArr1[j][1]/np.sqrt(nArr1[j][1])
	errArr1.append([err_p,err_i])
"""

#use fit mean of measured peaks and associated error to create calibration curve
coef_p = energyCalibFit(energyList, muArr1, errArr1, "Fit Mean [V]", saveDir)


#use calibration curve to calibrate a spectrum
file="110421_amp1/Americium_480min_combined.h5py"
settings=[homeDir+file, "Americium241-calib", 1, 59.54, savePlots]
popt, enRes, pcov, integ = enResFitting.enResPlot(settings,fit=fit,coef=coef_p,fitLow=50)
#popt, enRes, pcov = enResFitting.enResPlot_scale(settings,coef_p,fit,fitLow=50)
enResFitting.printParams(settings, -1, popt, enRes, pcov)
	

file="102021_amp1/cadmium109_45min.h5py"
settings=[homeDir+file,  "Cadmium109-calib", 1, 22.16, savePlots]
popt, enRes, pcov, integ = enResFitting.enResPlot(settings,coef=coef_p,fit=fit)
enResFitting.printParams(settings, -1, popt, enRes, pcov)

		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		