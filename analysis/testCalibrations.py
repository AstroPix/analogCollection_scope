import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
import sys,os,glob
sys.path.insert(1, 'energyCalibration/')
import enResFitting

###############################
# Global variables
###############################

savePlots=False #If false, get option to save each plot upon viewing. If true, plots are not displayed and all automatically saved
traceInteg=False
v1=False
chip=2 #Can be: 3, 4 [for v1], 1, 2 [for v2]
label = "integ" if traceInteg else "peaks"

#Choose what to plot
allPlots=False #overrides all other options
rawEnergyResolution=False #plot energy resolution 
rawMuSigma=False #plot components of energy calibration - mu and sigma
calibrationCurve=False #plot calibration curves
calibratedMuEnRes=False #compare calib functions - mean and enRes
#must have all calibration functions tested and saved
#AMANDA - make more robust in pulling fit info AND improve chi2 calculation
calibratedMuChi2=False #compare calib functions - ratio of calibrated/expected mean and fit chi2/ndof
calibratedSigma=False #compare calib functions - Gaussian sigma
calibratedMuEnResScan=True #For one calib function, plot calibrated vs expected energy for all sources with color bar indicating energy resolution. Optional ratio plot

###############################
# Helper functions
###############################

def sqrtFit(x,A,B):
	return A/np.sqrt(x)+B
	
def linFit(x,A,B):
	return A*x+B
	
def flattenIn(data):
	out=[]
	for evt in data:
		flatData=np.array(evt).flatten('F')
		out.append(flatData[0])
	return out
	
def calcChi2(true,meas,err,popt):
	sum_square=0
	A=popt[0]
	B=popt[1]
	for i in range(len(true)):
		fitVal=sqrtFit(meas[i],A,B)
		num = (true[i]-fitVal)**2
		denom = err[i]**2
		if denom==0:
			denom=1e-8
		sum_square += (num/denom)
		i+=1
	ndof = np.float64(len(true)-1)#loose 1 dof 
	res_var = sum_square/ndof	
	return res_var
	
	
def getChi2(indir):
	chiArr, ndofArr = [],[]
	fits=['linear','quad','tri','sqrt','spline1','spline3', 'piecewise']

	for fit in fits:
		os.chdir(indir+fit+'/')
		#AMANADA - make variable for peaks/integral with backwards compatibilithy
		fitFile = glob.glob('*fit*.txt') #ONLY POSSIBLE WHEN _ONE_ FIT FILE PRESENT (not both peaks and integral)
		openFile = open(fitFile[0],'r')
		lines=openFile.readlines()
		fitline=lines[1].split(' = ')
		vals=fitline[1].split('/')
		chiArr.append(float(vals[0]))
		ndofArr.append(float(vals[1]))
	return chiArr, ndofArr


def saveFromInput(saveFile): 
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	inp= 'a'
	while inp!='y' and inp!='n':
		print("Save plot? y/n")
		inp = input()
	savePlt=True if inp=='y' else False
	if savePlt: #save plt stored with plt.gcf()
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)
	
###############################
#Plotting functions - callable
###############################

def rawEnRes(dataDir,energyList, enResArr,energyList2=[],enResArr2=[]):
	plt.plot(energyList, enResArr, 'o', label="Amp1")
	plt.xlabel("True Energy [keV]")
	plt.ylabel(f"Energy Resolution [%]")
	if len(enResArr2)>0:
		plt.plot(energyList2, enResArr2, 'or', label="Amp2")	
		plt.legend(loc="best")
	plt.grid()
	saveto=f"{dataDir}{label}_comparePixels_enRes.pdf"
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)
	
def rawMuSig(dataDir,sigmaArr,sigErr,muArr,muErr,sigmaArr2=[],sigErr2=[],muArr2=[],muErr2=[]):
	sigmaArr=[sigmaArr,sigmaArr2]
	sigErr=[sigErr,sigErr2]
	muArr=[muArr,muArr2]
	muErr=[muErr,muErr2]
	colors=['--b','--r']
	
	for i,arr in enumerate(sigmaArr):
		if len(arr)==0:
			continue
		#sig=flattenIn(arr)
		x = np.linspace(np.min(muArr[i]), np.max(muArr[i]), 1000)
		plt.errorbar(muArr[i],sigmaArr[i], xerr=muErr[i], yerr=sigErr[i], fmt='o', label=f"Amp{i+1}")
		popt, pcov = curve_fit(sqrtFit, muArr[i], sigmaArr[i], absolute_sigma=True) #square root
		sqrt_fn = sqrtFit(x, *popt)
		chi2=calcChi2(sigmaArr[i],muArr[i],sigErr[i],popt)
		plt.plot(x,sqrt_fn,colors[i],label=f"{popt[0]:.3f}/sqrt(E)+{popt[1]:.3f}\n chi2/ndof={chi2:.4f}")
		#popt, pcov = curve_fit(linFit, mu1, sig1, absolute_sigma=True) #linear
		#lin_fn = linFit(x, *popt)
		#chi2=calcChi2(sig1,mu1,sigE1,popt)
		#plt.plot(x,lin_fn,':b',label=f"{popt[0]:.3f}E+{popt[1]:.3f}\n chi2/ndof={chi2:.4f}")
	plt.xlabel("Measured energy [V]")
	plt.ylabel(f"Measured sigma [V]")	
	plt.grid()
	plt.legend(loc="best")
	saveto=f"{dataDir}{label}_comparePixels_muVsSig.pdf"
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)

	
def calibCurve(dataDir,energyList,muArr,muErr,energyList2=[],muArr2=[],muErr2=[]):
	energyList=[energyList,energyList2]
	muArr=[muArr,muArr2]
	muErr=[muErr,muErr2]
	colors=["blue","red"]
	
	for i,arr in enumerate(muArr):
		if len(arr)==0:
			continue
		plt.errorbar(energyList[i], muArr[i], yerr=muErr[i], fmt='o', color=colors[i], label=f"Amp{i+1}")
		plt.xlim([-10,1.1*np.max(energyList[i])])
		
	plt.xlabel("True Energy [keV]")
	if traceInteg:
		plt.ylabel(f"Fit Mean [V*ns] (from integral)")
	else:
		plt.ylabel(f"Fit Mean [V] (from peak height)")	
	plt.legend(loc="best")
	plt.grid()
	saveto=f"{dataDir}{label}_comparePixels.pdf"
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)
	
def muEnResPlot(dir, ele, pix, calibPix):
	en, mu, sig, enRes, fits, muErr, sigErr, enresErr = enResFitting.getCalibVals_fromTxt(dir, ele, integral=traceInteg)
	fig, ax1 = plt.subplots()
	color = 'tab:blue'
	plt.xlabel('Fit function')
	plt.ylabel('Calibrated mean [keV]')
	plt.errorbar(fits,mu,yerr=muErr, fmt='o',color=color)
	plt.axhline(y=en, color='black', linestyle='-')
	ax1.tick_params(axis='y', labelcolor=color)
	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:red'
	ax2.set_ylabel('Energy resolution [%]')  # we already handled the x-label with ax1
	plt.errorbar(fits,enRes, yerr=enresErr, fmt='o',color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	plt.title(f"{ele} {en} keV line - amp{pix}")
	saveto=dir+ele+"_amp"+str(pix)+"_muEnRes_"+label+"_calib"+str(calibPix)+".pdf"
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)


def chi2RatioPlot(dir, ele, pix, calibPixel):
	chi2, ndof = getChi2(dir)
	en, mu, sig, enRes, fits, muErr, sigErr, enresErr = enResFitting.getCalibVals_fromTxt(dir, ele, integral=traceInteg)
	
	enRatio = [x/en for x in mu]
	for i,item in enumerate(ndof):
		if item==0:
			ndof[i]=1e-8
	goodness = [i / j for i, j in zip(chi2, ndof)]
	
	fig, ax1 = plt.subplots()
	color = 'tab:blue'
	plt.xlabel('Fit function')
	plt.ylabel('Calibrated / true energy')
	plt.plot(fits,enRatio,'o',color=color)
	plt.axhline(y=1., color='black', linestyle='-')
	ax1.tick_params(axis='y', labelcolor=color)
	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:red'
	ax2.set_ylabel('Chi2/ndof')  # we already handled the x-label with ax1
	plt.plot(fits,goodness,'o',color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	plt.title(f"{ele} {en} keV line - amp{pix}")
	saveto=dir+ele+"_amp"+str(pix)+"_chi2Ratio_"+label+"_calib"+str(calibPixel)+".pdf"
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)
	
def sigmaPlot(dir, ele, pix, calibPixel):
	en, mu, sig, enRes, fits, muErr, sigErr, enresErr = enResFitting.getCalibVals_fromTxt(dir, ele, integral=traceInteg)
	
	lim={22.16:2.13, 59.54:3.49, 122.06:5, 88.03:4.25} #sigma values necessary to achieve proposal sensitivity
	
	fig, ax1 = plt.subplots()
	color = 'tab:blue'
	plt.xlabel('Fit function')
	plt.ylabel('Sigma [keV]')
	plt.errorbar(fits,sig,yerr=sigErr, fmt='o',color=color)
	plt.axhline(y=lim[en], color='black', linestyle='-')
	ax1.tick_params(axis='y', labelcolor=color)
	plt.title(f"{ele} {en} keV line - amp{pix}")
	saveto=dir+ele+"_amp"+str(pix)+"_sigma_"+label+"_calib"+str(calibPixel)+".pdf"
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)

def muEnResScan(dir, fit, pix, calibPixel, ratioBool=False):
	plt.rcParams.update({'font.size': 18})
	en, mu, sig, enRes, muErr, sigErr, enresErr = enResFitting.getCalibVals_fromFit(dir, fit, integral=traceInteg)
	if ratioBool:
		ratio=[mu[i]/en[i] for i in range(len(mu))]
		#fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [3,1]})
		fig, axs = plt.subplots(2, sharex=True, figsize=(7,7), gridspec_kw={'height_ratios': [3,1]})

	#plot data
	pl=plt.scatter(en,mu, c=enRes, s=100, cmap='gist_rainbow')
	model = np.polyfit(en, mu, 1)
	linear_model=np.poly1d(model)
	x=np.arange(0,max(en))
	if ratioBool:
		plt.cla()
		#replot data on correct axes
		axs[0].scatter(en,mu, c=enRes, s=100, cmap='gist_rainbow')
		axs[0].plot(x,linear_model(x),'k--')
		cbar = fig.colorbar(pl, ax=axs.ravel().tolist())
		axs[0].plot([], [], ' ', label=f"Fit: y={model[0]:.2f}x+{model[1]:.2f}")
		axs[0].legend(loc='best')
		plt.setp(axs[0],ylabel=f"Calibrated Energy [keV]")
		#add ratio plot
		axs[1].plot(en,ratio,'bo')
		axs[1].hlines(1,1,max(en)*1.1,colors='black')
		plt.setp(axs[1],ylabel=f"Ratio \n(calibrated/true)")
		saveto=dir+fit+"_amp"+str(pix)+"_muEnResScan_"+label+"_calib"+str(calibPixel)+"_ratio.pdf"
	else:
		plt.plot(x,linear_model(x),'k--')
		cbar = plt.colorbar(pl)
		plt.plot([], [], ' ', label=f"Fit: y={model[0]:.2f}x+{model[1]:.2f}")
		plt.legend(loc='best')
		plt.ylabel("Calibrated Energy [keV]")
		plt.tight_layout() #reduce margin space	
		saveto=dir+fit+"_amp"+str(pix)+"_muEnResScan_"+label+"_calib"+str(calibPixel)+".pdf"
	plt.xlabel('True Energy [keV]')
	plt.suptitle(f"{fit} calibration with amp{calibPixel}, Amp{pix} data")
	cbar.set_label('Energy resolution [%]') 
	plt.savefig(saveto) if savePlots else saveFromInput(saveto)
	plt.clf()

############################################################################################
# Main
############################################################################################

####################
#raw data
if v1:
	dataDir=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp1_peaks/"
	dataDir2=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp2_peaks/"	
	energyList2, muArr2, sigmaArr2, enResArr2, muErr2, sigErr2, enresErr2 = enResFitting.getVals_fromTxt(dataDir2+"fitSpectra/", integral=traceInteg)	

else:
	dataDir=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v2_chip{chip}/amp1_peaks/"
	#only one pixel for v2
	
energyList, muArr, sigmaArr, enResArr, muErr, sigErr, enresErr = enResFitting.getVals_fromTxt(dataDir+"fitSpectra/", integral=traceInteg)

	
#short circuit function calls with boolean operators
#if first condition (boolean) is False, then second condition (function evaluation) not evaluated because overall condition is False
if v1:
	(allPlots or rawEnergyResolution) and rawEnRes(dataDir,energyList,enResArr,energyList2=energyList2, enResArr2=enResArr2) 
	(allPlots or rawMuSigma) and rawMuSig(dataDir,sigmaArr,sigErr,muArr,muErr,sigmaArr2=sigmaArr2,sigErr2=sigErr2,muArr2=muArr2,muErr2=muErr2)
	(allPlots or calibrationCurve) and calibCurve(dataDir,energyList,muArr,muErr,energyList2=energyList2,muArr2=muArr2,muErr2=muErr2)
else:
	(allPlots or rawEnergyResolution) and rawEnRes(dataDir,energyList,enResArr) 
	(allPlots or rawMuSigma) and rawMuSig(dataDir,sigmaArr,sigErr,muArr,muErr)
	(allPlots or calibrationCurve) and calibCurve(dataDir,energyList,muArr,muErr)



####################
#calibrated data
if v1:
	dataPixel=[1,2]
	calibPixel=[1,2]
	dataDirCalib=[f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp{dataPixel[0]}_peaks/amp{calibPixel[0]}Calib/",
		f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp{dataPixel[1]}_peaks/amp{calibPixel[1]}Calib/"]
	sources=["Cad","Am","Cobalt57-calib_122.06"] if chip==3 else ["Cad","Cobalt57-calib_122.06"]
	calibFits=["tri","spline3"]
else:
	#only one pixel for v2
	dataPixel=[1]
	calibPixel=[1]
	dataDirCalib=[f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v2_chip{chip}/amp{dataPixel[0]}_peaks/amp{calibPixel[0]}Calib/"]
	sources=["Cad","Cobalt57-calib_122.06"]
	calibFits=["tri"]


for source in sources:
	for i,p in enumerate(dataPixel):
		(allPlots or calibratedMuEnRes) and muEnResPlot(dataDirCalib[p-1],source,p,calibPixel[i])
		(allPlots or calibratedMuChi2) and chi2RatioPlot(dataDirCalib[p-1],source,p,calibPixel[i])
		(allPlots or calibratedSigma) and sigmaPlot(dataDirCalib[p-1],source,p,calibPixel[i])
		
for calib in calibFits:
	for i,p in enumerate(dataPixel):
			(allPlots or calibratedMuEnResScan) and muEnResScan(dataDirCalib[p-1],calib,p,calibPixel[i],ratioBool=True)
