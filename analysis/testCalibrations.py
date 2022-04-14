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

savePlots=False
traceInteg=False
v1=True
chip=3 #Can be: 3, 4 [for v1], 1, 2 [for v2]
label = "integ" if traceInteg else "peaks"

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
	
def saveFromInput():
	inp= 'a'
	while inp!='y' and inp!='n':
		print("Save plot? y/n")
		inp = input()
	
	result=True if inp=='y' else False
	return result
   
def tosave(saveFile): 
	plot=plt.gcf() #get current figure - saves fig in case savePlt==True
	plt.show() #creates new figure for display
	savePlt=saveFromInput()
	if savePlt: #save plt stored with plt.gcf()
		print(f"Saving {saveFile}")
		plot.savefig(saveFile)
	
###############################
#Plotting functions - callable
###############################

def rawEnRes(dataDir,energyList, enResArr,energyList2=[],enResArr2=[]):
	enresmu=flattenIn(enResArr)
	plt.plot(energyList, enresmu, 'o', label="Amp1")
	plt.xlabel("True Energy [keV]")
	plt.ylabel(f"Energy Resolution [%]")
	if len(enResArr2)>0:
		enres2mu=flattenIn(enResArr2)
		plt.plot(energyList2, enres2mu, 'or', label="Amp2")	
		plt.legend(loc="best")
	plt.grid()
	#plt.savefig(f"{dataDir}{label}_comparePixels_enRes.pdf") if savePlots else plt.show()
	saveto=f"{dataDir}{label}_comparePixels_enRes.pdf"
	plt.savefig(saveto) if savePlots else tosave(saveto)
	plt.clf()
	
def muEnResPlot(dir, ele, pix, savePlots):
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
	saveto=dir+ele+"_amp"+str(pix)+"_muEnRes_"+label+"_calib2.pdf"
	plt.savefig(saveto) if savePlots else plt.show()

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

def chi2RatioPlot(dir, ele, pix, savePlots):
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
	saveto=dir+ele+"_amp"+str(pix)+"_chi2Ratio_"+label+"_calib.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	
def sigmaPlot(dir, ele, pix, savePlots):
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
	saveto=dir+ele+"_amp"+str(pix)+"_sigma_"+label+"_calib.pdf"
	plt.savefig(saveto) if savePlots else plt.show()

def muEnResScan(dir, fit, pix, savePlots, ratioBool=False):
	plt.rcParams.update({'font.size': 18})
	en, mu, sig, enRes, muErr, sigErr, enresErr = enResFitting.getCalibVals_fromFit(dir, fit, integral=traceInteg)
	if ratioBool:
		ratio=[mu[i]/en[i] for i in range(len(mu))]
		fig, axs = plt.subplots(2, sharex=True, figsize=(9,4), gridspec_kw={'height_ratios': [3,1]})
	pl=plt.scatter(en,mu, c=enRes, s=100, cmap='gist_rainbow')
	model = np.polyfit(en, mu, 1)
	linear_model=np.poly1d(model)
	x=np.arange(0,max(en))
	plt.plot(x,linear_model(x),'k--')
	plt.xlabel('True Energy [keV]')
	plt.ylabel("Calibrated Energy [keV]")
	#if ratioBool:
	#	cbar = plt.colorbar(pl), ax=axs.ravel().tolist()
	#else:
	cbar = plt.colorbar(pl)
	cbar.set_label('Energy resolution [%]') 
	print(len(model))
	plt.plot([], [], ' ', label=f"Fit: y={model[0]:.2f}x+{model[1]:.2f}")
	plt.legend(loc='best')
	plt.tight_layout()
	if ratioBool:
		axs[1].plot(en,ratio,'bo')
		axs[1].hlines(1,1,max(en)*1.1,colors='black')
		plt.setp(axs[1],ylabel=f"Ratio (calibrated/true)")
	saveto=dir+fit+"_amp"+str(pix)+"_muEnResScan_"+label+"_calib.pdf"
	plt.savefig(saveto) if savePlots else plt.show()


############################################################################################
# Main
############################################################################################
if v1:
	dataDir=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp1_peaks/"
	dataDir2=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp2_peaks/"
else:
	dataDir=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v2_chip{chip}/amp1_peaks/"
	#only one pixel for v2
	
#raw data
energyList, muArr, sigmaArr, enResArr, muErr, sigErr, enresErr = enResFitting.getVals_fromTxt(dataDir+"fitSpectra/", integral=traceInteg)

if v1:
	energyList2, muArr2, sigmaArr2, enResArr2, muErr2, sigErr2, enresErr2 = enResFitting.getVals_fromTxt(dataDir2+"fitSpectra/", integral=traceInteg)
	

#plot energy resolution 
if v1:
	rawEnRes(dataDir,energyList,enResArr,energyList2=energyList2, enResArr2=enResArr2)
else:
	rawEnRes(dataDir,energyList,enResArr)


"""
#plot components of energy calibration - mu and sigma
#AMANDA - calculate chi2 and display on plot
sig1=flattenIn(sigmaArr1)
mu1=flattenIn(muArr1)
muE1=flattenIn(muErr1)
sigE1=flattenIn(sigErr1)
x = np.linspace(np.min(mu1), np.max(mu1), 1000)
plt.errorbar(mu1,sig1, xerr=muE1, yerr=sigE1, fmt='o', label="Amp1")
popt, pcov = curve_fit(sqrtFit, mu1, sig1, absolute_sigma=True) #square root
sqrt_fn = sqrtFit(x, *popt)
chi2=calcChi2(sig1,mu1,sigE1,popt)
plt.plot(x,sqrt_fn,'--b',label=f"{popt[0]:.3f}/sqrt(E)+{popt[1]:.3f}\n chi2/ndof={chi2:.4f}")
#popt, pcov = curve_fit(linFit, mu1, sig1, absolute_sigma=True) #linear
#lin_fn = linFit(x, *popt)
#chi2=calcChi2(sig1,mu1,sigE1,popt)
#plt.plot(x,lin_fn,':b',label=f"{popt[0]:.3f}E+{popt[1]:.3f}\n chi2/ndof={chi2:.4f}")
if v1:
	sig2=flattenIn(sigmaArr2)
	mu2=flattenIn(muArr2)
	muE2=flattenIn(muErr2)
	sigE2=flattenIn(sigErr2)
	plt.errorbar(mu2,sig2, xerr=muE2, yerr=sigE2, fmt='o', color="red", label="Amp2")
	popt2, pcov2 = curve_fit(sqrtFit, mu2, sig2, absolute_sigma=True) #square root
	sqrt_fn2 = sqrtFit(x, *popt2)
	chi22=calcChi2(sig2,mu2,sigE2,popt2)
	plt.plot(x,sqrt_fn2,'--r',label=f"{popt2[0]:.3f}/sqrt(E)+{popt2[1]:.3f}\n chi2/ndof={chi22:.4f}")
	#popt2, pcov2 = curve_fit(linFit, mu2, sig2, absolute_sigma=True) #linear
	#lin_fn2 = linFit(x, *popt2)
	#chi22=calcChi2(sig2,mu2,sigE2,popt2)
	#plt.plot(x,lin_fn2,':r',label=f"{popt2[0]:.3f}E+{popt2[1]:.3f}\n chi2/ndof={chi22:.4f}")
	plt.legend(loc="best")
plt.xlabel("Measured energy [V]")
plt.ylabel(f"Measured sigma [V]")	
plt.grid()
plt.savefig(f"{dataDir1}{label}_comparePixels_muVsSig.pdf") if savePlots else plt.show()
plt.clf()


#plot calibration curves
plt.errorbar(energyList1, mu1, yerr=muE1, fmt='o', label="Amp1")
plt.xlabel("True Energy [keV]")
if traceInteg:
	plt.ylabel(f"Fit Mean [V*ns] (from integral)")
else:
	plt.ylabel(f"Fit Mean [V] (from peak height)")	
if v1:
	plt.errorbar(energyList2, mu2, yerr=muE2, fmt='o', color="red", label="Amp2")
	plt.legend(loc="best")
plt.grid()
plt.xlim([-10,1.1*np.max(energyList1)])
plt.savefig(f"{dataDir1}{label}_comparePixels.pdf") if savePlots else plt.show()
plt.clf()


#calibrated data
if v1:
	dataDir1=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp1_peaks/amp1Calib/"
	dataDir2=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v1_chip00{chip}/amp2_peaks/amp2Calib/"
else:
	dataDir1=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v2_chip{chip}/amp1_peaks/amp1Calib/"
	dataDir2=f"/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/v2_chip{chip}/amp1_peaks/amp1Calib/"
	#only one pixel for v2

muEnResPlot(dataDir1, "Cad", 1, savePlots)
muEnResPlot(dataDir2, "Cad", 2, savePlots)
if chip==3:
	muEnResPlot(dataDir1, "Am", 1, savePlots)
	muEnResPlot(dataDir2, "Am", 2, savePlots)
muEnResPlot(dataDir1, "Cobalt57-calib_122.06", 1, savePlots)
muEnResPlot(dataDir2, "Cobalt57-calib_122.06", 2, savePlots)

chi2RatioPlot(dataDir1, "Cad", 1, savePlots)
chi2RatioPlot(dataDir2, "Cad", 2, savePlots)
if chip==3:
	chi2RatioPlot(dataDir1, "Am", 1, savePlots)
	chi2RatioPlot(dataDir2, "Am", 2, savePlots)
chi2RatioPlot(dataDir1, "Cobalt57-calib_122.06", 1, savePlots)
chi2RatioPlot(dataDir2, "Cobalt57-calib_122.06", 2, savePlots)

sigmaPlot(dataDir1, "Cad", 1, savePlots)
sigmaPlot(dataDir2, "Cad", 2, savePlots)
if chip==3:
	sigmaPlot(dataDir1, "Am", 1, savePlots)
	sigmaPlot(dataDir2, "Am", 2, savePlots)
sigmaPlot(dataDir1, "Cobalt57-calib_122.06", 1, savePlots)
sigmaPlot(dataDir2, "Cobalt57-calib_122.06", 2, savePlots)


muEnResScan(dataDir1,"tri",1,savePlots)
if v1:
	muEnResScan(dataDir1,"spline3",1,savePlots)
	muEnResScan(dataDir2,"spline3",2,savePlots)
"""