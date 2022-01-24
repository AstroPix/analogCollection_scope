import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy
from scipy.optimize import curve_fit
import sys,os,glob
sys.path.insert(1, 'energyCalibration/')
import enResFitting


savePlots=True
traceInteg=True
label = "integ" if traceInteg else "peaks"


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
		sum_square += (num/denom)
		i+=1
	ndof = np.float64(len(true)-1)#loose 1 dof 
	res_var = sum_square/ndof	
	return res_var
	
	
	
	
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
	saveto=dir+ele+"_amp"+str(pix)+"_chi2Ratio_"+label+"_calib2.pdf"
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
	saveto=dir+ele+"_amp"+str(pix)+"_sigma_"+label+"_calib2.pdf"
	plt.savefig(saveto) if savePlots else plt.show()



############################################################################################
############################################################################################
dataDir1="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/chip003/amp1_peaks/"
dataDir2="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/chip003/amp2_peaks/"

	
#raw data
energyList1, muArr1, sigmaArr1, enResArr1, muErr1, sigErr1, enresErr1 = enResFitting.getVals_fromTxt(dataDir1+"fitSpectra/", integral=traceInteg)

energyList2, muArr2, sigmaArr2, enResArr2, muErr2, sigErr2, enresErr2 = enResFitting.getVals_fromTxt(dataDir2+"fitSpectra/", integral=traceInteg)


#plot energy resolution 
enres1mu=flattenIn(enResArr1)
plt.plot(energyList1, enres1mu, 'o', label="Amp1")
enres2mu=flattenIn(enResArr2)
plt.plot(energyList2, enres2mu, 'or', label="Amp2")
plt.xlabel("True Energy [keV]")
plt.ylabel(f"Energy Resolution [%]")	
plt.legend(loc="best")
plt.grid()
plt.savefig(f"{dataDir1}{label}_comparePixels_enRes.pdf") if savePlots else plt.show()
plt.clf()

#plot components of energy calibration - mu and sigma

#AMANDA - calculate chi2 and display on plt
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
"""
popt, pcov = curve_fit(linFit, mu1, sig1, absolute_sigma=True) #linear
lin_fn = linFit(x, *popt)
chi2=calcChi2(sig1,mu1,sigE1,popt)
plt.plot(x,lin_fn,':b',label=f"{popt[0]:.3f}E+{popt[1]:.3f}\n chi2/ndof={chi2:.4f}")
"""
sig2=flattenIn(sigmaArr2)
mu2=flattenIn(muArr2)
muE2=flattenIn(muErr2)
sigE2=flattenIn(sigErr2)
plt.errorbar(mu2,sig2, xerr=muE2, yerr=sigE2, fmt='o', color="red", label="Amp2")
popt2, pcov2 = curve_fit(sqrtFit, mu2, sig2, absolute_sigma=True) #square root
sqrt_fn2 = sqrtFit(x, *popt2)
chi22=calcChi2(sig2,mu2,sigE2,popt2)
plt.plot(x,sqrt_fn2,'--r',label=f"{popt2[0]:.3f}/sqrt(E)+{popt2[1]:.3f}\n chi2/ndof={chi22:.4f}")
"""
popt2, pcov2 = curve_fit(linFit, mu2, sig2, absolute_sigma=True) #linear
lin_fn2 = linFit(x, *popt2)
chi22=calcChi2(sig2,mu2,sigE2,popt2)
plt.plot(x,lin_fn2,':r',label=f"{popt2[0]:.3f}E+{popt2[1]:.3f}\n chi2/ndof={chi22:.4f}")
"""
plt.xlabel("Measured energy [V]")
plt.ylabel(f"Measured sigma [V]")	
plt.legend(loc="best")
plt.grid()
plt.savefig(f"{dataDir1}{label}_comparePixels_muVsSig.pdf") if savePlots else plt.show()
plt.clf()


#plot calibration curves
plt.errorbar(energyList1, mu1, yerr=muE1, fmt='o', label="Amp1")
plt.errorbar(energyList2, mu2, yerr=muE2, fmt='o', color="red", label="Amp2")
plt.xlabel("True Energy [keV]")
if traceInteg:
	plt.ylabel(f"Fit Mean [V*ns] (from integral)")
else:
	plt.ylabel(f"Fit Mean [V] (from peak height)")	
plt.legend(loc="best")
plt.grid()
plt.xlim([-10,1.1*np.max(energyList1)])
plt.savefig(f"{dataDir1}{label}_comparePixels.pdf") if savePlots else plt.show()
plt.clf()


#calibrated data
dataDir1="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/chip003/amp1_peaks/amp2Calib/"
dataDir2="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/chip003/amp2_peaks/amp2Calib/"

muEnResPlot(dataDir1, "Cad", 1, savePlots)
muEnResPlot(dataDir2, "Cad", 2, savePlots)
#muEnResPlot(dataDir1, "Am", 1, savePlots)
#muEnResPlot(dataDir2, "Am", 2, savePlots)
muEnResPlot(dataDir1, "Cobalt57-calib_122.06", 1, savePlots)
muEnResPlot(dataDir2, "Cobalt57-calib_122.06", 2, savePlots)

chi2RatioPlot(dataDir1, "Cad", 1, savePlots)
chi2RatioPlot(dataDir2, "Cad", 2, savePlots)
#chi2RatioPlot(dataDir1, "Am", 1, savePlots)
#chi2RatioPlot(dataDir2, "Am", 2, savePlots)
chi2RatioPlot(dataDir1, "Cobalt57-calib_122.06", 1, savePlots)
chi2RatioPlot(dataDir2, "Cobalt57-calib_122.06", 2, savePlots)

sigmaPlot(dataDir1, "Cad", 1, savePlots)
sigmaPlot(dataDir2, "Cad", 2, savePlots)
#sigmaPlot(dataDir1, "Am", 1, savePlots)
#sigmaPlot(dataDir2, "Am", 2, savePlots)
sigmaPlot(dataDir1, "Cobalt57-calib_122.06", 1, savePlots)
sigmaPlot(dataDir2, "Cobalt57-calib_122.06", 2, savePlots)