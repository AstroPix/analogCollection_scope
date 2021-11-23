import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys,os,glob
sys.path.insert(1, 'energyCalibration/')
import enResFitting

def flattenIn(data):
	out=[]
	for evt in data:
		flatData=np.array(evt).flatten('F')
		out.append(flatData[0])
	return out
	
def muEnResPlot(dir,ele, pix, savePlots):
	en, mu, sig, n, enRes, fits = enResFitting.getCalibVals_fromTxt(dir, ele)
	fig, ax1 = plt.subplots()
	color = 'tab:blue'
	plt.xlabel('Fit function')
	plt.ylabel('Calibrated mean [keV]')
	plt.plot(fits,mu,'o',color=color)
	plt.axhline(y=en, color='black', linestyle='-')
	ax1.tick_params(axis='y', labelcolor=color)
	ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
	color = 'tab:red'
	ax2.set_ylabel('Energy resolution [%]')  # we already handled the x-label with ax1
	plt.plot(fits,enRes,'o',color=color)
	ax2.tick_params(axis='y', labelcolor=color)
	plt.title(f"{ele} {en} keV line - amp{pix}")
	saveto=dir+ele+"_amp"+str(pix)+"_muEnRes_testFits.pdf"
	plt.savefig(saveto) if savePlots else plt.show()

def getChi2(indir):
	chiArr, ndofArr = [],[]
	fits=['linear','quad','tri','sqrt','spline1','spline3', 'piecewise']

	for fit in fits:
		os.chdir(indir+fit+'/')
		fitFile = glob.glob('fit*.txt') #returns array with length 1
		openFile = open(fitFile[0],'r')
		lines=openFile.readlines()
		fitline=lines[1].split(' = ')
		vals=fitline[1].split('/')
		chiArr.append(float(vals[0]))
		ndofArr.append(float(vals[1]))
		
	return chiArr, ndofArr

def chi2RatioPlot(dir, ele, pix, savePlots):

	chi2, ndof = getChi2(dir)
	en, mu, sig, n, enRes, fits = enResFitting.getCalibVals_fromTxt(dir, ele)
	
	enRatio = [x/en for x in mu]
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
	saveto=dir+ele+"_amp"+str(pix)+"_chi2Ratio_testFits.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()



############################################################################################
############################################################################################
savePlots=True
dataDir1="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp1_peaks/"
dataDir2="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp2_peaks/"

	
#raw data
energyList1, muArr1, sigmaArr1, nArr1, enResArr1 = enResFitting.getVals_fromTxt(dataDir1+"fitSpectra/")
errArr1= enResFitting.calcError(sigmaArr1, nArr1)

energyList2, muArr2, sigmaArr2, nArr2, enResArr2 = enResFitting.getVals_fromTxt(dataDir2+"fitSpectra/")
errArr2= enResFitting.calcError(sigmaArr2, nArr2)

#plot energy resolution 
enres1mu=flattenIn(enResArr1)
plt.plot(energyList1, enres1mu, 'o', label="Amp1")
enres2mu=flattenIn(enResArr2)
plt.plot(energyList2, enres2mu, 'or', label="Amp2")
plt.xlabel("True Energy [keV]")
plt.ylabel(f"Energy Resolution [%]")	
plt.legend(loc="best")
plt.grid()
plt.savefig(f"{dataDir1}peaks_comparePixels_enRes.pdf") if savePlots else plt.show()
plt.clf()

#plot calibration curves
amp1mu=flattenIn(muArr1)
err1=flattenIn(errArr1)
plt.errorbar(energyList1, amp1mu, yerr=err1, fmt='o', label="Amp1")
amp2mu=flattenIn(muArr2)
err2=flattenIn(errArr2)
plt.errorbar(energyList2, amp2mu, yerr=err2, fmt='o', color="red", label="Amp2")
plt.xlabel("True Energy [keV]")
plt.ylabel(f"Fit Mean [V] (from peak height)")	
plt.legend(loc="best")
plt.grid()
plt.xlim([-10,1.1*np.max(energyList1)])
plt.savefig(f"{dataDir1}peaks_comparePixels.pdf") if savePlots else plt.show()
plt.clf()


#calibrated data
#dataDir1="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp1_peaks_noiseIncl/"
#dataDir2="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp2_peaks_noiseIncl/"

muEnResPlot(dataDir1, "Cad", 1, savePlots)
muEnResPlot(dataDir2, "Cad", 2, savePlots)
muEnResPlot(dataDir1, "Am", 1, savePlots)
muEnResPlot(dataDir2, "Am", 2, savePlots)

chi2RatioPlot(dataDir1, "Cad", 1, savePlots)
chi2RatioPlot(dataDir2, "Cad", 2, savePlots)
chi2RatioPlot(dataDir1, "Am", 1, savePlots)
chi2RatioPlot(dataDir2, "Am", 2, savePlots)
