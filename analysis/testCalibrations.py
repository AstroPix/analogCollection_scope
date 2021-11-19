import matplotlib.pyplot as plt
import numpy as np
import h5py
import enResFitting

def flattenIn(data):
	out=[]
	for evt in data:
		flatData=np.array(evt).flatten('F')
		out.append(flatData[0])
	return out
	
def muEnResPlot(dir,ele, pix, savePlots):
	plt.clf()
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
	saveto=dir+"muEnRes_testFits.pdf"
	plt.savefig(saveto) if savePlots else plt.show()
	plt.clf()

	
savePlots=False
dataDir1="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp1_peaks/"
dataDir2="/Users/asteinhe/AstroPixData/astropixOut_tmp/energyCalibration/amp2_peaks/"

	
#raw data
energyList1, muArr1, sigmaArr1, nArr1, enResArr1 = enResFitting.getVals_fromTxt(dataDir1+"fitSpectra/")
errArr1= enResFitting.calcError(sigmaArr1, nArr1)

energyList2, muArr2, sigmaArr2, nArr2, enResArr2 = enResFitting.getVals_fromTxt(dataDir2+"fitSpectra/")
errArr2= enResFitting.calcError(sigmaArr2, nArr2)


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
plt.savefig(f"peaks_comparePixels.pdf") if savePlots else plt.show()
plt.clf()


#calibrated DataLossWarning
muEnResPlot(dataDir1, "Cad", 1, savePlots)
muEnResPlot(dataDir2, "Cad", 2, savePlots)
muEnResPlot(dataDir1, "Am", 1, savePlots)
muEnResPlot(dataDir2, "Am", 2, savePlots)

	
