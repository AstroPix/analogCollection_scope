import matplotlib.pyplot as plt
import h5py


homeDir = "/Users/asteinhe/AstroPixData/astropixOut_tmp"



#display histogram to play with
def histDisplay(f_in,ds='run1_peaks'):
	f=h5py.File(f_in, 'r')

	data=f[ds]
	xBinWidth=0.002
	xMax=np.max(data)
	xMin=np.min(data)
	binEdges=np.arange(xMin,xMax+xBinWidth,xBinWidth)#use peakMax+xBinWidth to overshoot range and include all data

	#Create histogram of data
	plt.hist(data,bins=binEdges,label=r'Data', color='blue')
	plt.show()


#combine separate runs into one new file
def combineFiles(files,outFile):

	datap=np.array([])
	datai=np.array([])
	datas=np.array([])
	for fin in files:
		f=h5py.File(homeDir+fin,'r')
		datainp=np.array(f['run1_peaks'])
		dataini=np.array(f['run1_integral'])
		datains=np.array(f['run1_scaling'])
		datap=np.concatenate([datap,datainp])
		datai=np.concatenate([datai,dataini])
		datas=np.concatenate([datas,datains])
		f.close()

	h=h5py.File(outFile, 'w')
	h.create_dataset('run1_peaks', data=datap)
	h.create_dataset('run1_integral', data=datai)
	h.create_dataset('run1_scaling', data=datas)
	print(list(h.keys()))
	h.close()
	


#copy scaling dataset into main dataset
def copyScalingDS(f_in_scale, f_in):

	f=h5py.File(f_in_scale, 'r')

	my_array=f['run1_scaling']

	g = h5py.File(f_in, 'a')
	g.create_dataset('run1_scaling', data=my_array)
	print(list(g.keys()))
	f.close()
	g.close()



#################################################################
# main
#################################################################

if __name__ == "__main__":
	
	f="110421_amp1/Americium_480min_combined.h5py"
	filesIn=["102921_amp1/americium241_90min.h5py", "110421_amp1/Americium_120min.h5py","110421_amp1/test_Americium_30min.h5py","110421_amp1/_Americium_240min.h5py"]
	outFile="110421_amp1/Americium_480min_combined.h5py"
	
	histDisplay(homeDir+outFile)
	combineFiles(filesIn, outFile)
	copyScalingDS("102021_amp1/cobalt57_14h.h5py","102021_amp2/cobalt57_14h_scaling.h5py")

		
		
		
		
		


