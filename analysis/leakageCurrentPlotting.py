import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import sys,os,glob


def makeDF(fileList):
	#read in SEC and READ columns of first csv file
	df=pd.read_csv(fileList[0],usecols=['SEC', 'READ'])
	#append SEC and READ columns into original dataframe
	for i,f in enumerate(fileList):
		if i==0:
			continue
		df_tmp=pd.read_csv(f,usecols=['SEC', 'READ'])
		df=pd.concat([df,df_tmp])
		
	df=df.reset_index()
	#scale resulting df so first timestamp is at t=0
	time_scale=df['SEC'].min()
	df['SEC']=df['SEC']-time_scale
	#convert current to nA
	df['READ']=np.abs(df['READ']*1e9)
	return df
	

################################################################################
dataDir="/Users/asteinhe/AstroPixData/astropixOut_tmp/ps/"

campaign=["april/hv/","june/hv"]
naming=["inBeam","beam"]
dut1=["chip3","601"]
dut2=["chip2","603"]

for i,folder in enumerate(campaign):
	os.chdir(dataDir+folder)	
	allFiles = glob.glob(f"*{naming[i]}*.csv")
	dut1Files, dut2Files = [], []

	for f in allFiles:
		if ("chip2" in f) or (f.split('_')[2]=="603"):
			dut2Files.append(f)
		else:
			dut1Files.append(f)

	df1=makeDF(dut1Files)
	plt.scatter(df1['SEC'],df1['READ'], s=0.5)
	#plt.axhline() #LC in light after irradiation
	plt.xlabel('Time from first run (s)')
	plt.ylabel('|Leakage current| [nA]')
	plt.yscale('log')
	plt.title(f'DUT 1 ({dut1[i]})')
	plt.show()
	plt.clf()

	df2=makeDF(dut2Files)
	plt.scatter(df2['SEC'],df2['READ'], marker='x', s=0.5)
	#plt.axhline() #LC in light after irradiation
	plt.xlabel('Time from first run (s)')
	plt.ylabel('|Leakage current| [nA]')
	plt.title(f'DUT 2 ({dut2[i]})')
	plt.show()