import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy.stats
import pandas as pd


theDir = "/Users/hfleisc1/astroPix/astropixOut_tmp/"
#f = h5py.File( theDir + "102021_amp1/cadmium109_scalingTest.h5py" )
f = h5py.File( theDir + "102021_amp1/cadmium109_45min.h5py" )
print( f.keys() )


time = np.array(f["run1_t"])
traces = np.array(f["run1"])

peaks = np.array(f["run1_peaks"])
integral = np.array(f["run1_integral"])

print (peaks.shape)

exit()

plt.plot( peaks*1e3, integral, ".", markersize=3, alpha = 0.2 )
plt.xlabel( "Peak voltage [mV]" )
plt.ylabel( "Integrated voltage [V]" )

plt.savefig("Peak_vs_integral.pdf")



exit()

print( traces.shape )

for i in np.arange(0, 10):

    for j in np.arange(i*10, (i+1)*10):
    
        plt.plot( time*1e6, traces[j]*1e3, label = j)
    
    plt.xlabel( "time (μs)" )
    plt.legend(loc = "best")
    plt.ylabel( "signal (mV) ")
    plt.savefig(f"Cadmium2_{i}.png")
    plt.clf()
    
