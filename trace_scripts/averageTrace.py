import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy.stats
import pandas as pd

def get_average_trace( filename, dataset ):

    f = h5py.File(filename)
        

    traces = f[dataset]
#            time = np.array(f[k+"_t"])

    print( traces.shape )
    return np.mean(traces, axis = 0)
    
   

theDir = "/Users/hfleisc1/astroPix/astropixOut_tmp/"

otherF = h5py.File( theDir + "101821_amp1/0.05Vinj_HV25.h5py" )

time = np.array(otherF["run1_t"])


#for inj in np.linspace( 0.3, 1.8, 16 ):
for inj in [0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8 ] :

    av = get_average_trace( theDir + f"101521_amp1/{inj:.1f}Vinj.h5py", "testRun1" )
    
    baseline = np.mean( av[0:1900] )
    peak = np.max(av)

    #plt.plot( time*1e6, (av-baseline)*1e3, label = f"{inj:.1f} V" )
    #plt.plot( time*1e6, (av-baseline)/(peak - baseline), label = f"{inj:.1f} V" )
    plt.plot( time*1e6, (av-baseline)/inj, label = f"{inj:.1f} V" )


plt.legend(loc = "best")
plt.xlabel( "time (μs)" )
#plt.ylabel( "Normalized signal")
plt.ylabel( "Baseline-subtraced signal/injection voltage")
plt.title("Signal shape vs injection voltage")

plt.savefig( "shape_full_ratio.pdf")
plt.xlim( 0, 5)
plt.savefig( "shape_rise_ratio.pdf")
plt.xlim( -125, -80)
#plt.ylim( -0.01, 0.05)
plt.ylim( -0.001, 0.006)
#plt.ylim( -1, 12)
plt.savefig( "shape_prepulse_ratio.pdf" )
plt.show()
