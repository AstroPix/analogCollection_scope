import matplotlib.pyplot as plt
import numpy as np
import h5py
import scipy.stats
import pandas as pd

def get_noise_level( filename, p ):

    names = []
    means = []
    sigmas = []

    print(filename)

    f = h5py.File(filename)
    for k in f.keys():
            if k[-2:] in ["_t", "al", "ks"] :
                continue

            print(filename, k)

            traces = f[k]
            time = np.array(f[k+"_t"])


            noise = traces[:,0:5000]
            signal = traces[:,5000:]
            
            bl = np.unique([noise, signal])
            
            delta = np.min( bl[1:] - bl[:-1])
            
            bins = np.concatenate( (bl - delta/2, [bl[-1]+delta/2]) )

            p=p.replace("/", "_")

            plt.plot( time*1e6, np.mean(traces, axis=0) )
            plt.xlabel("time [μs]")
            plt.ylabel("signal [V]")
            plt.grid(True)
            plt.title(f"{p} {k}")
            plt.savefig(f"average_trace_{p}_{k}.png")
            plt.clf()

            mean = np.mean(noise)
            sig = np.std(noise)
            print( mean, sig)
            x = np.linspace(bins[0]-delta, bins[-1]+delta, 100)
            y = scipy.stats.norm(loc=mean, scale=sig).pdf(x)
            #bins=np.linspace(mean-6*sig, mean+8*sig, 100)
            plt.hist( noise.flatten(), bins=bins, histtype="step", label = "Noise" )
            plt.hist( signal.flatten(), bins=bins, histtype="step", label = "Pulse")
            plt.plot(x, 5000*500*delta*y, label = "Gaussian fit (noise)")
            plt.legend(loc="best")
            plt.xlabel("scope signal [V]")
            plt.ylabel("number of occurrences")
            plt.yscale("log")
            plt.grid(True)
            plt.ylim(ymin = 0.01, ymax = None)
            plt.title(f"{p} {k}")
            plt.savefig(f"noise_hist_{p}_{k}.png")
            plt.clf()
            names.append(k)
            means.append(mean)
            sigmas.append(sig)
            
            
            sp = np.fft.rfft(noise)
            freq = np.fft.rfftfreq(noise.shape[-1])
            freq = freq[1:]
            sp = sp[:,1:]
            plt.plot(freq, np.mean(sp.real.T, axis=1), label = "Noise" )
        
            sp = np.fft.rfft(signal)
            sp = sp[:,1:]
            plt.plot(freq, np.mean(sp.real.T, axis=1), label = "Pulse")
            
            plt.plot(freq, 1e-4/freq, "--", label = "1e-4/f")
        
            plt.legend(loc="best")
            plt.yscale("log")
            plt.xscale("log")
            plt.title(f"{p} {k}")
            plt.xlabel("noise frequency???")
            plt.ylabel("coefficient???")

            plt.savefig(f"noise_freq_{p}_{k}.png")
            plt.clf()

           
    return names, means, sigmas
    
HV = [25, 40, 60, 80]
theDir = "/Users/hfleisc1/astroPix/astropixOut_tmp/"
for pixel in [1, 2]:
    for add in ["/", "/configEachRun/"]:
        theMeans = []
        theSigmas = []
        for hv in HV:
    
            try:
                names, means, sigmas = get_noise_level(f"{theDir}/101821_amp{pixel}/{add}/0.05Vinj_HV{hv}.h5py", f"amp{pixel}{add}HV{hv}")
        
                theMeans.append(means[0])
                theSigmas.append(sigmas[0])
            except:
                pass
        
        if len(theMeans) == 0:
            continue

        plt.plot( HV, (np.array(theMeans)-0.537)*1e3, "-o", label="mean-0.537 V")
        plt.plot( HV, np.array(theSigmas)*1e3, "-x", label="sdev")

        plt.xlabel("Bias voltage (V)")
        plt.ylabel("Noise voltage (mV)")
        plt.grid(True)
        add = add.replace("/", "")

        plt.title(f"pixel {pixel} {add}")
        plt.legend(loc="best")
        plt.savefig(f"noise_vs_HV_amp{pixel}_{add}.png")
        plt.clf()


exit()

plt.show()

exit()


exit()




mean = np.mean(signal)
sig = np.std(signal)
print( mean, sig)
y = scipy.stats.norm(loc=mean, scale=sig).pdf(x)
plt.hist( signal.flatten(), bins=bins, histtype="step" )
plt.plot(x, 1500*y)


plt.yscale("log")
plt.show()

exit()

