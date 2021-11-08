"""
    This function pulls SiPM pulses from the scope and stores them in arrays.
    
    Parameters
    ----------
    
    no_of_traces: float
        However many traces or SiPM pulses you want to record.
    data_set_name: string
        Name of the data set within the file. The naming scheme for data sets is 
        filename_run#, ex. "BurstCube_PostVibe_Cs137_061419_run1."
        
    Returns
    -------
    
    Data sets. Two data sets should be created, the scope scaling information and 
    The scope data is stored in an array (a h5py data set). h5py data sets are nice 
    because you can splice into them. In this case, the size of the array is determined by the 
    number of traces you want from the scope and number of data points the scope collects for 
    each trace; the size of the array is number of traces by amount of scope points. Each row 
    is a scope trace, so plotting/data analysis is done in a for loop that looks at each row 
    one at a time.
    
    """
import numpy as np
import time
from datetime import datetime
import h5py
import pyvisa as visa

def get_data(file, scope, run_time, data_set_name, no_of_traces = 100, noise_range = (0, 2000), signal_range = (2000,10000), overwrite=False, useTraceLimit=False):

    # Determines number of points from each scope trace and creates an empty array
    _, last_trace = scope.read_triggered_event()
    last_trace = np.array([int(s) for s in last_trace.split(',')])
    
    # Scaling dictionary is used to scale the scope traces to account for scope settings
    scaling_dict = scope.read_scaling_config()
    
    arr = []
    peaks = []
    integrals = []
    peakTime = []
    trigTime = []
    baseline = []
    
    n_bad_comm = 0
    i = 0
    t_start = time.time()
    t_end = t_start + run_time * 60 #run time in minutes
    n_dup = 0
    single = True
    
    time_axis = None
    
    bla = time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime(t_end) )
    
    print(f"Starting {run_time} minute run; will be done at {bla}")
    
    while time.time() < t_end:
        # Code to discount duplicates (when the scope gets stuck on a trigger):
        
        try: 
            _, trace = scope.read_triggered_event()
            trace = np.array([int(s) for s in trace.split(',')])
            if np.sum(trace - last_trace) == 0:
                i -= 1
                n_dup +=1
                print("%d duplicates" %(n_dup))
            else:
                last_trace = trace
                ttime=datetime.now() #returns local time
                trigTime.append(ttime.strftime('%d %b %Y %H:%M:%S.%f'))
                time_scaled, trace_scaled = scope.scale_data(scaling_dict, trace)
                
                if (no_of_traces < 0) or (i < no_of_traces):
                    arr.append(trace_scaled)
                    time_axis = time_scaled
                
                noise_sample = np.mean(trace_scaled[noise_range[0]:noise_range[1]])
                baseline.append( noise_sample )
                trace_scaled -= noise_sample
                peaks.append( np.max(trace_scaled)  )
                peakIndex = np.argmax(trace_scaled, axis=0) #approximate location of signal pulse from highest measurement
                peakTime = time_scaled[peakIndex]
                integrals.append( np.sum(trace_scaled[signal_range[0]:signal_range[1]] ) )  
                
        
        # Code to override Visa errors:
        except visa.VisaIOError:
            n_bad_comm += 1
            print("Communication timeout... %d" %(n_bad_comm))
            i -= 1 
        
        i += 1       
        if useTraceLimit and i >= no_of_traces:
            break

    t_stop = time.time()
    
    run_len = t_stop - t_start
    run_min = run_len / 60
    print(f"Recorded {i} traces in {run_min:0.3f} minutes. Average rate: {i/run_len:.2f} Hz")

    if overwrite:
        for name in ["", "_t", "_peaks", "_integral", "_baseline", "_peakTime", "_trigTime"]: 
            if (data_set_name + name) in file:
                del file[data_set_name + name]

    file.create_dataset(data_set_name, data=np.array(arr))    
    file.create_dataset(data_set_name+"_t", data=np.array(time_axis))   
    file.create_dataset(data_set_name+"_peaks", data=np.array(peaks))   
    file.create_dataset(data_set_name+"_integral", data=np.array(integrals))   
    file.create_dataset(data_set_name+"_baseline", data=np.array(baseline)) 
    file.create_dataset(data_set_name+"_peakTime", data=np.array(peakTime))   
    file.create_dataset(data_set_name+"_trigTime", data=trigTime)
    
    
    return 


