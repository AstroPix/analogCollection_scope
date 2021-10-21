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

def get_data(run_time, data_set_name, no_of_traces = 100, noise_range = (0, 2000), signal_range = (2000,10000)):

    # Determines number of points from each scope trace and creates an empty array
    _, last_trace = scope.read_triggered_event()
    last_trace = np.array([int(s) for s in last_trace.split(',')])
    
    curve_length = len(last_trace)
    #arr = np.zeros((no_of_traces, curve_length))
    arr = [] 
    
    # Scaling dictionary is used to scale the scope traces to account for scope settings
    scaling_dict = scope.read_scaling_config()
    
    peaks = []
    integrals = []
    
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
                time_scaled, trace_scaled = scope.scale_data(scaling_dict, trace)
                
                if (no_of_traces < 0) or (i < no_of_traces):
                    arr.append(trace_scaled)
                    time_axis = time_scaled
                #if i % 1000 == 0:
                #    print("At {0:d} / {1:d}".format(i, no_of_traces))
                
                noise_sample = np.mean(trace_scaled[noise_range[0]:noise_range[1]])
                trace_scaled -= noise_sample
                peaks.append( np.max(trace_scaled)  )
                integrals.append( np.sum(trace_scaled[signal_range[0]:signal_range[1]] ) )    
                
        
        # Code to override Visa errors:
        except visa.VisaIOError:
            n_bad_comm += 1
            print("Communication timeout... %d" %(n_bad_comm))
            i -= 1 
        
        i += 1       
    
        #if i > 0 and i %1000 == 0 and single:
        #    t_now = time.time()
        #    single = False 
        #    elapsed = time.time() - t_start
        #    rate = float(i)/float(t_now - t_start)
        #    print("{2:s} ; At {0:d}/{1:d}".format(i, no_of_traces, time.strftime('%a, %d %b %Y %H:%M:%S GMT', time.localtime())))
        #    print("\tRate: {0:6.3f} Hz\t Elapsed: {1:6.2f} s\t Estimated total run length: {2:6.2f} s\t Estimated time remaining: {3:6.2f}".format(rate, elapsed, no_of_traces/(rate), no_of_traces/(rate) - elapsed))

        #if i%1000 == 1:
        #    single = True
    
    t_stop = time.time()
    
    run_len = t_stop - t_start
    run_min = run_len / 60
    print(f"Recorded {i} traces in {run_min:0.3f} minutes. Average rate: {i/run_len:.2f} Hz")

    
    dset = f.create_dataset(data_set_name, data=np.array(arr))    
    dset = f.create_dataset(data_set_name+"_t", data=np.array(time_axis))   

    dset = f.create_dataset(data_set_name+"_peaks", data=np.array(peaks))   
    dset = f.create_dataset(data_set_name+"_integral", data=np.array(integrals))   

    return


