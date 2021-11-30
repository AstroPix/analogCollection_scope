RUNNING ENERGY CALIBRATION SCRIPT

1) Populate `amp1_inputs.txt`(`amp2_inputs.txt`) with paths to data files from amp1(amp2). If full path not given, assumed to be in `homeDir` (defined in runOptions)
	a) Fits for Gaussian photopeaks
2) Populate `amp1_edge_inputs.txt`(`amp2_edge_inputs.txt`) with paths to data files from amp1(amp2). If full path not given, assumed to be in `homeDir` (defined in runOptions)
	a) Fits for integrated Gaussian Compton edge
3) Check the run settings in `runOptions.txt` and amend as necessary
4) Run the code with :
	python energyCalibration.py runOptions.txt
	a) Terminal output will show what file is being fit, and print out fit parameters
5) If `savePlots=True` in runOptions.txt, output plots are saved in `saveDir`




AMANDA
- distinguish between gaussian and edge fit in code, not in inputs
- class to store info rather than arrays 