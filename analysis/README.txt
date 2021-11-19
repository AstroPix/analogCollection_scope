RUNNING ENERGY CALIBRATION SCRIPT

1) Populate `amp1_inputs.txt` with paths to data files from amp1. Use FULL path!
2) Populate `amp2_inputs.txt` with paths to data files from amp2. Use FULL path!
3) Check the run settings in `runOptions.txt` and amend as necessary
4) Run the code with :
	python energyCalibration.py runOptions.txt
5) If `savePlots=True` in runOptions.txt, output plots are saved in `saveDir`




AMANDA
- read in amp1_inputs
- read in amp2_ inputs
- design runOptions
- read in runOptions
- separate useful scripts into own dir
- clean up enResFitting
- class to store info rather than arrays 