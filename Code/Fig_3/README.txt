SinusoidalModulationPrediction.m: Generates the curves shown in Fig. 3 one at a time. Adjust 'wm' and the 'FileName' parameters accordingly for each cure you wish to plot.   

Notes on running the script above: Designed to use a full 128 GB of DDR5 RAM on a 24-core Intel i9-13900KF workstation and will run for an hour or so. Depending on the CPU speed and RAM constraints on your system, you may have to intelligently allocate the number of workers in your parpool, the 'Fineness' (integration step) metric and the sweep points in the 'm' variable.   

Plot_SinModulation.m: Plots the curve data obtained from 'SinusoidalModulationPrediction'. Included are three sets of generated curve data:
wm_0_001.mat : w=0.001
wm_0_01.mat  : w=0.01
wm_0_1.mat   : w=0.1



