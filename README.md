# TRDNoisePadCalib
Calibration of the TRD noise and pad status in LHC Run3

This project is a part of my CERN Service Task. 
This includes the calibration of the TRD noise and pad status by regularly taking and using the TRD standalone NOISE RUNS. 
##variance.C : This macro basically loops over all the digits and all the 1.2M channels and then stores the desired digit and channel information. Thereby this info is used to calculate the adcSum, adcSumSquared, adcMean, variance and adcRMS values. 
##checkNoise.C : This macro goes through all the 1.2M channels per PadRow per Layer per Stack per Sector and then plotting them shows in general where the noisy pads are actually located.
##pulseheightplot.C : As the name suggest, the macro plots the pulse height plots yet as a means of diagonizing the NOISY Channels.
##noiseandpad.C : This macro uses the adcMean and adcRMS information along with the global Channel and global Pad numbers in order to obtain the NOISE MAP per sector per stack per layer.
