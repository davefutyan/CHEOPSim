#!/bin/bash

#This script can be used to generate a 2D array of time to saturation as a function of stellar magnitude and spectral type
#The configuration file should be runCHEOPSim_saturation.xml
#The call to FullWellSimulator::calculateTimeToSaturation() must be uncommented

for mag in 5.0 5.5 6.0 6.5 7.0 7.5 8.0 8.5 9.0 9.5 10.0 10.5 11.0 11.5 12.0 12.5 13.0
do
    for specType in O7 O8 O9 B0 B1 B2 B3 B4 B5 B6 B7 B8 B9 A0 A1 A2 A3 A4 A5 A6 A7 A8 A9 F0 F1 F2 F3 F4 F5 F6 F7 F8 F9 G0 G1 G2 G3 G4 G5 G6 G7 G8 G9 K0 K1 K2 K3 K4 K5 K6 K7 K8 K9 M0 M1 M2 M3 M4 M5 M6 M7 M8 M9
    do
	echo $mag $specType
	echo "00:00:00.0 000:00:00.0 "$mag" "$specType > mystars.txt
	rm -r CH_PR900001_TG000001
	bin/runCHEOPSim
    done
done
