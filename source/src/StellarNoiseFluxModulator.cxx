/*
 * StellarNoiseFluxProducer.cxx
 *
 *  Created on: 3 Feb 2014
 *      Author: futyand
 */

#include <iostream>
#include <fstream>

#include "StellarNoiseFluxModulator.hxx"

void StellarNoiseFluxModulator::doBegin(Data* data, bool fullFrame) {

	if (data->getFieldOfView()->getStars().size() == 0) {
		throw runtime_error("Error in StellarNoiseFluxModulator::doBegin: No stars in FOV. Check that the pointing direction corresponds to the target star coordinates.");
	}
	Star::spectralType spectralType = data->getFieldOfView()->getStars()[m_starIndex]->getSpectralType();

	string granulationFile;
	switch (spectralType) {
	case Star::F0:
		granulationFile = "granulation_F0_M1.4_R1.8_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F1:
		granulationFile = "granulation_F1_M1.4_R1.6_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F2:
		granulationFile = "granulation_F2_M1.4_R1.6_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F3:
		granulationFile = "granulation_F3_M1.4_R1.6_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F4:
		granulationFile = "granulation_F4_M1.4_R1.5_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F5:
		granulationFile = "granulation_F5_M1.4_R1.5_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F6:
		granulationFile = "granulation_F6_M1.2_R1.4_Teff6500_Tstep15_fft.dat";
		break;
	case Star::F7:
		granulationFile = "granulation_F7_M1.2_R1.2_Teff6000_Tstep15_fft.dat";
		break;
	case Star::F8:
		granulationFile = "granulation_F8_M1.2_R1.2_Teff6000_Tstep15_fft.dat";
		break;
	case Star::F9:
		granulationFile = "granulation_F9_M1.2_R1.2_Teff6000_Tstep15_fft.dat";
		break;
	case Star::G0:
		granulationFile = "granulation_G0_M1.2_R1.1_Teff6000_Tstep15_fft.dat";
		break;
	case Star::G1:
		granulationFile = "granulation_G1_M1.2_R1.1_Teff6000_Tstep15_fft.dat";
		break;
	case Star::G2:
		granulationFile = "granulation_G2_M1.1_R1.0_Teff6000_Tstep15_fft.dat";
		break;
	case Star::G3:
		granulationFile = "granulation_G3_M1.0_R1.0_Teff5500_Tstep15_fft.dat";
		break;
	case Star::G4:
		granulationFile = "granulation_G4_M1.0_R1.0_Teff5500_Tstep15_fft.dat";
		break;
	case Star::G5:
		granulationFile = "granulation_G5_M1.0_R1.0_Teff5500_Tstep15_fft.dat";
		break;
	case Star::G6:
		granulationFile = "granulation_G6_M1.0_R1.0_Teff5500_Tstep15_fft.dat";
		break;
	case Star::G7:
		granulationFile = "granulation_G7_M1.0_R1.0_Teff5500_Tstep15_fft.dat";
		break;
	case Star::G8:
		granulationFile = "granulation_G8_M0.9_R0.9_Teff5500_Tstep15_fft.dat";
		break;
	case Star::G9:
		granulationFile = "granulation_G9_M0.9_R0.9_Teff5500_Tstep15_fft.dat";
		break;
	case Star::K0:
		granulationFile = "granulation_K0_M0.9_R0.8_Teff5500_Tstep15_fft.dat";
		break;
	case Star::K1:
		granulationFile = "granulation_K1_M0.8_R0.8_Teff5000_Tstep15_fft.dat";
		break;
	case Star::K2:
		granulationFile = "granulation_K2_M0.8_R0.8_Teff5000_Tstep15_fft.dat";
		break;
	case Star::K3:
		granulationFile = "granulation_K3_M0.8_R0.7_Teff5000_Tstep15_fft.dat";
		break;
	case Star::K4:
		granulationFile = "granulation_K4_M0.7_R0.7_Teff4500_Tstep15_fft.dat";
		break;
	case Star::K5:
		granulationFile = "granulation_K5_M0.7_R0.7_Teff4500_Tstep15_fft.dat";
		break;
	case Star::K6:
		granulationFile = "granulation_K6_M0.6_R0.6_Teff4000_Tstep15_fft.dat";
		break;
	case Star::K7:
		granulationFile = "granulation_K7_M0.6_R0.6_Teff4000_Tstep15_fft.dat";
		break;
	case Star::K8:
		granulationFile = "granulation_K8_M0.6_R0.6_Teff4000_Tstep15_fft.dat";
		break;
	case Star::K9:
		granulationFile = "granulation_K9_M0.6_R0.5_Teff4000_Tstep15_fft.dat";
		break;
	case Star::M0:
		granulationFile = "granulation_M0_M0.6_R0.5_Teff4000_Tstep15_fft.dat";
		break;
	case Star::M1:
		granulationFile = "granulation_M1_M0.5_R0.5_Teff3500_Tstep15_fft.dat";
		break;
	case Star::M2:
		granulationFile = "granulation_M2_M0.5_R0.5_Teff3500_Tstep15_fft.dat";
		break;
	case Star::M3:
		granulationFile = "granulation_M3_M0.4_R0.4_Teff3500_Tstep15_fft.dat";
		break;
	case Star::M4:
		granulationFile = "granulation_M4_M0.2_R0.3_Teff3000_Tstep15_fft.dat";
		break;
	case Star::M5:
		granulationFile = "granulation_M5_M0.2_R0.3_Teff3000_Tstep15_fft.dat";
		break;
	case Star::M6:
		granulationFile = "granulation_M6_M0.2_R0.3_Teff3000_Tstep15_fft.dat";
		break;
	case Star::M7:
		granulationFile = "granulation_M7_M0.2_R0.3_Teff3000_Tstep15_fft.dat";
		break;
	case Star::M8:
		granulationFile = "granulation_M8_M0.2_R0.3_Teff3000_Tstep15_fft.dat";
		break;
	case Star::M9:
		granulationFile = "granulation_M9_M0.2_R0.3_Teff3000_Tstep15_fft.dat";
		break;
	default:
		throw runtime_error("ERROR: Invalid spectral type in StellarNoiseFluxModulator::doBegin. Supported spectral types are [F,G,K,M][0-9], e.g. K5.");
	}

	cout << "granulation file: " << granulationFile << endl;

	ifstream granulation_in(string(getenv("CHEOPS_SW"))+"/resources/"+granulationFile, ios::in);
	if (!granulation_in.good()) throw runtime_error("Error opening incident flux granulation file: "+granulationFile);
	int i=0;
	string dummy;
	while(granulation_in.good()) {
		granulation_in >> dummy;
		if(granulation_in.eof()) break;
		granulation_in >> m_granulation[i];
		//cout << "granulation[" << i << "]=" << m_granulation[i] << endl;
		i++;
	}
	granulation_in.close();

}

void StellarNoiseFluxModulator::process(Data* data, int timeStep, bool fullFrame) const {

	unsigned currentTime = static_cast<unsigned>(lround(data->getTimeConfiguration().getTimeSinceStart(timeStep)));
    double fluxFactor = 1.+1E-6*m_granulation[(currentTime/15)%11518]; //input file contains 11518 time steps at 15 second intervals
	data->getFieldOfView()->getStars()[m_starIndex]->getTimeSeries()[timeStep]->setNoiseFluxFactor(fluxFactor);

}
