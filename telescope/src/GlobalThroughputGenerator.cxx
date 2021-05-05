/*
 * GlobalThroughputGenerator.cxx
 *
 *  Created on: 3 Feb 2016
 *      Author: futyand
 */

#include <fstream>

#include "GlobalThroughputGenerator.hxx"
#include "PSFGenerator.hxx"
#include "source/include/FluxConverter.hxx"

void GlobalThroughputGenerator::doBegin(Data* data, bool fullFrame) {

	//Set the black body effective temperature of the target star needed for the calculation of the integrated throughput*QE
	m_effectiveTemperature = 5660.; //Default to effective temperature for G5 star if there are no stars in the FOV
	m_spectralType = "G5";
	if (data->getFieldOfView()->getStars().size()>0) {
		m_effectiveTemperature = (data->getFieldOfView()->getStars()[0])->getEffectiveTemperature();
		m_spectralType = (data->getFieldOfView()->getStars()[0])->getSpectralTypeString();
	}

	//The following file is needed for calculating the limb darkening coefficients (see comments in Star.cxx)
//	ofstream globalThroughputFile("throughputQE.dat", ios::app);
//	WavelengthDependence * wavelengthDependence = data->getWavelengthDependence();
//	for (int wavelength = 330; wavelength<=1100; wavelength++) {
//		globalThroughputFile << wavelength*10 << " " << wavelengthDependence->getThroughput(wavelength)*wavelengthDependence->getQE(wavelength) << endl;
//	}
//	globalThroughputFile.close();

}

void GlobalThroughputGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	Image * image = data->getImages().back();

	//Get the integrated throughput*QE weighted by the black body distribution defined by the effective temperature of the target star
	//taking into account the variation of QE with the CCD temperature
	double ccdTemperature = data->getSatelliteData(timeStep)->getCcdTemperature();
	double integratedThroughputQE = data->getWavelengthDependence()->getIntegratedThroughputQE(m_effectiveTemperature,ccdTemperature,true);

	if (timeStep == 0 || ccdTemperature != SatelliteData::kDefaultCcdTemperature) {
		cout << "Integrated throughput*QE (spectral type " << m_spectralType << ", Teff= " << m_effectiveTemperature << ") = " << integratedThroughputQE << endl;
//		if (timeStep == 0) {
//			//This code can be uncommented for use in combination with photonflux.sh to fill a file with values for photon flux in the CHEOPS pass band for each spectral type
//			ofstream globalThroughput_file("photonFlux_aperture30cm_GaiaPassband.txt", ios::app);
//			globalThroughput_file << m_spectralType << " " << m_effectiveTemperature << " " << data->getFieldOfView()->getStars()[0]->getMeanPhotonFlux()*PSFGenerator::kTelescopeArea
//					<< " " << data->getFieldOfView()->getStars()[0]->getMeanPhotonFlux()*PSFGenerator::kTelescopeArea*integratedThroughputQE << (m_effectiveTemperature>7200 ? " BB" : " SED") << endl;
//			globalThroughput_file.close();
//		}
//		if (timeStep == 0) {
//			//This code can be uncommented for use in combination with globalThroughput.sh to fill a file with global throughput values for each spectral type
//			ofstream globalThroughput_file("globalThroughput.txt", ios::app);
//			globalThroughput_file << m_spectralType << " " << m_effectiveTemperature << " " << data->getFieldOfView()->getStars()[0]->getMeanPhotonFlux()*PSFGenerator::kTelescopeArea << " " << integratedThroughputQE << endl;
//			globalThroughput_file.close();
//		}
//		ofstream globalThroughput_file("globalThroughput.txt", ios::app);
//		double nPho_Vega = 4161462975; //Obtained by running StarProducer with m_VegaSpectrum = true in StarProducer.cxx and magnitude 0.035 in mystars.txt
//		globalThroughput_file << data->getFieldOfView()->getStars()[0]->getSpectralTypeString() << "	" << m_effectiveTemperature << "	"
//				<< setprecision(10) << integratedThroughputQE << "	" << nPho_Vega*integratedThroughputQE << "	"
//				<< setprecision(9) << nPho_Vega*integratedThroughputQE*0.5111 << endl;
//		globalThroughput_file.close();
	}

	//Determine the correction to the predicted flux to account for the empirical discrepancy between observed and predicted fluxes
	//See https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion
	FluxConverter * fluxConv = new FluxConverter();
	double fluxCorrection = fluxConv->fluxCorrectionFactor(m_effectiveTemperature);
	delete fluxConv;

	//Rescale the image by the wavelength integrated throughput*QE multiplied by the empirical correction factor
	int ymin = data->doChargeTransferEoL() ? 0. : image->getYOffset(); //Include pixels below the subarray if end of life CTI is to be simulated
	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=ymin; iy<image->getYOffset()+image->getYDim(); iy++) {

			double nPhotons = image->getPixelValue(ix,iy);
			double nElectrons = nPhotons * integratedThroughputQE * fluxCorrection;
			image->setPixelValue(ix,iy,nElectrons);

		}
	}

	//Also rescale the snapshots of the PSF at the start and end of the exposure to be used for frame transfer smearing
	if (data->doFrameTransferSmearing() && data->getImagesToSmearUp().size() == data->getImages().size()) {

		Image * imageToSmearUp = data->getImagesToSmearUp().back();
		Image * imageToSmearDown = data->getImagesToSmearDown().back();

		int xOffset_smear = imageToSmearUp->getXOffset();
		int yOffset_smear = imageToSmearUp->getYOffset();
		int xDim_smear = imageToSmearUp->getXDim();
		int yDim_smear = imageToSmearUp->getYDim();

		for (int ix=xOffset_smear; ix<xOffset_smear+xDim_smear; ix++) {
			for (int iy=yOffset_smear; iy<yOffset_smear+yDim_smear; iy++) {
				imageToSmearUp->setPixelValue(ix,iy,imageToSmearUp->getPixelValue(ix,iy) * integratedThroughputQE);
				imageToSmearDown->setPixelValue(ix,iy,imageToSmearDown->getPixelValue(ix,iy) * integratedThroughputQE);
			}
		}

	}

	//Set the truth information
	image->getTruthData()->setGlobalThroughput(integratedThroughputQE);

}
