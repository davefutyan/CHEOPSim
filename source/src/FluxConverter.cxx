/*
 * FluxConverter.cxx
 *
 *  Created on: 2 Jun 2020
 *      Author: futyand
 */

#include "data/include/SatelliteData.hxx"
#include "StarProducer.hxx"

#include "FluxConverter.hxx"

FluxConverter::FluxConverter(double blackBodyOffset, double fitIntercept, double fitSlope) : m_blackBodyOffset(blackBodyOffset), m_fitIntercept(fitIntercept), m_fitSlope(fitSlope) {}

double FluxConverter::cheopsMinusRefBandMagnitude(WavelengthDependence::REFERENCE_BAND refBand, double effectiveTemperature, const WavelengthDependence * wavelengthDependence, bool target) {

	//Implement the formulae described in https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion

	double integratedTransmission_refBand_Vega = wavelengthDependence->getIntegratedRefBandTransmission(refBand,0.,target);
	double integratedThroughputQE_CHEOPS_Vega = wavelengthDependence->getIntegratedThroughputQE(0.,SatelliteData::kDefaultCcdTemperature,target);
	double F0_G = StarProducer::kNPhoVega * integratedTransmission_refBand_Vega;
	double F0_X = StarProducer::kNPhoVega * integratedThroughputQE_CHEOPS_Vega;
	//cout << "cheopsMinusRefBandMagnitude: " << integratedThroughput_refBand_Vega << " " << integratedThroughput_CHEOPS_Vega << " " << F0_G << " " << F0_X << endl;

	double integratedTransmission_refBand = wavelengthDependence->getIntegratedRefBandTransmission(refBand,double(effectiveTemperature),target);
	double integratedThroughputQE_CHEOPS = wavelengthDependence->getIntegratedThroughputQE(double(effectiveTemperature),SatelliteData::kDefaultCcdTemperature,target);

	double F_G = StarProducer::kNPhoVega * integratedTransmission_refBand;
	double F_X = StarProducer::kNPhoVega * integratedThroughputQE_CHEOPS;
	double deltaMag = -2.5 * log10((F_X * F0_G) / (F0_X * F_G));
	//cout << "deltaMag(X-V) for TEff = " << effectiveTemperature << ": " << deltaMag << endl;

	//Decrease deltaMag by blackBodyOffset for Teff>7200K in order to ensure a smooth transition between the use of SEDs (Teff<7200K) and the use of black bodies (Teff>7200K) to model the stellar spectra
	if (effectiveTemperature>7200.) deltaMag += m_blackBodyOffset;

	//Scale the fit parameters so that the fit value is 1 at the effective temperature of Vega. See https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion for details.
	double fitScaleFactor = m_fitIntercept + m_fitSlope * 9600;
	double scaledFitIntercept = m_fitIntercept/fitScaleFactor;
	double scaledFitSlope = m_fitSlope/fitScaleFactor;

	//Correct for the empirically observed flux deficit. See https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion for details.
	deltaMag -= 2.5 * log10(scaledFitIntercept + scaledFitSlope * effectiveTemperature);

	return deltaMag;

}
