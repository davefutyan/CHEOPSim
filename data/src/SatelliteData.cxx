/*
 * SatelliteData.cxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

#include "SatelliteData.hxx"

using namespace std;

constexpr double SatelliteData::kDefaultCcdTemperature;
constexpr double SatelliteData::kDefaultFeeBiasTemperature;
constexpr double SatelliteData::kDefaultFeeAdcTemperature;
constexpr double SatelliteData::kCcdTemperatureSigma[2];
constexpr double SatelliteData::kFeeBiasTemperatureSigma[2];
constexpr double SatelliteData::kFeeAdcTemperatureSigma[2];
constexpr double SatelliteData::kVoltageSigma[2][4];
constexpr double SatelliteData::kVoltageDrift[2][4];

SatelliteData::SatelliteData() {
	m_ccdTemperature = kDefaultCcdTemperature;
	m_feeBiasTemperature = kDefaultFeeBiasTemperature;
	m_telescopeTemperature = kDefaultTelescopeTemperature;
	m_fluctuatedCcdTemperature = kDefaultCcdTemperature;
	m_fluctuatedFeeBiasTemperature = kDefaultFeeBiasTemperature;
	m_fluctuatedFeeAdcTemperature = kDefaultFeeAdcTemperature;
	for (unsigned i=0; i<4; i++) m_fluctuatedVoltage[i] = 0.;
	m_dpuTemperature = kDefaultDpuTemperature;
	m_moonAngle = 0.;
	m_sunAngle = 0.;
	m_earthLimbAngle = 0.;
	m_earthOccultationFlag = false;
	m_SAAFlag = false;
	m_strayLightFlux = 0.;
	m_strayLightFlag = false;
	m_discardFlag = false;
}

void SatelliteData::addJitterAPE(APE ape, bool validAocs, bool validScience) {
	m_jitterAPEs.push_back(ape);
	m_validAocsFlags.push_back(validAocs);
	m_validScienceFlags.push_back(validScience);
}

bool SatelliteData::validAocs() const {
	bool validAocs = true;
	for (unsigned i=0; i<m_validAocsFlags.size(); i++) {
		if (!m_validAocsFlags[i]) validAocs = false;
	}
	return validAocs;
}

bool SatelliteData::validScience() const {
	bool validScience = true;
	for (unsigned i=0; i<m_validScienceFlags.size(); i++) {
		if (!m_validScienceFlags[i]) validScience = false;
	}
	return validScience;
}
