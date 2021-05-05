/*
 * TruthData2.cxx
 *
 *  Created on: 4 Dec 2014
 *      Author: futyand
 */

#include <memory>
#include <iostream>

#include "TruthData.hxx"

TruthData::~TruthData() {
	m_psfs.clear();
	m_hotPixels.clear();
	m_deadPixels.clear();
	m_cosmicPixels.clear();
	m_smearTrailRow.clear();
	m_globalThroughput.clear();
	m_gain.clear();
}

TruthData& TruthData::operator +=(const TruthData& truthData) {

	//Combine the saturation information
	this->m_fullWellSaturated = (this->m_fullWellSaturated || truthData.m_fullWellSaturated);
	this->m_adcSaturated = (this->m_adcSaturated || truthData.m_adcSaturated);

	//Add the global throughput and gain values to the list (to be averaged in the case of image stacking)
	this->m_globalThroughput.insert(this->m_globalThroughput.end(),truthData.m_globalThroughput.begin(),truthData.m_globalThroughput.end());
	this->m_gain.insert(this->m_gain.end(),truthData.m_gain.begin(),truthData.m_gain.end());

	//Increment the zodiacal light and stray light
	this->m_zodiacalLight += truthData.m_zodiacalLight;
	this->m_strayLight += truthData.m_strayLight;

	//Combine the PSF information
	if (this->m_psfs.size()==0 && truthData.m_psfs.size()>0) {//for image stacking: stackedImages is initially an empty vector
		this->m_psfs.insert(this->m_psfs.end(),truthData.m_psfs.begin(),truthData.m_psfs.end());
	} else if (this->m_psfs.size() == truthData.m_psfs.size()) {
		for (unsigned i=0; i<this->m_psfs.size();  i++) {
			this->m_psfs[i] += truthData.m_psfs[i];
		}
	} else {
		throw runtime_error("TruthData::operator+=: mismatch in number of PSFs");
	}

	//Increment the hot pixel fluxes
	if (this->m_hotPixels.size()==0 && truthData.m_hotPixels.size()>0) {
		this->m_hotPixels.insert(this->m_hotPixels.end(),truthData.m_hotPixels.begin(),truthData.m_hotPixels.end());
	} else if (this->m_hotPixels.size() == truthData.m_hotPixels.size()) {
		for (unsigned i=0; i<this->m_hotPixels.size();  i++) {
			this->m_hotPixels[i].m_nElectrons += truthData.m_hotPixels[i].m_nElectrons;
			//If the hot pixel being added is a telegraphic pixel in the active state, set the combined state to active
			if (truthData.m_hotPixels[i].m_type == telegraphicActive) this->m_hotPixels[i].m_type = telegraphicActive;
		}
	} else {
		throw runtime_error("TruthData::operator+=: mismatch in number of hot pixels");
	}

	//Propagate the dead pixel information
	if (this->m_deadPixels.size()==0 && truthData.m_deadPixels.size()>0) {
		this->m_deadPixels.insert(this->m_deadPixels.end(),truthData.m_deadPixels.begin(),truthData.m_deadPixels.end());
	} else if (this->m_deadPixels.size() != truthData.m_deadPixels.size()) {
		throw runtime_error("TruthData::operator+=: mismatch in number of dead pixels");
	}

	//Concatenate the cosmic pixel vectors of the images being added
	this->m_cosmicPixels.insert(this->m_cosmicPixels.end(),
								truthData.m_cosmicPixels.begin(),
								truthData.m_cosmicPixels.end());

	//Increment the smear trail row contents
	if (this->m_smearTrailRow.size()==0 && truthData.m_smearTrailRow.size()>0) {
		this->m_smearTrailRow.insert(this->m_smearTrailRow.end(),truthData.m_smearTrailRow.begin(),truthData.m_smearTrailRow.end());
	} else if (this->m_smearTrailRow.size() == truthData.m_smearTrailRow.size()) {
		for (unsigned i=0; i<this->m_smearTrailRow.size();  i++) {
			this->m_smearTrailRow[i] += truthData.m_smearTrailRow[i];
		}
	} else {
		throw runtime_error("TruthData::operator+=: mismatch in smear trail row length");
	}

	return *this;

}

TruthData& TruthData::operator =(const TruthData& truthData) {

	this->m_fullWellSaturated = truthData.m_fullWellSaturated;
	this->m_adcSaturated = truthData.m_adcSaturated;
	this->m_zodiacalLight = truthData.m_zodiacalLight;
	this->m_strayLight = truthData.m_strayLight;

	m_psfs.clear();
	m_hotPixels.clear();
	m_deadPixels.clear();
	m_cosmicPixels.clear();
	m_smearTrailRow.clear();
	m_globalThroughput.clear();
	m_gain.clear();

	this->m_psfs.insert(this->m_psfs.begin(),truthData.m_psfs.begin(),truthData.m_psfs.end());
	this->m_hotPixels.insert(this->m_hotPixels.begin(),truthData.m_hotPixels.begin(),truthData.m_hotPixels.end());
	this->m_deadPixels.insert(this->m_deadPixels.begin(),truthData.m_deadPixels.begin(),truthData.m_deadPixels.end());
	this->m_cosmicPixels.insert(this->m_cosmicPixels.begin(),truthData.m_cosmicPixels.begin(),truthData.m_cosmicPixels.end());
	this->m_smearTrailRow.insert(this->m_smearTrailRow.begin(),truthData.m_smearTrailRow.begin(),truthData.m_smearTrailRow.end());
	this->m_globalThroughput.insert(this->m_globalThroughput.begin(),truthData.m_globalThroughput.begin(),truthData.m_globalThroughput.end());
	this->m_gain.insert(this->m_gain.begin(),truthData.m_gain.begin(),truthData.m_gain.end());

	return *this;

}
