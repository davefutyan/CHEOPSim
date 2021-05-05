/*
 * TimeConfiguration.cxx
 *
 *  Created on: 8 May 2018
 *      Author: futyand
 */

#include "TimeConfiguration.hxx"

using namespace std;

TimeConfiguration::TimeConfiguration(boost::posix_time::ptime startTime, unsigned numberOfStackedImages,
		unsigned exposuresPerStack, double exposureTimeAsDouble, unsigned maxImagesPerCube, bool doFullFrame, double repetitionPeriod) :
		m_visitStartTime(startTime),m_numberOfStackedImages(numberOfStackedImages),m_exposuresPerStack(exposuresPerStack),m_exposureTimeAsDouble(exposureTimeAsDouble),
		m_maxImagesPerCube(maxImagesPerCube),m_numberOfStackedImages_PITL(numberOfStackedImages),m_imagetteStackingNumber(1),m_doFullFrame(doFullFrame),m_numberOfPreSubArrayTimeSteps(0) {

	// Delay the start time (defined as the start time of the sub-array image sequence) by one exposure duration (corresponding to the full frame image) + 35s
	// This is to simulate a typical scenario in which a full frame image is taken 30s after the start of the visit,
	// and there is a 5s gap between the end of the full frame exposure and the start of the sub-array image sequence
	// If no there is full frame image, then the delay is fixed at 30s.
	double delay = m_doFullFrame ? m_exposureTimeAsDouble+35. : 30.;
	long seconds = static_cast<long>(lround((floor(delay))));
	long milliseconds = static_cast<long>(lround((delay - floor(delay))*1000.));
	m_startTime = m_visitStartTime + boost::posix_time::seconds(seconds) + boost::posix_time::milliseconds(milliseconds);
	if (repetitionPeriod > 0. ) { // User has specified a non default repetition period in the configuration
		m_repetitionPeriod = repetitionPeriod;
	} else { // Default repetition period based on exposure duration
		if (m_exposureTimeAsDouble < 1.) {
			//For exposure durations <1s, the repetition period is equal to the exposure duration plus 1s
			m_repetitionPeriod = m_exposureTimeAsDouble + 1.;
		} else {
			//For exposure durations >1s, the value of the repetition period is always an integer in order to sync with jitter which has 1s resolution.
			//In reality, for exposure durations >1s, the repetition period should match the exposure duration.
			//However exposure durations longer than 1s are expected to also always be integers.
			m_repetitionPeriod = ceil(exposureTimeAsDouble);
		}
	}
}

void TimeConfiguration::setImagetteStackingNumber() {

	//Check that exposuresPerStack is valid (appears in CHEOPSim_web/LUT_image_stacking.txt)
	unsigned validImageStackingOrder[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,24,26,28,30,33,36,39,40,44};
	if (find(begin(validImageStackingOrder), end(validImageStackingOrder), m_exposuresPerStack) == end(validImageStackingOrder)) {
		throw runtime_error("Error in TimeConfiguration::setImagetteStackingNumber: invalid value for exposuresPerStack: "+to_string(m_exposuresPerStack));
	}

	//Assign imagette stacking order according to the mapping to image stacking order in CHEOPSim_web/LUT_image_stacking.txt
    if (m_exposuresPerStack<16) {
        m_imagetteStackingNumber = 1;
    } else if (m_exposuresPerStack<30) {
        m_imagetteStackingNumber = 2;
    } else if (m_exposuresPerStack<40) {
        m_imagetteStackingNumber = 3;
    } else {
        m_imagetteStackingNumber = 4;
    }

}

unsigned TimeConfiguration::getNumberOfStackedImageCubes() const {
	unsigned numberOfStackedImages = m_numberOfStackedImages_PITL;
	unsigned nCubes = numberOfStackedImages/m_maxImagesPerCube;
	if (numberOfStackedImages%m_maxImagesPerCube != 0) nCubes+=1;
	return nCubes;
}

unsigned TimeConfiguration::getNumberOfUnstackedImageCubes() const {
	unsigned numberOfTimeSteps = getNumberOfTimeStepsPITL();
	unsigned nCubes = numberOfTimeSteps/m_maxImagesPerCube;
	if (numberOfTimeSteps%m_maxImagesPerCube != 0) nCubes+=1;
	return nCubes;
}
