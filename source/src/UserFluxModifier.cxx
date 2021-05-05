/*
 * UserFluxModifier.cxx
 *
 *  Created on: 16 Dec 2016
 *      Author: futyand
 */

#include <fstream>

#include "Logger.hxx"

#include "UserFluxModifier.hxx"

void UserFluxModifier::initialize(const ModuleParams& params) {
	m_fluxFilename = params.GetAsString("userFluxFilename");
}

void UserFluxModifier::doBegin(Data* data, bool fullFrame) {

	//user uploaded flux time seris file, to be read from execution directory
	ifstream flux_in(m_fluxFilename, ios::in);
	if (!flux_in.good()) throw runtime_error("Error in UserFluxModifier::doBegin: Error opening file "+m_fluxFilename);
	string line;
	while (getline(flux_in, line)) {
		if(line[0] != '#') {
			double fluxFactor, time;
			stringstream(line) >> time >> fluxFactor;
			m_time.push_back(time);
			m_fluxFactor.push_back(fluxFactor);
		}
	}
	flux_in.close();

	unsigned fileDuration = (m_time.size()-1)*(m_time[1]-m_time[0]);
	if (data->getTimeConfiguration().getDuration() > fileDuration) {
	        throw runtime_error("User defined stellar flux time series is shorter than the simulation duration.");
	}

}

void UserFluxModifier::process(Data* data, int timeStep, bool fullFrame) const {

        TimeConfiguration timeConf = data->getTimeConfiguration();
        double timeSinceStart = timeConf.getTimeSinceStart(timeStep) + (timeConf.getStartTimeUTC()-timeConf.getVisitStartTimeUTC()).getSeconds();

	double fluxFactor = 1.;
	for (unsigned i=1; i<m_time.size(); i++) {
		if (m_time[i] >= timeSinceStart) {
			double mu = (timeSinceStart - m_time[i-1])/(m_time[i]-m_time[i-1]);
			fluxFactor = m_fluxFactor[i-1]*(1.-mu)+m_fluxFactor[i]*mu;
			break;
		}
	}

	if (data->getFieldOfView()->getStars().size() == 0) {
		throw runtime_error("Error in UserFluxModifier: No stars in FOV. Check that the pointing direction corresponds to the target star coordinates.");
	}
	data->getFieldOfView()->getStars()[m_starIndex]->getTimeSeries()[timeStep]->setUserFluxFactor(fluxFactor);

}
