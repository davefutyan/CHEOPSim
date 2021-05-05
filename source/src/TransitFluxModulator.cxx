/*
 * TransitFluxModulator.cxx
 *
 *  Created on: Dec 18, 2013
 *      Author: futyand
 */

#include <cmath>
#include <iostream>
#include <fstream>

#include "boost/filesystem.hpp"

#include "BarycentricOffset.hxx"

#include "TransitFluxModulator.hxx"

using namespace std;

void TransitFluxModulator::initialize(const ModuleParams & params) {

	m_firstTransitTimeFraction = params.GetAsDouble("firstTransitTime");
	m_planetRadius = params.GetAsDouble("planetRadius");
	string planetScale = params.GetAsString("planetScale");
	m_orbitPeriod = params.GetAsDouble("orbitPeriod");
	string orbitUnits = params.GetAsString("orbitUnits");
	m_impactParameter = params.GetAsDouble("impactParameter");
	m_doLimbDarkening = params.GetAsBool("doLimbDarkening");
	m_doModification = params.GetAsBool("doModification");

	if (orbitUnits == "hours") {
		m_orbitPeriod *= 3600.;
	} else if (orbitUnits == "days") {
		m_orbitPeriod *= (3600.*24.);
	} else {
		throw runtime_error("Error in TransitFluxModulator::initialize: orbitUnits must be hours or days");
	}
	if (planetScale == "Jupiter") {
		m_planetRadius *= kJupiterRadius;
	} else if (planetScale == "Neptune") {
		m_planetRadius *= kNeptuneRadius;
	} else if (planetScale == "Earth") {
		m_planetRadius *= kEarthRadius;
	} else {
		throw runtime_error("Error in TransitFluxModulator::initialize: planetScale must be Jupiter or Earth");
	}

}

void TransitFluxModulator::doBegin(Data * data, bool fullFrame) {

	if (data->getFieldOfView()->getStars().size() == 0) {
		throw runtime_error("Error in TransitFluxModulator::doBegin: No stars in FOV. Check that the pointing direction corresponds to the target star coordinates.");
	}
	Star * star = data->getFieldOfView()->getStars()[m_starIndex];

	m_starRadius = star->getRadius()*kSunRadius;
	if (m_starRadius<1.E-6) throw runtime_error("Error in TransitFluxModulator::initializeTransitModel: star radius is zero or negative");
	double radiusRatio = m_planetRadius/m_starRadius;
	pair<double,double> limbDarkeningCoeffs = star->getLimbDarkeningCoefficients();

	m_transitModel = new TransitModel(radiusRatio,m_doLimbDarkening,limbDarkeningCoeffs);

	double mu_star = star->getMass()*kMuSun;
	if (mu_star<1.E-6) throw runtime_error("Error in TransitFluxModulator::initializeTransitModel: star mass is zero or negative");

	m_semiMajorAxis = pow(mu_star/pow(2.*M_PI/m_orbitPeriod,2),1./3.);

	if (m_semiMajorAxis < m_starRadius) {
		throw runtime_error("Error in TransitFluxModulator configuration for star "+to_string(m_starIndex)+": semi major axis of orbit is smaller than star radius");
	}

	double inclination = acos(m_impactParameter*m_starRadius/m_semiMajorAxis);
	double transitDuration = m_orbitPeriod/M_PI*asin(m_starRadius*sqrt(pow((1.+radiusRatio),2)-pow(m_impactParameter,2))/(m_semiMajorAxis*sin(inclination)));

	TimeConfiguration timeConf = data->getTimeConfiguration();
	DeltaTime startTimeToTransit = m_firstTransitTimeFraction * timeConf.getDuration();
	UTC transitTime = timeConf.getStartTimeUTC() + startTimeToTransit;
	BarycentricOffset offset(star->getSkyPosition().getRightAscension()/3600.,star->getSkyPosition().getDeclination()/3600.);
	m_transitTime_BJD = transitTime.getMjd() + offset;

	cout << "Transit parameters:" << endl;
	cout << "   transitDuration = " << transitDuration/3600. << endl;
	cout << "   m_semiMajorAxis = " << m_semiMajorAxis << endl;
	cout << "   inclination = " << inclination << endl;
	cout << "   radiusRatio = " << radiusRatio << endl;
	cout << "   limb darkening coefficients: " << limbDarkeningCoeffs.first << ", " << limbDarkeningCoeffs.second << endl;
	cout << "   transit depth = " << (1.0 - m_transitModel->fluxFactor(m_impactParameter)) << endl;

	// Modify the analytic transit curve to add structure defined in an input file
	if (m_doModification) {
		unsigned modification[148][30];
		ifstream modification_in(string(getenv("CHEOPS_SW"))+"/resources/transitModification.txt", ios::in);
		if (!modification_in.good()) throw runtime_error("Error opening incident flux modification file");
		int j=0;
		int i=0;
		while(!modification_in.eof()) {
			modification_in >> modification[i][j];
			//cout << "modification[" << i << "][" << j << "]=" << modification[i][j] << endl;
			i++;
			if (i==148) {
				i=0;
				j++;
			}
		}
		modification_in.close();

		for (int i=0; i<148; i++) {
			for (int j=0; j<30; j++) {
				if (modification[i][j]==1) {
					m_transitModification[i] = (((15.-j)/15.)/5000.)+1.;
					continue;
				}
			}
		}
	}

}

void TransitFluxModulator::process(Data * data, int timeStep, bool fullFrame) const {

	TimeConfiguration timeConf = data->getTimeConfiguration();
    double currentTime = timeConf.getTimeSinceStart(timeStep);
    double firstTransitTime = m_firstTransitTimeFraction * timeConf.getDuration();

    int nOrbit = floor(double( (currentTime-firstTransitTime)/m_orbitPeriod ));
	double transitTime = firstTransitTime + nOrbit*m_orbitPeriod;
	double orbitFraction = (currentTime-transitTime)/m_orbitPeriod;

	double starRadiusFraction = std::sqrt(pow(m_semiMajorAxis/m_starRadius*sin(2*M_PI*orbitFraction),2) +
			                              pow(m_impactParameter*cos(2*M_PI*orbitFraction),2));

	double fluxFactor=1.;
	if (starRadiusFraction < (m_starRadius+m_planetRadius/m_starRadius) && fabs(orbitFraction-0.5)>0.25 ) {
		fluxFactor *= m_transitModel->fluxFactor(starRadiusFraction);
	}

	if (m_doModification) {
		int nStackedImages = static_cast<int>(timeConf.getNumberOfStackedImages());
		int exposuresPerStack = static_cast<int>(timeConf.getExposuresPerStack());
		if (timeStep>(nStackedImages/2-148/2)*exposuresPerStack && timeStep<(nStackedImages/2+148/2)*exposuresPerStack) {
			fluxFactor *= m_transitModification[timeStep/exposuresPerStack-(nStackedImages/2-148/2)];
		}
	}

	//cout << "starRadiusFraction = " << starRadiusFraction << ", fluxFactor = " << fluxFactor << endl;
	Star * star = data->getFieldOfView()->getStars()[m_starIndex];
	star->getTimeSeries()[timeStep]->setTransitFluxFactor(fluxFactor);

	//Although not time-step dependent, set the transit time and transit period data members for the star here
	//rather than in doBegin because in the case of full frame on, star information gets reset in second call
	//to StarProducer::doBegin (1st call for full frame, 2nd call for sub-array)
	//whereas TransitFluxModulator::doBegin is only called for the full frame case
	star->setTransitTime(m_transitTime_BJD);
	star->setTransitPeriod(m_orbitPeriod/(3600.*24.));

//	ofstream transitFluxFile("transitTimeseries.txt", ios::app);
//	if (timeStep%10==0) transitFluxFile << currentTime << "	" << fluxFactor << endl;
//	transitFluxFile.close();

}
