/*
 * StrayLightGenerator.cxx
 *
 *  Created on: 20 Feb 2016
 *      Author: futyand
 */

#include <fstream>

#include "REF_APP_StrayLight.hxx"

#include "StrayLightGenerator.hxx"

void StrayLightGenerator::initialize(const ModuleParams& params) {

	m_strayLightFilename = params.GetAsString("strayLightFile");

	cout << "Reading stray light time series from " << m_strayLightFilename << endl;

	m_visitConstraints = nullptr;

	m_strayLightFromVisitConstraints = m_strayLightFilename.find("MPS_PRE_Visits") != string::npos;

	if (m_strayLightFromVisitConstraints) return;

	if (m_strayLightFilename.substr(m_strayLightFilename.find_last_of(".")+1) == "txt") {

		//user uploaded file, to be read from execution directory
		ifstream sl_in(m_strayLightFilename, ios::in);
	    if (!sl_in.good()) throw runtime_error("Error in StrayLightGenerator::initialize: Error opening file "+m_strayLightFilename);
		string line;
		while (getline(sl_in, line)) {
			if(line[0] != '#') {
				double time, flux;
				stringstream(line) >> time >> flux;
	    		setStrayLightValues(time,flux);
			}
		}
		sl_in.close();

	} else {

    	//Open the stray light fits file
    	string filename = string(getenv("CHEOPS_SW"))+"/resources/"+m_strayLightFilename;
    	RefAppStraylight * sl_file = new RefAppStraylight(filename,"READONLY");

    	//Read in the stray light data
    	while (sl_file->ReadRow()) {
    		setStrayLightValues(sl_file->getCellTime(),sl_file->getCellFlux());
    	}
    	delete sl_file;

    }

}

void StrayLightGenerator::setStrayLightValues(double time, double flux) {

	m_timeFromFile.push_back(time*60); //convert to seconds
	m_strayLightFromFile.push_back(flux);
	cout << "stray light at time " << m_timeFromFile.back() << "s = " << m_strayLightFromFile.back() << " photons/s" << endl;

}

void StrayLightGenerator::doBegin(Data* data, bool fullFrame) {

	data->setPreSubArrayTimeSteps();

	//Set the integrated throughput weighted by the black body distribution defined by the effective temperature of the target star
	double effectiveTemperature = 5660.; //Default to effective temperature for G5 star if there are no stars in the FOV
	if (data->getFieldOfView()->getStars().size()>0) {
		effectiveTemperature = (data->getFieldOfView()->getStars()[0])->getEffectiveTemperature();
	}
	m_integratedThroughput = data->getWavelengthDependence()->getIntegratedThroughput(effectiveTemperature);

	TimeConfiguration timeConf = data->getTimeConfiguration();

	if (m_strayLightFromVisitConstraints) {

		//Open the MPS_PRE_VisitConstraints extension of the MPS_PRE_Visits fits file
		m_visitConstraints = new MpsPreVisitconstraints(m_strayLightFilename + "[MPS_PRE_VisitConstraints]");
		if (m_visitConstraints == nullptr) {
			throw runtime_error("Error in StrayLightGenerator::doBegin: Error opening MPS_PRE_VisitConstraints extension of the MPS_PRE_Visits fits file");
		}
		while (m_visitConstraints->ReadRow()) {
			double timeSinceStart = (m_visitConstraints->getCellUtcTime() - timeConf.getVisitStartTimeUTC()).getSeconds();
			if (timeSinceStart > -65. && timeSinceStart < timeConf.getDuration()+65.) { // VisitConstraints cadence is expected to be 60s
				if (m_visitConstraints->isNullStrayLight()) {
					throw runtime_error("Error in StrayLightGenerator::doBegin: MPS_PRE_VisitConstraints file contains NULL values for stray light.");
				}
				setStrayLightValues(timeSinceStart/60.,m_visitConstraints->getCellStrayLight());
			}
		}
		delete m_visitConstraints;
		if (m_timeFromFile.size() == 0) throw runtime_error("Error in StrayLightGenerator::doBegin: No entries found in MPS_PRE_VisitConstraints file for which the time corresponds to the visit");

	} else {
	  
		if (m_timeFromFile.back() != OrbitSimulator::kCheopsOrbitPeriod && m_timeFromFile.back() < timeConf.getDuration()+60.) {
			throw runtime_error("error in StrayLightGenerator::doBegin: Time of last entry in stray light file must either be equal to the CHEOPSim orbit period ("
					    +to_string(OrbitSimulator::kCheopsOrbitPeriod/60.)+" minutes), or must be larger than the duration of the simulation ("
					    +to_string((timeConf.getDuration()+60.)/60.)+" minutes). Time of last entry is "+to_string(m_timeFromFile.back()/60.)+" minutes.");
		}

	}


	//Set the stray light flag for each time step of the simulation
	for (int timeStep=-int(timeConf.getNumberOfPreSubArrayTimeSteps()); timeStep<int(timeConf.getNumberOfTimeSteps());  timeStep++) {

		SatelliteData * satData = data->getSatelliteData(timeStep,!data->moduleIsActive("JitterProducer"));
		double strayLightFlux = getStrayLightFlux(timeStep, data->getTimeConfiguration());
		//cout << strayLightFlux << " " <<  data->getVisit().m_strayLightThreshold << " " << (strayLightFlux > data->getVisit().m_strayLightThreshold) << endl;
		satData->setStrayLightFlux(strayLightFlux);
		satData->setStrayLightFlag(strayLightFlux > data->getVisit().m_strayLightThreshold);

	}

}

double StrayLightGenerator::getStrayLightFlux(int timeStep, TimeConfiguration timeConf) const {

	unsigned currentTime = timeConf.getTimeSinceStart(timeStep) + (timeConf.getStartTimeUTC()-timeConf.getVisitStartTimeUTC()).getSeconds();
	if (m_timeFromFile.back() == OrbitSimulator::kCheopsOrbitPeriod) currentTime = currentTime%OrbitSimulator::kCheopsOrbitPeriod;

	int ifile = 1;

	//Identify the first row in the stray light file for which the time is later than the current time relative to the start of the orbit
	for (unsigned i=1; i<m_timeFromFile.size(); i++) {
		if (m_timeFromFile[i] >= currentTime) {
			ifile = i;
			break;
		}
	}

	//Perform a linear interpolation between the stray light values in the rows on either side of the current time
	double timeFrac = static_cast<double>(currentTime - m_timeFromFile[ifile-1])/static_cast<double>(m_timeFromFile[ifile]-m_timeFromFile[ifile-1]);
	double strayLightFlux = (1.-timeFrac)*m_strayLightFromFile[ifile-1] + timeFrac*m_strayLightFromFile[ifile];
	//cout << ifile << " " << currentTime << " " << m_timeFromFile[ifile-1] << " " << m_timeFromFile[ifile] << " "
	//	 << m_strayLightFromFile[ifile-1] << " " << m_strayLightFromFile[ifile] << " " << timeFrac << " " << strayLightFlux << endl;

	return strayLightFlux;

}

void StrayLightGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	TimeConfiguration timeConf = data->getTimeConfiguration();
	double strayLightFlux = getStrayLightFlux(timeStep, timeConf);

	//Set the truth information
	Image* image = data->getImages().back();
	image->getTruthData()->setStrayLight(strayLightFlux*timeConf.getExposureTimeAsDouble());

	//Stray light has optical throughput already applied, so in order to avoid it being applied a second time by GlobalThroughputGenerator,
	//scale the flux here to the value it would have before applying optical throughput
	strayLightFlux /= m_integratedThroughput;

	//Add the stray light uniformly to each pixel of the image
	int ymin = data->doChargeTransferEoL() ? 0. : image->getYOffset(); //Include pixels below the subarray if end of life CTI is to be simulated
	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=ymin; iy<image->getYOffset()+image->getYDim(); iy++) {
			image->incrementPixelValue(ix,iy,strayLightFlux*timeConf.getExposureTimeAsDouble());
		}
	}

}
