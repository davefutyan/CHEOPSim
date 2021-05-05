/*
 * OrbitSimulator.cxx
 *
 *  Created on: 14 Feb 2014
 *      Author: futyand
 */

#include <fstream>

#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"
#include "boost/math/special_functions/sign.hpp"

#include "CreateFitsFile.hxx"
#include "REF_APP_Temperature.hxx"
#include "AUX_REF_Orbit.hxx"
#include "AUX_RES_Orbit.hxx"

#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "BarycentricOffset.hxx"

#include "source/include/TransitFluxModulator.hxx"
#include "OrbitSimulator.hxx"

OrbitSimulator::~OrbitSimulator() {
	if (m_attitude != nullptr) delete m_attitude;
	if (m_MPSVisitConstraints != nullptr) delete m_MPSVisitConstraints;
}

void OrbitSimulator::initialize(const ModuleParams& params) {

	m_ccdTemperatureVariation = params.GetAsBool("ccdTemperatureVariation"); //Boolean to indicate whether or not to apply CCD temperature variation
	m_ccdMeanTemperature = params.GetAsDouble("ccdMeanTemperature"); //Units: Kelvin
	m_ccdTemperatureAmplitude = params.GetAsDouble("ccdTemperatureAmplitude"); //Amplitude of sinusoidal variation
	m_ccdTemperaturePeriod = params.GetAsDouble("ccdTemperaturePeriod")*60; //Period of sinusoidal variation in seconds
	m_feeTemperatureVariation = params.GetAsBool("feeTemperatureVariation"); //Boolean to indicate whether or not to apply FEE temperature variation
	m_feeMeanTemperature = params.GetAsDouble("feeMeanTemperature"); //Units: Kelvin
	m_feeTemperatureAmplitude = params.GetAsDouble("feeTemperatureAmplitude"); //Amplitude of sinusoidal variation
	m_feeTemperaturePeriod = params.GetAsDouble("feeTemperaturePeriod")*60; //Period of sinusoidal variation in seconds
	m_telescopeTemperatureVariation = params.GetAsBool("telescopeTemperatureVariation"); //Boolean to indicate whether or not to apply telescope temperature variation
	m_telescopeMeanTemperature = params.GetAsDouble("telescopeMeanTemperature"); //Units: Kelvin
	m_telescopeTemperatureAmplitude = params.GetAsDouble("telescopeTemperatureAmplitude"); //Amplitude of sinusoidal variation
	m_telescopeTemperaturePeriod = params.GetAsDouble("telescopeTemperaturePeriod")*60; //Amplitude of sinusoidal variation in seconds
	m_telescopeTemperatureFromFile = params.GetAsBool("telescopeTemperatureFromFile"); //Boolean to indicate temperature from file or sinusoid
	m_telescopeTemperatureFilename = params.GetAsString("telescopeTemperatureFilename"); //Filename for telescope temperature
	m_rotateFOV = params.GetAsBool("rotateFOV"); //Boolean to indicate whether or not to rotate the FOV
	m_orbitFilename = params.GetAsString("orbitFilename"); //Filename for orbit data
	m_minAngleToOrbitalPlane = params.GetAsDouble("minAngleToOrbitalPlane"); //Minimum angle (degrees) between pointing direction and orbital plane for roll angle calculation
	m_attitudeCadence = params.GetAsInt("attitudeCadence"); //Cadence in seconds for writing out spacecraft attitude data to SCI_RAW_Attitude data structure
	m_orbitCadence = params.GetAsInt("orbitCadence"); //Cadence in minutes for writing out spacecraft orbit data to AUX_RES_Orbit data structure
	m_SAAMapFilename = params.GetAsString("SAAMapFilename"); //Filename for FITS file containing the SAA map
	m_saaFlagFromVisitConstraints = params.GetAsBool("saaFlagFromVisitConstraints"); //Flag to indicate whether to read the SAA from MPS_PRE_VisitConstraints rather than calculating it from the latitude/longitude in the orbit file

	if (m_telescopeTemperatureFromFile) readTemperatures(telescope);

}

void OrbitSimulator::doBegin(Data* data, bool fullFrame) {

	data->setPreSubArrayTimeSteps();

	//Set the maximum and minimum values for the range of telescope temperatures, used for the PSF breathing
	double max = -999.;
	double min = 999.;
	if (m_telescopeTemperatureFromFile) {
		for (unsigned i=0; i<m_telescopeTemperature.size(); i++) {
			if (m_telescopeTemperature[i] > max) max = m_telescopeTemperature[i];
			if (m_telescopeTemperature[i] < min) min = m_telescopeTemperature[i];
		}
	} else {
		max = m_telescopeMeanTemperature + m_telescopeTemperatureAmplitude;
		min = m_telescopeMeanTemperature - m_telescopeTemperatureAmplitude;
	}
	data->setBreathingTemperatureRange(min,max);

	//Define the line of sight vector in the inertial frame
	double ra = data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.*M_PI/180.;
	double dec = data->getFieldOfView()->getPointingDirection().getDeclination()/3600.*M_PI/180.;
	m_lineOfSight(0) = cos(dec)*cos(ra);
	m_lineOfSight(1) = cos(dec)*sin(ra);
	m_lineOfSight(2) = sin(dec);
	m_lineOfSight.normalize();
	//cout << "line of sight: " << endl << m_lineOfSight << endl << " " << ra*180./M_PI << " " << dec*180./M_PI << endl;
	cout << "line of sight (RA,dec): (" << ra*180./M_PI << "," << dec*180./M_PI << ")" << endl;

	//Read the orbit data
	readOrbit(data);

	//Print the angle between the pointing direction and the orbital plane
	//There must be at least 3 points in the orbit data in order to define the orbital plane
	if (m_orbitData.size()>=3) printAngleToOrbitalPlane(data->getFieldOfView()->getPointingDirection());

	//Read the SAA map
	readSAAMap();

	//Open the FITS file to contain the spacecraft attitude information
	Data::Visit visit = data->getVisit();
	string dirname = data->getOutputDirectory()+"/data";
	if (!boost::filesystem::exists(dirname)) boost::filesystem::create_directory(dirname);
	TimeConfiguration timeConf = data->getTimeConfiguration();
	UTC startTime = timeConf.getVisitStartTimeUTC();
	UTC endTime = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration());
	m_attitude = createFitsFile<SciRawAttitude>(dirname, startTime,
			                           VisitId(visit.m_progType,visit.m_progId,visit.m_reqId,visit.m_visitCtr), PassId());

	//Set the header keywords for the SciRawAttitude FITS file
	string specType = "G5"; //Default to spectral type to G5 if there are no stars in the FOV
	double GbandMag = 9.; //Default to 9th magnitude if there are no stars in the FOV
	double cheopsMag = 9.; //Default to 9th magnitude if there are no stars in the FOV
	double GbandMagErr = 0.;
	double cheopsMagErr = 0.;
	double TEff = 0.;
	if (data->getFieldOfView()->getStars().size()>0) {
		Star * targetStar = data->getFieldOfView()->getStars()[0];
		specType = targetStar->getSpectralTypeString();
		TEff = targetStar->getEffectiveTemperature();
		GbandMag = targetStar->getGaiaMagnitude();
		cheopsMag = targetStar->getCheopsMagnitude();
		GbandMagErr = targetStar->getGaiaMagnitudeError();
		cheopsMagErr = targetStar->getCheopsMagnitudeError();
	}
	m_attitude->setKeyProcChn("CHEOPSim");
	m_attitude->setKeyProcNum(visit.m_versionNum);
	m_attitude->setKeyArchRev(0);
	m_attitude->setKeyVStrtU(startTime);
	m_attitude->setKeyVStrtM(startTime.getMjd());
	m_attitude->setKeyVStopU(endTime+DeltaTime(1)); //extra second added for DFS input, since the OBT time of the last CE is slightly delayed w.r.t. the end of the last exposure
	m_attitude->setKeyVStopM((endTime+DeltaTime(1)).getMjd());
	m_attitude->setKeyTargname(visit.m_targName);
	m_attitude->setKeySpectype(specType);
	m_attitude->setKeyTEff(TEff);
	m_attitude->setKeyMagG(GbandMag);
	m_attitude->setKeyMagChps(cheopsMag);
	m_attitude->setKeyMagGerr(GbandMagErr);
	m_attitude->setKeyMagCerr(cheopsMagErr);
	m_attitude->setKeyPiName(visit.m_piName);
	m_attitude->setKeyPiUid(visit.m_piUid);
	m_attitude->setKeyObsCat(visit.m_obsCat);
	m_attitude->setKeyProgtype(visit.m_progType);
	m_attitude->setKeyProgId(visit.m_progId);
	m_attitude->setKeyReqId(visit.m_reqId);
	m_attitude->setKeyVisitctr(visit.m_visitCtr);
	m_attitude->setKeyObsid(visit.m_obsId);
	m_attitude->setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	m_attitude->setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
	m_attitude->setKeyPrpVst1(visit.m_prpFirst);
	m_attitude->setKeyPrpVstn(visit.m_prpLast);

	//Append the visit data fits file with the visit constraints table and set the header keywords
	m_MPSVisitConstraints = Append_MpsPreVisitconstraints(data->getMpsPreVisits());
	m_MPSVisitConstraints->setKeyProcNum(visit.m_versionNum);
	m_MPSVisitConstraints->setKeyArchRev(0);
	m_MPSVisitConstraints->setKeyVStrtU(data->getMpsPreVisits()->getKeyVStrtU());
	m_MPSVisitConstraints->setKeyVStrtM(data->getMpsPreVisits()->getKeyVStrtM());
	m_MPSVisitConstraints->setKeyVStopU(data->getMpsPreVisits()->getKeyVStopU());
	m_MPSVisitConstraints->setKeyVStopM(data->getMpsPreVisits()->getKeyVStopM());

	unsigned exposuresPerStack = timeConf.getExposuresPerStack();

	MpsPreVisitconstraints * visitConstraints_in = nullptr;
	if (m_saaFlagFromVisitConstraints) {
		visitConstraints_in = data->getMpsVisitConstraints();
		if (visitConstraints_in == nullptr) {
			throw runtime_error("Error in OrbitSimulator::doBegin: saaFlagFromVisitConstraints requested but null pointer to MPS_PRE_VisitsConstraints. Was an MPS_PRE_Visits file uploaded?");
		}
	}

	//double sl = 0.;
	for (int timeStep=-int(timeConf.getNumberOfPreSubArrayTimeSteps()); timeStep<int(timeConf.getNumberOfTimeSteps());  timeStep++) {
		//sl += 0.005; //dummy stray light values to generate test VisitConstraints with stray light filled

		SatelliteData * satData = data->getSatelliteData(timeStep,(!data->moduleIsActive("JitterProducer")&&!data->moduleIsActive("StrayLightGenerator")));
		if (satData->getRollAngles().size() != 0 && timeStep != 0) {
			throw runtime_error("Error in OrbitSimulator::doBegin: roll angles already defined for current timestep");
		}

		double timeSinceStart = data->getTimeConfiguration().getTimeSinceStart(timeStep);
		double timeSinceStartOfVisit = data->getTimeConfiguration().getTimeSinceStart(timeStep+int(timeConf.getNumberOfPreSubArrayTimeSteps()));

		//Set the telescope and DPU temperatures for the current time step
		if (m_ccdTemperatureVariation) satData->setCcdTemperature(getSinusoidTemperature(timeSinceStart,ccd));
		if (m_feeTemperatureVariation) satData->setFeeBiasTemperature(getSinusoidTemperature(timeSinceStart,fee));
		if (m_telescopeTemperatureVariation) {
			if (m_telescopeTemperatureFromFile) {
				satData->setTelescopeTemperature(getTemperatureFromFile(timeSinceStartOfVisit,telescope));
			} else {
				satData->setTelescopeTemperature(getSinusoidTemperature(timeSinceStart,telescope));
			}
		}
		satData->setDpuTemperature(getSinusoidTemperature(timeSinceStart,dpu));

		//Identify the two closest points to the current time in the OrbitData array
		double repetitionPeriod = data->getTimeConfiguration().getRepetitionPeriod();
		UTC currentTime = data->getTimeConfiguration().getStartTimeUTC() + DeltaTime(timeSinceStart+repetitionPeriod/2.);
		unsigned index;
		double mu;
		bool found;
		for (unsigned i=0; i<m_orbitData.size()-1; i++) {
			if (m_orbitData[i+1].utc() >= currentTime) {
				found = true;
				index = i;
				mu = (currentTime-m_orbitData[i].utc()).getSeconds()/(m_orbitData[i+1].utc()-m_orbitData[i].utc()).getSeconds();
				break;
			}
		}
		if (!found) throw runtime_error("Error in OrbitSimulator::doBegin: entry in orbit array corresponding to current time not found");

		//Calculate the angle between the line of sight and Earth limb,
		//interpolating linearly between the two closest points in time in the OrbitData array
		Eigen::Vector3d orbitPos;
		orbitPos(0) = m_orbitData[index].x()*(1.-mu)+m_orbitData[index+1].x()*mu;
		orbitPos(1) = m_orbitData[index].y()*(1.-mu)+m_orbitData[index+1].y()*mu;
		orbitPos(2) = m_orbitData[index].z()*(1.-mu)+m_orbitData[index+1].z()*mu;
		//cout << orbitPos(0) << " " << m_orbitData[index].x() << " " << m_orbitData[index+1].x() << " " << m_orbitData[index].utc().getUtc() << " " << m_orbitData[index+1].utc().getUtc() << " " << mu << endl;
		double angleToNormalToEarthLimb = acos((orbitPos.normalized()).dot(m_lineOfSight))*180./M_PI;
		//Assume atmosphere thickness is 100km
		double angleToEarthLimb = 90. - angleToNormalToEarthLimb + acos((TransitFluxModulator::kEarthRadius+100.)/orbitPos.norm())*180./M_PI;
		satData->setEarthLimbAngle(angleToEarthLimb);
		//cout << setprecision(5) << "Orbit position: " << orbitPos(0) << "\t" <<orbitPos(1) << "\t" <<orbitPos(2) << "\t" << currentTime.getUtc() << "\t" << angleToEarthLimb << endl;
		//cout << "Angle between line of sight and Earth limb: " << satData->getEarthLimbAngle() << endl;

		bool earthOccultationFlag = false;
		bool saaFlag = false;
		bool strayLightIsNull = true;
		double strayLight = 0.;
		double moonAngle,sunAngle;

		//Calculate the latitude and longitude of the satellite,
		//interpolating linearly between the two closest points in time in the OrbitData array
		//Use special treatment when the two points are on either side of the boundary between 180 and -180.
		double longitude1 = m_orbitData[index].longitude();
		double longitude2 = m_orbitData[index+1].longitude();
		if (fabs(longitude1)>150. && fabs(longitude2)>150. && longitude1*longitude2<0.) {
			if (longitude1<0.) longitude1 += 360.;
			if (longitude2<0.) longitude2 += 360.;
		}
		double longitude = longitude1*(1.-mu)+longitude2*mu;
		double latitude = m_orbitData[index].latitude()*(1.-mu)+m_orbitData[index+1].latitude()*mu;

		if (!m_saaFlagFromVisitConstraints) {

			//Define the moon and sun angles, interpolating linearly between the two closest points in time in the OrbitData array
			moonAngle = m_orbitData[index].moonAngle()*(1.-mu)+m_orbitData[index+1].moonAngle()*mu;
			sunAngle = m_orbitData[index].sunAngle()*(1.-mu)+m_orbitData[index+1].sunAngle()*mu;
			//cout << "Angle between line of sight and moon: " << satData->getMoonAngle() << endl;
			//cout << "Angle between line of sight and sun: " << satData->getSunAngle() << endl;

			//Determine the latitude and longitude indices of the closest point in the SAA map
			int ilat,ilong;
			for (ilat=0; ilat<90; ilat++) {
				if (m_SAAMap_lat[ilat]>latitude) {
					if (ilat>0 && m_SAAMap_lat[ilat]-latitude > latitude-m_SAAMap_lat[ilat-1]) ilat-=1;
					break;
				}
			}
			for (ilong=0; ilong<121; ilong++) {
				if (m_SAAMap_long[ilong]>longitude) {
					if (ilong>0 && m_SAAMap_long[ilong]-longitude > longitude-m_SAAMap_long[ilong-1]) ilong-=1;
					break;
				}
			}

			//Set the SAA flag
			saaFlag = m_SAAMap_value[ilat][ilong];

			//Set the Earth occultation flag
			earthOccultationFlag = angleToEarthLimb<0.;

			if (data->moduleIsActive("StrayLightGenerator")) {
				strayLight = satData->getStrayLightFlux();
				strayLightIsNull = false;
			}

		} else {

			bool found = false;
			double deltaTime_prev = 999.;
			UTC visitConstraintsUtc;
			visitConstraints_in->SetReadRow(1); //reset to the start of the file
			while (visitConstraints_in->ReadRow()) {
				double deltaTime = fabs((visitConstraints_in->getCellUtcTime() - currentTime).getSeconds());
				if (deltaTime < 300 && deltaTime < deltaTime_prev) {
					moonAngle = visitConstraints_in->getCellLosToMoonAngle();
					sunAngle = visitConstraints_in->getCellLosToSunAngle();
					saaFlag = visitConstraints_in->getCellSaaFlag();
					earthOccultationFlag = visitConstraints_in->getCellEarthOccultation();
					if (!visitConstraints_in->isNullStrayLight()) {
						strayLightIsNull = false;
						strayLight = visitConstraints_in->getCellStrayLight();
					}
					found = true;
					deltaTime_prev = deltaTime;
					visitConstraintsUtc = visitConstraints_in->getCellUtcTime();
				}
			}
			if (!found) {
				throw runtime_error("Error in OrbitSimulator::doBegin: No entry found in MPS_PRE_VisitConstraints file for UTC time "+currentTime.getUtc());
			}
			//cout << currentTime.getUtc() << " " << visitConstraintsUtc.getUtc() << " " << saaFlag << " " << earthOccultationFlag << " " << strayLight << endl;
			satData->setStrayLightFlux(strayLight);
			satData->setStrayLightFlag(strayLight > data->getVisit().m_strayLightThreshold);

		}

		satData->setMoonAngle(moonAngle);
		satData->setSunAngle(sunAngle);
		satData->setSAAFlag(saaFlag);
		satData->setEarthOccultationFlag(earthOccultationFlag);
		satData->setLatitude(latitude);
		satData->setLongitude(longitude);

		//Write out the mission planning visit constraint data once per stacked image and for the first timeStep corresponding to the start of the visit
		if (fabs(timeSinceStart/(exposuresPerStack*repetitionPeriod) - round(timeSinceStart/(exposuresPerStack*repetitionPeriod))) < 1E-6 || timeStep == -int(timeConf.getNumberOfPreSubArrayTimeSteps())) {

			Data::Visit visit = data->getVisit();

			m_MPSVisitConstraints->setCellProgrammeType(visit.m_progType);
			m_MPSVisitConstraints->setCellProgrammeId(visit.m_progId);
			m_MPSVisitConstraints->setCellRequestId(visit.m_reqId);
			m_MPSVisitConstraints->setCellObsid(visit.m_obsId);
			m_MPSVisitConstraints->setCellUtcTime(currentTime);
			m_MPSVisitConstraints->setCellMjdTime(currentTime.getMjd());
			//BarycentricOffset offset(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.,data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
			//m_MPSVisitConstraints->setCellBjdTime(currentTime.getMjd() + offset);
			m_MPSVisitConstraints->setCellLosToMoonAngle(moonAngle);
			m_MPSVisitConstraints->setCellLosToSunAngle(sunAngle);
			//m_MPSVisitConstraints->setCellScLatitude(latitude);
			//m_MPSVisitConstraints->setCellScLongitude(longitude);
			//m_MPSVisitConstraints->setCellLosToEarthlimbAngle(angleToEarthLimb);
			if (strayLightIsNull) {
				m_MPSVisitConstraints->setNullStrayLight();
			} else {
				m_MPSVisitConstraints->setCellStrayLight(strayLight);
			}
			m_MPSVisitConstraints->setCellEarthOccultation(earthOccultationFlag);
			m_MPSVisitConstraints->setCellSaaFlag(saaFlag);
			//m_MPSVisitConstraints->setCellStrayLight(sl); //to generate test VisitConstraints with stray light filled
			//m_MPSVisitConstraints->setCellSaaFlag(sl>.5); //to generate test VisitConstraints with saa sometimes true
			m_MPSVisitConstraints->WriteRow();
		}

		//Calculate the jittered roll angle and pointing direction for each second of the exposure
		setRollAnglesAndPointingDirections(data,timeStep);

	}

}

void OrbitSimulator::setRollAnglesAndPointingDirections(Data * data, int timeStep) const {

	TimeConfiguration timeConf = data->getTimeConfiguration();
	double timeSinceStart = timeConf.getTimeSinceStart(timeStep);

	int istart = 0;
	//If the timestep is the first one of the simulation, determine the number of seconds between the start of the visit
	//and the start of the first timestep to ensure that the entire visit is covered
	if (timeStep == -int(timeConf.getNumberOfPreSubArrayTimeSteps())) {
		int additionalSeconds = lround((timeConf.getStartTimeUTC() + DeltaTime(timeSinceStart) - timeConf.getVisitStartTimeUTC()).getSeconds());
		istart = -additionalSeconds;
		timeSinceStart -= additionalSeconds;
	}

	unsigned nSeconds = timeConf.getJitterSamplingsPerExposure();
    //If the timestep is the last one of the simulation, add two more seconds, to correspond to the end of the last second of the last exposure,
    //and to one second after that (required for DFS input, since the OBT time of the last CE is slightly delayed w.r.t. the end of the last exposure)
    if (timeStep == int(timeConf.getNumberOfTimeSteps())-1) nSeconds += 2;

	//Get the jitter APEs if available (i.e. if JitterProducer has been run)
	unsigned jitterSize = data->getSatelliteData(timeStep)->getJitterAPEs().size();
	bool doJitter = (jitterSize > 0);
	if (doJitter && jitterSize != nSeconds) {
		throw runtime_error("Error in OrbitSimulator::setRollAnglesAndPointingDirections: size of getJitterAPEs vector ("+
				to_string(jitterSize)+") does not match the exposure duration ("+to_string(nSeconds)+")");
	}

	//double rollAngle_old = -1.;
	for (int i=istart; i<int(nSeconds); i++) {

		UTC currentTime = timeConf.getStartTimeUTC() + DeltaTime(timeSinceStart);

		//Calculate the rotation matrix from the inertial frame to the satellite frame
		//If FOV rotation is deactivated, use a static matrix according to the orbit parameters at the start of the simulation.
		//This will result in constant roll angle but will ensure that the roll angle and pointing direction
		//(which can still vary with time if jitter is switched on) are assigned consistently
		Eigen::Matrix3d M_inertialToSat = inertialToSatMatrix(m_rotateFOV ? currentTime : timeConf.getStartTimeUTC());

		//Get the jitter APEs if available
		SatelliteData::APE ape(0.,0.,0.);
		if (doJitter) ape = data->getSatelliteData(timeStep)->getJitterAPEs()[i<0 ? 0:i];

		//Calculate the roll angle if requested
		double rollAngle = calculateRollAngle(M_inertialToSat,ape.getRollOffset());
		//cout << timeSinceStart << " " << rollAngle << " " << endl;

		//Calculate the jittered pointing direction
		SkyPosition pointingDirection = calculatePointingDirection(M_inertialToSat,ape);

		//			//Roll rate calculation: to allow plotting of roll rate, uncomment this and replace rollAngle with RollRate for m_attitude->setCellScRollAngle
		//			double rollRate = 0.;
		//			if (rollAngle_old >= 0.) {
		//				rollRate = fabs(rollAngle - rollAngle_old);
		//			} else if (timeStep > 0) {
		//				vector<double> rollAngles = data->getSatelliteData(timeStep-1,false,false)->getRollAngles();
		//				rollRate = fabs(rollAngle - rollAngles[rollAngles.size()-1]);
		//			}
		//			cout << timeSinceStart << " " << imageMidTime_utc.getMjd() << " " << rollAngle << " " << rollRate << endl;
		//			rollAngle_old = rollAngle;

		if (i >= 0 && i < int(timeConf.getJitterSamplingsPerExposure())) {
			//Add the roll angle to the list of jittered roll angles for the current exposure
			data->getSatelliteData(timeStep)->addRollAngle(rollAngle);
			//Add the pointing direction to the list of jittered pointing directions for the current exposure
			data->getSatelliteData(timeStep)->addPointingDirection(pointingDirection);
		}

		//Set the spacecraft attitude data with the user defined cadence for exposure durations > 1s
		//and with cadence corresponding to the repetition period for exposure durations <1s
		if ((timeConf.getExposureTimeAsDouble() < 1. || static_cast<unsigned>(lround(timeSinceStart)) % m_attitudeCadence == 0) &&
				currentTime>=timeConf.getVisitStartTimeUTC()) {
			m_attitude->setCellUtcTime(currentTime);
			m_attitude->setCellMjdTime(currentTime.getMjd());
			m_attitude->setCellObtTime(currentTime.getObt());
			m_attitude->setCellScRollAngle(360.-rollAngle); // In real data, the roll angle is defined looking down onto the CCD rather than up towards the target.
			m_attitude->setCellScRa(pointingDirection.getRightAscension()/3600.);
			m_attitude->setCellScDec(pointingDirection.getDeclination()/3600.);
			m_attitude->WriteRow();
		}

		timeSinceStart += 1.;

	}

}

void OrbitSimulator::readTemperatures(COMPONENT component) {

	if (component == telescope) {

		cout << "Reading telescope temperature from " << m_telescopeTemperatureFilename << endl;

	    if (m_telescopeTemperatureFilename.substr(m_telescopeTemperatureFilename.find_last_of(".")+1) == "txt") {

	    	//user uploaded temperature file, to be read from execution directory
	    	ifstream temperature_in(m_telescopeTemperatureFilename, ios::in);
	    	if (!temperature_in.good()) throw runtime_error("Error in OrbitSimulator::readTemperatures: Error opening file "+m_telescopeTemperatureFilename);
	    	string line;
	    	while (getline(temperature_in, line)) {
	    		if(line[0] != '#') {
	    			double temperature, time;
	    			stringstream(line) >> time >> temperature;
	    			m_telescopeTime.push_back(time);
	    			m_telescopeTemperature.push_back(temperature+273.15);
	    		}
	    	}
	    	temperature_in.close();

	    } else {

	    	//Open the temperature fits file
	    	string filename = string(getenv("CHEOPS_SW"))+"/resources/"+m_telescopeTemperatureFilename;
	    	RefAppTemperature * temperature_file = new RefAppTemperature(filename,"READONLY");

	    	m_telescopeTime.clear();
	    	m_telescopeTemperature.clear();
	    	//Read in the temperature data
	    	while (temperature_file->ReadRow()) {
	    		m_telescopeTime.push_back(temperature_file->getCellTime());
	    		m_telescopeTemperature.push_back(temperature_file->getCellTemperature()+273.15);
	    	}
	    	delete temperature_file;

	    }

	} else {
		throw runtime_error("Error in OrbitSimulator::readTemperature: only reading of telescope temperature is supported");
	}

}

double OrbitSimulator::getTemperatureFromFile(double timeSinceStart, COMPONENT component) const {

	if (component == telescope) {

		unsigned fileDuration = (m_telescopeTime.size()-1)*(m_telescopeTime[1]-m_telescopeTime[0]);
		for (unsigned i=0; i<m_telescopeTime.size(); i++) {
			if (m_telescopeTime[i] >= static_cast<unsigned>(lround(timeSinceStart))%fileDuration) return m_telescopeTemperature[i];
		}

	} else {
		throw runtime_error("Error in OrbitSimulator::getTemperatureFromFile: only reading of telescope temperature is supported");
	}

	throw runtime_error("Error in OrbitSimulator::getTemperatureFromFile: no temperature found for time step");

}

double OrbitSimulator::getSinusoidTemperature(double timeSinceStart, COMPONENT component) const {

	double mean, amplitude, period;
	if (component == telescope) {
		mean = m_telescopeMeanTemperature;
		amplitude = m_telescopeTemperatureAmplitude;
		period = m_telescopeTemperaturePeriod;
	} else if (component == fee) {
		mean = m_feeMeanTemperature;
		amplitude = m_feeTemperatureAmplitude;
		period = m_feeTemperaturePeriod;
	} else if (component == ccd) {
		mean = m_ccdMeanTemperature;
		amplitude = m_ccdTemperatureAmplitude;
		period = m_ccdTemperaturePeriod;
	} else if (component == dpu) {
		mean = SatelliteData::kDefaultDpuTemperature;
		amplitude = kDpuTemperatureAmplitude;
		period = kCheopsOrbitPeriod;
	} else {
		throw runtime_error("Error in OrbitSimulator::getSinusoidTemperature: OrbitSimulator::COMPONENT invalid");
	}

	double numberOfPeriods = timeSinceStart/period;
	double periodFraction = numberOfPeriods - trunc(numberOfPeriods);
	//cout << component << " " << period << " " << timeSinceStart << " " << periodFraction << " " << (mean + amplitude*sin(periodFraction*2*M_PI)) << endl;
	return (mean + amplitude*sin(periodFraction*2*M_PI));

}

void OrbitSimulator::readOrbit(Data * data) {

	//Determine the UTC start and end time of the simulation, including margins of m_orbitCadence
	//before and after the start and end of the simulation, to ensure that the full frame image
	//and the time offsets for HK data are covered
	TimeConfiguration timeConf = data->getTimeConfiguration();
	UTC startTime = timeConf.getVisitStartTimeUTC() - DeltaTime(static_cast<double>(m_orbitCadence*60.));
	UTC endTime = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration()) + DeltaTime(static_cast<double>(m_orbitCadence*60.));

	//Open the input fits file containing the orbit data
	string filename = string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_orbitFilename;
	if (m_orbitFilename == "UserOrbit.fits") filename = m_orbitFilename;

	//Read the orbit data from the file
	cout << "Reading orbit data" << endl;
	vector<OrbitData> orbitDataFromFile;

	if (m_orbitFilename.find("AUX_REF_Orbit") != string::npos) {
		//AUX_REF_Orbit files are provided in a directory containing multiple files. Need to identify the one or two files corresponding to the visit

		string orbitDirectory = string(getenv("CHEOPS_SW"))+"/resources/"+m_orbitFilename;
		if (!boost::filesystem::exists(orbitDirectory)) {
			throw runtime_error("Error in OrbitSimulator::readOrbit: directory "+orbitDirectory+" does not exist");
		}

		//Obtain a sorted list of AUX_REF_Orbit files in the directory
		vector<string> orbitFiles;
		boost::filesystem::directory_iterator it(orbitDirectory), eod;
		BOOST_FOREACH(boost::filesystem::path const &fs_path, std::make_pair(it, eod)) {
			if(is_regular_file(fs_path)) {
				string filename = fs_path.string();
				if(filename.substr(filename.find_last_of(".")+1) == "fits" && filename.find("AUX_REF_Orbit") != string::npos) {
					orbitFiles.push_back(filename);
				}
			}
		}
		std::sort(orbitFiles.begin(), orbitFiles.end());

		for (vector<string>::const_iterator file = orbitFiles.begin(); file!=orbitFiles.end(); ++file)  {

			AuxRefOrbit * orbit_in = new AuxRefOrbit(*file,"READONLY");

			//Determine if either the start or end of the visit lies within the validity time of the file
			//cout << startTime.getUtc() << " " << endTime.getUtc() << " " << orbit_in->getKeyVStrtU().getUtc() << " " << orbit_in->getKeyVStopU().getUtc() << endl;
			if ((startTime >= orbit_in->getKeyVStrtU() && startTime < orbit_in->getKeyVStopU()) ||
					(endTime >= orbit_in->getKeyVStrtU() && endTime < orbit_in->getKeyVStopU())) {
				cout << *file << endl;

				//Read in the data for times covering the visit
				if (file != orbitFiles.begin()) orbit_in->ReadRow(); //If this is not the first file, skip the first row because the first row of each file duplicates the last row of the previous (midnight entry)
				while (orbit_in->ReadRow()) {
					if (orbit_in->getCellEpoch() >= startTime-DeltaTime(60) &&
							orbit_in->getCellEpoch() <= endTime+DeltaTime(m_orbitCadence*60)) {
						//cout << (orbit_in->getCellEpoch()).getUtc() << " " << (startTime-DeltaTime(60)).getUtc() << " " << (endTime+DeltaTime(60)).getUtc() << endl;
						orbitDataFromFile.push_back(OrbitData(orbit_in->getCellEpoch(),orbit_in->getCellX(),orbit_in->getCellY(),orbit_in->getCellZ(),
								orbit_in->getCellLatitude(),orbit_in->getCellLongitude()));
					}
				}

			}

			delete orbit_in;

		}


	} else {
		//Older AUX_RES_Orbit files are single files corresponding to the calendar year. Re-use the same file for other years.

		AuxResOrbit * orbit_in = new AuxResOrbit(filename,"READONLY");

		//Read the year from the file and determine the year offset w.r.t. the simulation year
		orbit_in->ReadRow();
		DeltaTime yearOffset = DeltaTime(static_cast<double>((orbit_in->getCellEpoch().getYear() - startTime.getYear()) * 3600*24*365));

		do {
			if (orbit_in->getCellEpoch()-yearOffset >= startTime-DeltaTime(60) &&
					orbit_in->getCellEpoch()-yearOffset <= endTime+DeltaTime(m_orbitCadence*60)) {
				//cout << (orbit_in->getCellEpoch()-yearOffset).getUtc() << " " << (startTime-DeltaTime(60)).getUtc() << " " << (endTime+DeltaTime(60)).getUtc() << endl;
				orbitDataFromFile.push_back(OrbitData(orbit_in->getCellEpoch()-yearOffset,orbit_in->getCellX(),orbit_in->getCellY(),orbit_in->getCellZ(),
						orbit_in->getCellLatitude(),orbit_in->getCellLongitude()));
			}
		} while (orbit_in->ReadRow());

		//Close the input file
		delete orbit_in;

	}

	if (orbitDataFromFile.size() == 0) {
		throw runtime_error("Error in OrbitSimulator::readOrbit: no entries found in orbit file compatible with the simulation start and end time");
	}

	//Initialize the orbit data vector as a list of UTC times at 1 minute intervals between the start and end of the simulation,
	//plus margins of m_orbitCadence before and after the start and end of the simulation
	UTC utcTime = startTime;
	unsigned numberOfMinutes = static_cast<unsigned>(lround(timeConf.getDuration()))/kCheopsOrbitTimeStep + 2*m_orbitCadence;
	for (unsigned i=0; i<=numberOfMinutes; i++) {
		m_orbitData.push_back(OrbitData(utcTime));
		utcTime += DeltaTime(kCheopsOrbitTimeStep);
	}

	//Set the spacecraft orbit position for each minute by linear interpolation
	//from the time resolution of the input file to that of the simulation
	setOrbitPositions(orbitDataFromFile);

	//Set the spacecraft orbit velocity for each minute
	setOrbitVelocities();

	if (!m_saaFlagFromVisitConstraints) {
		//Set the angle between the moon and the pointing direction for each minute
		setMoonAngles(startTime,endTime);

		//Set the angle between the sun and the pointing direction for each minute
		setSunAngles(startTime,endTime);
	}

	//Open the output file and write the orbit data at intervals corresponding to the user defined cadence
	Data::Visit visit = data->getVisit();
	string dirname = data->getOutputDirectory()+"/aux";
	if (!boost::filesystem::exists(dirname)) boost::filesystem::create_directory(dirname);
	AuxResOrbit * orbit_out = new AuxResOrbit(buildFileName(dirname, startTime, VisitId(),
			PassId(), std::string(), AuxResOrbit::getExtName()), "CREATE");
	orbit_out->setKeyArchRev(0);
	orbit_out->setKeyProcNum(visit.m_versionNum);
	orbit_out->setKeyProcChn("CHEOPSim");
	orbit_out->setKeyVStrtU(startTime);
	orbit_out->setKeyVStopU(endTime);
	orbit_out->setKeyUseSttm(startTime);
	orbit_out->setKeyUseEntm(endTime);
	for (unsigned i=0; i<m_orbitData.size(); i+=m_orbitCadence) { //Write out with user defined cadence
		orbit_out->setCellEpoch(m_orbitData[i].utc());
		double orbitRadius = pow(m_orbitData[i].x(),2) + pow(m_orbitData[i].y(),2) + pow(m_orbitData[i].z(),2);
		if (orbitRadius<1.) {//Can happen if the simulation does not start and finish in the same year (orbitDataFromFile corresponds to calendar year)
			orbit_out->setNullX();
			orbit_out->setNullY();
			orbit_out->setNullZ();
			orbit_out->setNullLatitude();
			orbit_out->setNullLongitude();
			orbit_out->setNullXDot();
			orbit_out->setNullYDot();
			orbit_out->setNullZDot();
		} else {
			orbit_out->setCellX(m_orbitData[i].x());
			orbit_out->setCellY(m_orbitData[i].y());
			orbit_out->setCellZ(m_orbitData[i].z());
			orbit_out->setCellLatitude(m_orbitData[i].latitude());
			orbit_out->setCellLongitude(m_orbitData[i].longitude());
			orbit_out->setCellXDot(m_orbitData[i].vx());
			orbit_out->setCellYDot(m_orbitData[i].vy());
			orbit_out->setCellZDot(m_orbitData[i].vz());
		}
		orbit_out->WriteRow();
	}
	delete orbit_out;

}

void OrbitSimulator::setOrbitPositions(const vector<OrbitData> orbitData_ref) {

	for (unsigned i=0; i<m_orbitData.size(); i++) {
		for (unsigned j=1; j<orbitData_ref.size(); j++) {
			if (orbitData_ref[j].utc() >= m_orbitData[i].utc()) {
				double mu = (m_orbitData[i].utc()-orbitData_ref[j-1].utc()).getSeconds()/(orbitData_ref[j].utc()-orbitData_ref[j-1].utc()).getSeconds();
				double x = orbitData_ref[j-1].x()*(1.-mu)+orbitData_ref[j].x()*mu;
				double y = orbitData_ref[j-1].y()*(1.-mu)+orbitData_ref[j].y()*mu;
				double z = orbitData_ref[j-1].z()*(1.-mu)+orbitData_ref[j].z()*mu;
				//interpolated altitude will be up to 90km smaller than true altitude for 5 minute sampling interval: maximum altitude reduction (half way between samplings) is given by r(1-cos(theta/2)) where r is the true orbit radius (rEarth+altitude) and theta is the angle swept between samplings.
				//cout << sqrt(x*x+y*y+z*z)-6371. << " " << sqrt(orbitData_ref[j].x()*orbitData_ref[j].x()+orbitData_ref[j].y()*orbitData_ref[j].y()+orbitData_ref[j].z()*orbitData_ref[j].z())-6371. << " " << sqrt(orbitData_ref[j-1].x()*orbitData_ref[j-1].x()+orbitData_ref[j-1].y()*orbitData_ref[j-1].y()+orbitData_ref[j-1].z()*orbitData_ref[j-1].z())-6371. << " " << mu << endl;
				double latitude = orbitData_ref[j-1].latitude()*(1.-mu)+orbitData_ref[j].latitude()*mu;
				double longitude = orbitData_ref[j-1].longitude()*(1.-mu)+orbitData_ref[j].longitude()*mu;
				m_orbitData[i].setPosition(x,y,z);
				m_orbitData[i].setLatLong(latitude,longitude);
				break;
			}
		}
	}

}

void OrbitSimulator::setOrbitVelocities() {

	for (unsigned i=0; i<m_orbitData.size(); i++) {

		double vx,vy,vz;
		if (i==0) {
			vx = (m_orbitData[i+1].x() - m_orbitData[i].x())/kCheopsOrbitTimeStep;
			vy = (m_orbitData[i+1].y() - m_orbitData[i].y())/kCheopsOrbitTimeStep;
			vz = (m_orbitData[i+1].z() - m_orbitData[i].z())/kCheopsOrbitTimeStep;
		} else if (i==m_orbitData.size()-1) {
			vx = (m_orbitData[i].x() - m_orbitData[i-1].x())/kCheopsOrbitTimeStep;
			vy = (m_orbitData[i].y() - m_orbitData[i-1].y())/kCheopsOrbitTimeStep;
			vz = (m_orbitData[i].z() - m_orbitData[i-1].z())/kCheopsOrbitTimeStep;
		} else {
			vx = (m_orbitData[i+1].x() - m_orbitData[i-1].x())/(2*kCheopsOrbitTimeStep);
			vy = (m_orbitData[i+1].y() - m_orbitData[i-1].y())/(2*kCheopsOrbitTimeStep);
			vz = (m_orbitData[i+1].z() - m_orbitData[i-1].z())/(2*kCheopsOrbitTimeStep);
		}
		m_orbitData[i].setVelocity(vx,vy,vz);

		//double radius = sqrt(m_orbitData[i].x()*m_orbitData[i].x() + m_orbitData[i].y()*m_orbitData[i].y() + m_orbitData[i].z()*m_orbitData[i].z());
		//double speed = sqrt(vx*vx + vy*vy + vz*vz);
		//cout << m_orbitData[i].utc().getFileNamePattern() << " " << radius << " " << speed << " " << m_orbitData[i].x() << " " << m_orbitData[i].y() << " " << m_orbitData[i].z() << " " << m_orbitData[i].vx() << " " << m_orbitData[i].vy() << " " << m_orbitData[i].vz() << endl;

	}

}

void OrbitSimulator::setMoonAngles(UTC startTime, UTC endTime) {

	//Open the input fits file containing the moon coordinates as a function of time
	string filename = string(getenv("CHEOPS_TESTDATA"))+"/resources/CHEOPS-ESOC-TRA-013_moon.fits";
	AuxResOrbit * moonPos_in = new AuxResOrbit(filename,"READONLY");

	//Read the year from the file and determine the year offset w.r.t. the simulation year
	cout << "Reading Moon position data" << endl;
	moonPos_in->ReadRow();
	DeltaTime yearOffset = DeltaTime(static_cast<double>((moonPos_in->getCellEpoch().getYear() - (startTime+DeltaTime(static_cast<double>(m_orbitCadence*60.))).getYear()) * 3600*24*365));

	//Read the Moon coordinates from the file
	vector<OrbitData> moonData;
	do {
		if (moonPos_in->getCellEpoch()-yearOffset >= startTime-DeltaTime(60) &&
			moonPos_in->getCellEpoch()-yearOffset <= endTime+DeltaTime(60)) {
			moonData.push_back(OrbitData(moonPos_in->getCellEpoch()-yearOffset,
					moonPos_in->getCellX(),moonPos_in->getCellY(),moonPos_in->getCellZ()));
		}
	} while (moonPos_in->ReadRow());

	delete moonPos_in;

	if (moonData.size() == 0) {
		throw runtime_error("Error in OrbitSimulator::setMoonAngles: no entries found in Moon coordinates file compatible with the simulation start and end time");
	}

	for (unsigned i=0; i<m_orbitData.size(); i++) {
		Eigen::Vector3d orbitPos, moonPos;
		orbitPos(0) = m_orbitData[i].x();
		orbitPos(1) = m_orbitData[i].y();
		orbitPos(2) = m_orbitData[i].z();
		for (unsigned j=1; j<moonData.size(); j++) {
			if (moonData[j].utc() >= m_orbitData[i].utc()) {
				double mu = (m_orbitData[i].utc()-moonData[j-1].utc()).getSeconds()/(moonData[j].utc()-moonData[j-1].utc()).getSeconds();
				moonPos(0) = moonData[j-1].x()*(1.-mu)+moonData[j].x()*mu;
				moonPos(1) = moonData[j-1].y()*(1.-mu)+moonData[j].y()*mu;
				moonPos(2) = moonData[j-1].z()*(1.-mu)+moonData[j].z()*mu;
				break;
			}
		}
		//Take into account parallax due to orbit of satellite around Earth
		m_orbitData[i].setMoonAngle(acos(((moonPos - orbitPos).normalized()).dot(m_lineOfSight))*180./M_PI);
	}

}

void OrbitSimulator::setSunAngles(UTC startTime, UTC endTime) {

	//Open the input fits file containing the sun coordinates as a function of time
	string filename = string(getenv("CHEOPS_TESTDATA"))+"/resources/CHEOPS-ESOC-TRA-014_sun.fits";
	AuxResOrbit * sunPos_in = new AuxResOrbit(filename,"READONLY");

	//Read the year from the file and determine the year offset w.r.t. the simulation year
	cout << "Reading Sun position data" << endl;
	sunPos_in->ReadRow();
	DeltaTime yearOffset = DeltaTime(static_cast<double>((sunPos_in->getCellEpoch().getYear() - (startTime+DeltaTime(static_cast<double>(m_orbitCadence*60.))).getYear()) * 3600*24*365));

	//Read the Sun coordinates from the file
	vector<OrbitData> sunData;
	do {
		if (sunPos_in->getCellEpoch()-yearOffset >= startTime-DeltaTime(60) &&
			sunPos_in->getCellEpoch()-yearOffset <= endTime+DeltaTime(60)) {
			sunData.push_back(OrbitData(sunPos_in->getCellEpoch()-yearOffset,
					sunPos_in->getCellX(),sunPos_in->getCellY(),sunPos_in->getCellZ()));
		}
	} while (sunPos_in->ReadRow());

	delete sunPos_in;

	if (sunData.size() == 0) {
		throw runtime_error("Error in OrbitSimulator::setSunAngles: no entries found in Sun coordinates file compatible with the simulation start and end time");
	}

	for (unsigned i=0; i<m_orbitData.size(); i++) {
		Eigen::Vector3d sunPos;
		for (unsigned j=1; j<sunData.size(); j++) {
			if (sunData[j].utc() >= m_orbitData[i].utc()) {
				double mu = (m_orbitData[i].utc()-sunData[j-1].utc()).getSeconds()/(sunData[j].utc()-sunData[j-1].utc()).getSeconds();
				sunPos(0) = sunData[j-1].x()*(1.-mu)+sunData[j].x()*mu;
				sunPos(1) = sunData[j-1].y()*(1.-mu)+sunData[j].y()*mu;
				sunPos(2) = sunData[j-1].z()*(1.-mu)+sunData[j].z()*mu;
				break;
			}
		}
		m_orbitData[i].setSunAngle(acos((sunPos.normalized()).dot(m_lineOfSight))*180./M_PI);
	}

}

void OrbitSimulator::readSAAMap() {

	ofstream referenceFilesList("reference_files.txt", ios::app);
	referenceFilesList << m_SAAMapFilename << endl;
	referenceFilesList.close();

	//Open the SAA map input fits file and read in the data
	m_SAAMap = new ExtAppSaamap(string(getenv("CHEOPS_SW"))+"/resources/"+m_SAAMapFilename,"READONLY");
	int oldlat = -999;
	int ilat = -1;
	int ilong = 0;
	while (m_SAAMap->ReadRow()) {
		if (m_SAAMap->getCellLatitude() != oldlat) {
			ilat++;
			m_SAAMap_lat[ilat] = m_SAAMap->getCellLatitude();
			ilong = 0;
			oldlat = m_SAAMap->getCellLatitude();
		}
		if (ilat>89 || ilong>120) throw runtime_error("Error in OrbitSimulator::doBegin: SAA map dimension exeeds 90(lat)x121(long)");
		m_SAAMap_value[ilat][ilong] = m_SAAMap->getCellSaaFlag();
		if (ilat == 0) m_SAAMap_long[ilong] = m_SAAMap->getCellLongitude();
		ilong++;
	}
	delete m_SAAMap;

}

void OrbitSimulator::printAngleToOrbitalPlane(SkyPosition pointingDirection) {

	if (m_orbitData.size()<3) throw runtime_error("Error in OrbitSimulator::printAngleToOrbitalPlane: there must be at least 3 points in the orbit data in order to define the orbital plane");

	//Compute the vector perpendicular to the orbital plane
	Eigen::Vector3d pos1(m_orbitData[0].x(),m_orbitData[0].y(),m_orbitData[0].z());
	Eigen::Vector3d pos2(m_orbitData[m_orbitData.size()/2].x(),m_orbitData[m_orbitData.size()/2].y(),m_orbitData[m_orbitData.size()/2].z());
	Eigen::Vector3d pos3(m_orbitData[m_orbitData.size()-1].x(),m_orbitData[m_orbitData.size()-1].y(),m_orbitData[m_orbitData.size()-1].z());
	Eigen::Vector3d normalToOrbitalPlane = ((pos3-pos2).cross(pos2-pos1)).normalized();
	double normalToOrbitalPlane_ra = atan2(static_cast<double>(normalToOrbitalPlane(1)),static_cast<double>(normalToOrbitalPlane(0)))*180./M_PI;
	double normalToOrbitalPlane_dec = asin(static_cast<double>(normalToOrbitalPlane(2)))*180./M_PI;
	if (normalToOrbitalPlane_ra<0) normalToOrbitalPlane_ra += 360;
	cout << "Normal to orbital plane (RA, dec): (" << normalToOrbitalPlane_ra << "," << normalToOrbitalPlane_dec << ")" << endl;

	//Calculate the angle between the line of sight and the normal to the orbital plane
	double angleToPlaneNormal = acos(m_lineOfSight.dot(normalToOrbitalPlane))*180./M_PI;
	if (angleToPlaneNormal>90) angleToPlaneNormal = 180 - angleToPlaneNormal;
	double angleToPlane = 90 - angleToPlaneNormal;
	cout << "Angle between line of sight and orbital plane: " << angleToPlane << endl;

}

Eigen::Matrix3d OrbitSimulator::inertialToSatMatrix(UTC currentTime) const {

	//cout << "UTC time: " << currentTime.getUtc() << endl;

	//Identify the two closest points to the current time in the OrbitData array
	unsigned index;
	double mu;
	bool found;
	for (unsigned i=0; i<m_orbitData.size()-1; i++) {
		if (m_orbitData[i+1].utc() >= currentTime) {
			found = true;
			index = i;
			mu = (currentTime-m_orbitData[i].utc()).getSeconds()/(m_orbitData[i+1].utc()-m_orbitData[i].utc()).getSeconds();
			break;
		}
	}
	if (!found) throw runtime_error("Error in OrbitSimulator::inertialToSatMatrix: entry in orbit array corresponding to current time not found");

	//Define the spacecraft position vector in the inertial frame,
	//interpolating linearly between the two closest points in time in the OrbitData array
	Eigen::Vector3d pos_J2000;
	pos_J2000(0) = m_orbitData[index].x()*(1.-mu)+m_orbitData[index+1].x()*mu;
	pos_J2000(1) = m_orbitData[index].y()*(1.-mu)+m_orbitData[index+1].y()*mu;
	pos_J2000(2) = m_orbitData[index].z()*(1.-mu)+m_orbitData[index+1].z()*mu;
	//cout << "pos_J2000: " << endl << pos_J2000 << endl;

	//Define the spacecraft velocity vector in the inertial frame,
	//interpolating linearly between the two closest points in time in the OrbitData array
	Eigen::Vector3d vel_J2000;
	vel_J2000(0) = m_orbitData[index].vx()*(1.-mu)+m_orbitData[index+1].vx()*mu;
	vel_J2000(1) = m_orbitData[index].vy()*(1.-mu)+m_orbitData[index+1].vy()*mu;
	vel_J2000(2) = m_orbitData[index].vz()*(1.-mu)+m_orbitData[index+1].vz()*mu;
	//cout << "vel_J2000: " << endl << vel_J2000 << endl;

	//The following code was converted to C++ from psudo-code provided by Carlos Corral Van Damme

	//Calculate the elements of the rotation matrix to convert
	//from the inertial frame to the LVLH frame
	Eigen::Vector3d uZ = -pos_J2000.normalized();
	Eigen::Vector3d uY = (uZ.cross(vel_J2000)).normalized();
	Eigen::Vector3d uX = uY.cross(uZ);

	//Matrix from inertial to LVLH frame
	Eigen::Matrix3d LVLH_J2000;
	LVLH_J2000 << uX,uY,uZ;
	//LVLH_J2000.transposeInPlace();
	//if (timeSinceStart==4490) {
	//	cout << "uX: " << endl << uX << endl;
	//	cout << "uY: " << endl << uY << endl;
	//	cout << "uZ: " << endl << uZ << endl;
	//	cout << "LVLH_J2000: " << endl << LVLH_J2000 << endl;
	//}

	//Calculate the elements of the rotation matrix to convert
	//from the inertial frame to the orbital frame
	Eigen::Vector3d Yop(0,-1,0);
	Eigen::Vector3d NEq(0,0,1);
	Eigen::Vector3d Yop_J2000 = LVLH_J2000 * Yop;
	Eigen::Vector3d Xop_J2000 = (Yop_J2000.cross(NEq)).normalized();
	Eigen::Vector3d Zop_J2000 = Xop_J2000.cross(Yop_J2000);

	//Matrix from inertial to orbital frame
	Eigen::Matrix3d M_inertialToOrbit;
	M_inertialToOrbit << Xop_J2000,Yop_J2000,Zop_J2000;
	M_inertialToOrbit.transposeInPlace();
	//if (timeSinceStart==4490) cout << "M_inertialToOrbit: " << endl << M_inertialToOrbit << endl;

	//Compute line of sight and Earth->S/C vectors in orbital frame
	Eigen::Vector3d x = M_inertialToOrbit * m_lineOfSight;
	Eigen::Vector3d p = M_inertialToOrbit * pos_J2000.normalized();

	//Ensures a min value of the lineOfSight Y component to ensure the Y axes always has
	//non-zero norm (lineOfSight Y maps into Y axes X/Z components in orb. frame)
	Eigen::Vector3d xm = x;
	double xm1 = static_cast<double>(xm(1));
	xm(1) = boost::math::sign(xm1) * max(fabs(xm1), sin(m_minAngleToOrbitalPlane*M_PI/180.));
	//cout << x << endl << endl << xm << endl << endl << endl;
	Eigen::Vector3d ym = p.cross(xm.normalized());

	//Z is computed from real lineOfSight and modified Y
	Eigen::Vector3d z = (x.cross(ym)).normalized();

	//Y is computed from real X and Z
	Eigen::Vector3d y = z.cross(x);

	//Rotation matrix from orbit to satellite frame
	Eigen::Matrix3d M_OrbitToSat;
	M_OrbitToSat << x,y,z;
	M_OrbitToSat.transposeInPlace();

	//Rotation matrix from inertial to satellite frame
	Eigen::Matrix3d M_inertialToSat = M_OrbitToSat * M_inertialToOrbit;
	//if (timeSinceStart==4490)	cout << "M_OrbitToSat: " << endl << M_OrbitToSat << endl << endl;
	//if (timeSinceStart==4490)	cout << "M_inertialToSat: " << endl << M_inertialToSat << endl << endl;
	//cout << currentTime.getUtc() << " " << M_inertialToSat(2,1) << " " << M_inertialToSat(1,1) << " " << xm1 << " " << xm(1) << endl;

	return M_inertialToSat;
}

double OrbitSimulator::calculateRollAngle(Eigen::Matrix3d M_inertialToSat, double rollOffset) const {

	//Calculate the roll angle from the satellite frame rotation matrix elements
	double rollAngle = atan2(static_cast<double>(M_inertialToSat(1,2)), static_cast<double>(M_inertialToSat(2,2)));

	//Eigen::Vector3d ea = M_inertialToSat.eulerAngles(3, 2, 1);
	//double rollAngle2 = ea(1);
	//cout << (rollAngle2 - rollAngle) << endl;

	//Convert from radians to degrees
	rollAngle *= (180./M_PI);
	if (rollAngle < 0.) rollAngle += 360.;

	//Add the offset due to jitter APE (converted from arcseconds to degrees)
	rollAngle += rollOffset/3600.;

	return rollAngle;

}

SkyPosition OrbitSimulator::calculatePointingDirection(Eigen::Matrix3d M_inertialToSat, SatelliteData::APE ape) const {

	// Note: If we simply add the APE values on to the RA/DEC coordinates we will
	// be distorting the jitter near the celestial poles due to non-linear
	// coordinate warping. Hence, we instead shift the ref vector with a
	// rotation matrix in satellite frame.

	//Convert APE values from arcseconds to radians
	double ape_X = (ape.m_X/3600.)*M_PI/180.;
	double ape_Y = (ape.m_Y/3600.)*M_PI/180.;
	double ape_Z = (ape.m_Z/3600.)*M_PI/180.;
	//Compute approximation of rotation matrix to shift a vector by APE
	Eigen::Matrix3d dR_ape;
	dR_ape <<  1,     -ape_Z,   ape_Y,
	           ape_Z,  1,      -ape_X,
	          -ape_Y,  ape_X,   1;
	//std::cout << dR_ape;

	//Shift pointing direction vector by APE
	Eigen::Vector3d m_ref = M_inertialToSat.transpose() * dR_ape * M_inertialToSat * m_lineOfSight;
	m_ref.normalize(); //correct for normalization error

	//Convert perturbed J2000 unit vector back to RA/DEC
	Eigen::Vector2d xy;
	xy << m_ref(0), m_ref(1);
	xy.normalize();
	double ra  = acos(xy(0))*(xy(1)<0.?-1.:1.);
	double dec = asin(m_ref(2));

	//Convert from radians to degrees and set RA in range [0,360]
	dec *= (180./M_PI);
	ra *= (180./M_PI);
	if (ra < 0.) ra += 360.;
	if (ra > 360.) ra -= 360.;

	return SkyPosition(ra*3600.,dec*3600.);

}
