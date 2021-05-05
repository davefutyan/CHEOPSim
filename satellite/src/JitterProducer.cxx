/*
 * JitterProducer.cxx
 *
 *  Created on: 13 Feb 2014
 *      Author: futyand
 */

#include <iostream>
#include <fstream>

#include "boost/math/tools/stats.hpp"

#include "REF_APP_Jitter.hxx"

#include "JitterProducer.hxx"

void JitterProducer::initialize(const ModuleParams& params) {

  boost::math::tools::stats<double> jitterStats,jitterRollStats,jitterPitchStats,jitterYawStats;

	m_jitterFilename = params.GetAsString("jitterFilename");
	m_jitterFileOffset = params.GetAsDouble("jitterFileOffset");
	m_jitterScale = params.GetAsDouble("jitterScale");

	m_jitterSize = 0;

    RefAppJitter * jitter_file;

    //If input jitter file is a user uploaded ascii file, convert it to FITS format
    if (m_jitterFilename.substr(m_jitterFilename.find_last_of(".")+1) == "txt") {

    	jitter_file = new RefAppJitter("jitter.fits", "CREATE");

		//user uploaded jitter file, to be read from execution directory
		ifstream jitter_in(m_jitterFilename, ios::in);
	    if (!jitter_in.good()) throw runtime_error("Error in JitterProducer::initialize: Error opening file "+m_jitterFilename);
	    string line;
	    double time,jitterX,jitterY,jitterZ;
	    int pitlFlag1,pitlFlag2;
	    while (getline(jitter_in, line)) {
	    	if(line[0] != '#') {
	    		stringstream(line) >> time >> pitlFlag1 >> pitlFlag2 >> jitterX >> jitterY >> jitterZ;
		    	if (time>48.*3600.) continue;
	    		if (time>=m_jitterFileOffset) {
	    			jitter_file->setCellTime(time-m_jitterFileOffset);
	    			jitter_file->setCellValidAocs(pitlFlag1);
	    			jitter_file->setCellValidScience(pitlFlag2!=255);
	    			jitter_file->setCellRoll(jitterX);
	    			jitter_file->setCellPitch(jitterY);
	    			jitter_file->setCellYaw(jitterZ);
	    			jitter_file->WriteRow();
	    		}
	    	}
	    }
		jitter_in.close();
		delete jitter_file;
		jitter_file = new RefAppJitter("jitter.fits","READONLY");

	} else {

		//Open the jitter fits file
		jitter_file = new RefAppJitter(string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_jitterFilename,"READONLY");

	}

	//Read in the jitter data
    int i_old = -1;
	while (jitter_file->ReadRow()) {
		if (m_jitterSize>kNJitter) break; //Truncate the time series at 48 hours in order not to exceed array size limit
		int i = jitter_file->getCellTime();
		if (i > i_old) {
			i_old = i;
			m_APE[i].m_X = jitter_file->getCellRoll()*m_jitterScale;
			m_APE[i].m_Y = jitter_file->getCellPitch()*m_jitterScale;
			m_APE[i].m_Z = jitter_file->getCellYaw()*m_jitterScale;
			m_validAocs[i] = jitter_file->getCellValidAocs();
			m_validScience[i] = jitter_file->getCellValidScience();
			//cout << i << " " << m_APE[i].m_X << " " << m_APE[i].m_Y << " " << m_APE[i].m_Z << endl;
			jitterRollStats.add(m_APE[i].m_X);
			jitterPitchStats.add(m_APE[i].m_Y);
			jitterYawStats.add(m_APE[i].m_Z);
			jitterStats.add(sqrt(m_APE[i].m_Y*m_APE[i].m_Y + m_APE[i].m_Z*m_APE[i].m_Z));
			m_jitterSize++;
		}
	}

	delete jitter_file;

	cout << "jitter file: " << m_jitterFilename << " : " << m_jitterSize << " lines read" << endl;
	cout << "jitter rms roll (arcseconds): " << jitterRollStats.rms() << endl;
	cout << "jitter rms pitch (arcseconds): " << jitterPitchStats.rms() << endl;
	cout << "jitter rms yaw (arcseconds): " << jitterYawStats.rms() << endl;
	cout << "jitter rms sqrt(pitch^2+yaw^2) (arcseconds): " << jitterStats.rms() << endl;

}

void JitterProducer::doBegin(Data* data, bool fullFrame) {

	data->setPreSubArrayTimeSteps();

	TimeConfiguration timeConf = data->getTimeConfiguration();

	for (int timeStep=-int(timeConf.getNumberOfPreSubArrayTimeSteps()); timeStep<int(timeConf.getNumberOfTimeSteps());  timeStep++) {

		SatelliteData * satData = data->getSatelliteData(timeStep,true);

		unsigned nSeconds = timeConf.getJitterSamplingsPerExposure();
		//If the timestep is the last one of the simulation, add two more seconds, to correspond to the end of the last second of the last exposure,
        //and to one second after that (required for SCI_RAW_Attitude file used for DFS input, filled in OrbitSimulator,
        //since the OBT time of the last CE is slightly delayed w.r.t. the end of the last exposure)
        if (timeStep == int(timeConf.getNumberOfTimeSteps()-1)) nSeconds += 2;

		unsigned secondsSinceStart = static_cast<unsigned>(lround(timeConf.getTimeSinceVisitStart(timeStep)));

		for (unsigned i=0; i<nSeconds; i++) {
			unsigned ijitter = secondsSinceStart%m_jitterSize + i;
			//cout << secondsSinceStart << " " << ijitter << " " << m_APE[ijitter].m_X << " " << m_APE[ijitter].m_Y << " " << m_APE[ijitter].m_Z << endl;
			satData->addJitterAPE(m_APE[ijitter],m_validAocs[ijitter],m_validScience[ijitter]);
		}

	}

}

