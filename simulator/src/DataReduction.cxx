/*
 * DataReduction.cxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

#include <fstream>

#include "boost/filesystem.hpp"
#include "boost/foreach.hpp"
#include "boost/algorithm/string/replace.hpp"

#include "BarycentricOffset.hxx"
#include "REF_APP_GainCorrection.hxx"

#include "data/include/Star.hxx"
#include "satellite/include/OrbitSimulator.hxx"
#include "telescope/include/PSFGenerator.hxx"
#include "Photometry.hxx"
#include "DataReduction.hxx"

void DataReduction::initialize(const ModuleParams & params) {

	m_writeIncidentLightCurve = params.GetAsBool("writeIncidentLightCurve");
	m_extractLightCurve = params.GetAsBool("extractLightCurve");
	m_imageDirectory = params.GetAsString("imageDirectory");
	m_generateNoiseCurve = params.GetAsBool("generateNoiseCurve");
	m_subtractFlatField = params.GetAsBool("subtractFlatField");

	m_outputDirectory = "light_curves";

	m_truthBarycentre = params.GetAsBool("truthBarycentre");

}

void DataReduction::doBegin(Data* data, bool fullFrame) {

//	data->getIdealLightCurveParams()->m_photonNoise = false;
//	data->getIdealLightCurveParams()->m_biasOffset = 553.8;
//	data->getIdealLightCurveParams()->m_readNoise = 7.2;

	//Initialize the photometric extraction
	m_photometry = new Photometry(data->getSubarrayDimensions(),data->getPhotometryParams());

	boost::mt19937 randomNumberEngine(3579);
    boost::poisson_distribution<int,double> poissonDistribution(1.);
    m_poissonNoiseGenerator = new RANDOM_POISSON(randomNumberEngine,poissonDistribution);

	boost::mt19937 randomNumberEngine2(3579+1);
    boost::normal_distribution<double> normalDistribution(data->getIdealLightCurveParams()->m_biasOffset,data->getIdealLightCurveParams()->m_readNoise);
    m_gaussianNoiseGenerator = new RANDOM_GAUSS(randomNumberEngine2,normalDistribution);

}

void DataReduction::doEnd(Data * data) const {

	if (m_writeIncidentLightCurve || m_extractLightCurve) {
		boost::filesystem::create_directory(data->getOutputDirectory() + "/" + m_outputDirectory);
	}

	if (m_writeIncidentLightCurve) {
		if (data->getFieldOfView()->getStars().size()>0) {
			writeIncidentLightCurve(data);
		} else {
			cout << "INFO: Incident light curve will not be generated (requires StarProducer to be run)" << endl;
		}
	}

	if (m_extractLightCurve) {

		//Read the full set of stacked images from files on disk and extract the flux from each image
		vector<Flux> extractedFlux = extractFluxVector(data);
		//for (unsigned i=0; i<extractedFlux.size(); i++) cout << setprecision(9) << extractedFlux[i] << endl;
		if (extractedFlux.size() != data->getTimeConfiguration().getNumberOfStackedImagesPITL()) {
			throw runtime_error("Error in DataReduction::process: size of extracted flux vector not consistent with time configuration");
		}
		writeOutputLightCurves(data,extractedFlux);

		//Uncomment the following line to append the extracted flux from the first image to an ascii file.
		//This can be used together with the extractedFlux.sh script to generate an array of extracted flux
		//values for different magnitudes and spectral types
		//writeExtractedFlux(data, extractedFlux[0].m_flux);

	}

}

vector<DataReduction::Flux> DataReduction::extractFluxVector(Data* data) const {

	//Write any remaining images in memory to disk
	data->resetFitsImageCube(Data::stacked);

	//Get the input directory path
	boost::filesystem::path imageDirectory;
	if (m_imageDirectory=="null" || m_imageDirectory=="") {
		imageDirectory = data->getOutputDirectory()+"/data";
	} else {
		imageDirectory = m_imageDirectory.c_str();
	}
	if (!boost::filesystem::exists(imageDirectory)) throw runtime_error("Error in DataReduction::readStackedImages: input directory does not exist");

	//For the case of flat field subtraction, set the flat field to be subtracted
	if (m_subtractFlatField) {
		Image * flatField = data->getFlatFieldToSubtract();
		if (flatField != nullptr) {
			m_photometry->setFlatField(flatField);
		} else {
			throw runtime_error("Error in DataReduction::extractFluxVector: Flat field subtraction requested, but flat field is undefined");
		}
	}

	vector<Flux> extractedFlux;
	string startTime="";

	//Loop over all files in the specified directory
	bool doublePrecision = false;
	int nFiles = 0;
	int nImages = 0;
	boost::filesystem::directory_iterator it(imageDirectory), eod;
	BOOST_FOREACH(boost::filesystem::path const &fs_path, std::make_pair(it, eod)) {
		if(is_regular_file(fs_path)) {

			//Check that the filename ends with .fits, contains "SCI_RAW_SubArray" or "SIM_RAW_DoublePrecisionSubArray"
			//and contains PR99 or PR90 (signifying stacked simulated images)
			string filename = fs_path.string();
			if(filename.substr(filename.find_last_of(".")+1) == "fits" &&
			   (filename.find("PR"+to_string(Data::kProgramType)) != string::npos || filename.find("PR"+to_string(data->getVisit().m_progType)) != string::npos)) {

				if (filename.find("SCI_RAW_SubArray") != string::npos) {
					extractFluxesFromCube<SciRawSubarray>(data,filename,nFiles,nImages,extractedFlux);
				} else if (filename.find("SIM_RAW_DoublePrecisionSubArray") != string::npos) {
					doublePrecision = true;
					extractFluxesFromCube<SimRawDoubleprecisionsubarray>(data,filename,nFiles,nImages,extractedFlux);
				}

			}

		}

	}

	if (nFiles==0) {
		throw runtime_error("Error in DataReduction::readStackedImages: no image files found");
	} else {
		if (doublePrecision) {
			setTimeConfigurationFromImageCube<SimRawDoubleprecisionsubarray>(data,nImages);
		} else {
			setTimeConfigurationFromImageCube<SciRawSubarray>(data,nImages);
		}
	}

	return extractedFlux;

}

template<class T> void DataReduction::extractFluxesFromCube(Data* data, string filename, int & nFiles, int & nImages, vector<Flux> & extractedFlux) const {

	cout << "Reading images from " << filename << endl;

	//Get the image cube and truth data
	T * imageCube = new T(filename,"READONLY");
	SciRawImagemetadata * metaData = new SciRawImagemetadata(filename,"READONLY");
	SimTruSubarray * truthData = nullptr;
	if (m_truthBarycentre) {
		boost::replace_last(filename, "data", "truth");
		boost::replace_first(filename, "SCI_RAW_SubArray", "SIM_TRU_SubArray");
		boost::replace_first(filename, "SIM_RAW_DoublePrecisionSubArray", "SIM_TRU_SubArray");
		truthData = new SimTruSubarray(filename,"READONLY");
	}

	//Extract the fluxes from all the images in the cube
	long cubeSize[3];
	imageCube->GetSize(cubeSize);

	//Loop over image cube layers
	for (int imageCount = 0; imageCount<cubeSize[2]; imageCount++) {

		//Fill the image array
		Image * image = new Image(cubeSize[0],cubeSize[1],imageCube->getKeyXWinoff(),imageCube->getKeyYWinoff());
		for (int j = 0; j<cubeSize[1]; j++) {
			for (int i = 0; i<cubeSize[0]; i++) {
				image->setPixelValue(image->getXOffset()+i,
						image->getYOffset()+j,
						((*imageCube)[imageCount])[j][i]);
			}
		}

		//Get the image mid time from the metadata
		if (!metaData->ReadRow()) throw runtime_error("Error reading metadata in DataReduction::extractFluxesFromCube");
		UTC midTime = metaData->getCellUtcTime();

		//Get the jitter AOCS and science validity flags from the truth metadata
		bool validAocs = true;
		bool validScience = true;
		double rollAngle = -1.;
		if (truthData != nullptr) {
			if (!truthData->ReadRow()) throw runtime_error("Error reading truth data in DataReduction::extractFluxesFromCube");
			validAocs = truthData->getCellValidAocs();
			validScience = truthData->getCellValidScience();
			rollAngle = truthData->getCellRollAngle();
		}

		//define the metadata status flag based on the jitter validity flags
		int status;
		if (!validAocs && !validScience) {
			status = 0;
		} else if (validAocs && validScience) {
			status = 1;
		} else if (validAocs && !validScience) {
			status = 2;
		} else {
			status = 3;
		}

		//Determine PSF barycentre using either truth information or photometry
		pair<double,double> barycentre = make_pair(0.,0.);
		if (m_truthBarycentre) {
			if (truthData != nullptr) {
				if (!isnan(truthData->getCellPsfMeanX()[0])) {
					Data::SubarrayDimensions subarray = data->getSubarrayDimensions();
					double targetLocationX_subarray = subarray.m_targetLocationX - image->getXOffset();
					double targetLocationY_subarray = subarray.m_targetLocationY - image->getYOffset();
					barycentre = make_pair(truthData->getCellPsfMeanX()[0]-targetLocationX_subarray,
										   truthData->getCellPsfMeanY()[0]-targetLocationY_subarray);
				}
			}
		} else {
			barycentre = m_photometry->getBarycentre(image);
		}

		//Get the centroid location
		Data::SubarrayDimensions subarray = data->getSubarrayDimensions();
		double centroidX = subarray.m_targetLocationX - subarray.m_xOffset + barycentre.first;
		double centroidY = subarray.m_targetLocationY - subarray.m_yOffset + barycentre.second;

		//Extract the flux from the image array
    	double flux = m_photometry->extractFlux(image,barycentre);
    	delete image;
    	extractedFlux.push_back(Flux(flux,midTime,status,rollAngle,centroidX,centroidY));

	}

	//Store the first image cube in the data for access to header keywords etc
	if (nFiles==0) {
		data->setFitsImageCube(imageCube, nullptr, nullptr);
	} else {
		delete imageCube;
	}
	if (m_truthBarycentre) delete truthData;

	nFiles++;
	nImages += cubeSize[2];

}

template <class T> void DataReduction::setTimeConfigurationFromImageCube(Data* data, unsigned nImages) const {

	T * imageCube = nullptr;
	data->getFitsImageCube(&imageCube);
	if (imageCube == nullptr) throw runtime_error("Error in DataReduction::setTimeConfigurationFromImageCube: image cube pointer is null");

	//Get the start time of the simulation, subtracting half of the image duration from the TStrtU keyword
	//since the latter corresponds to the mid time of the first image
	//Also subtract 30s because the TimeConfiguration constructor will add a 30s delay to the specified start time
	UTC startTime_utc = imageCube->getKeyTStrtU() - DeltaTime(static_cast<double>(imageCube->getKeyTexptime()/2.)) - DeltaTime(30);

	//Convert from UTC to boost::posix_time::ptime
	string startTime_str = *(&startTime_utc);
	int delimiter = startTime_str.find_first_of("T");
	string date = startTime_str.substr(0,delimiter);
	string time = startTime_str.substr(delimiter+1);
	boost::posix_time::ptime startTime = boost::posix_time::time_from_string(date+" "+time);

	unsigned exposuresPerStack = imageCube->getKeyNexp();
	double exposureTimeAsDouble = imageCube->getKeyExptime();

	TimeConfiguration timeConf(startTime,nImages,exposuresPerStack,exposureTimeAsDouble,0,false);

	data->setTimeConfiguration(timeConf);

}

void DataReduction::writeIncidentLightCurve(Data * data) const {

	if (data->getFieldOfView()->getStars().size() == 0) throw runtime_error("Error in DataReduction::writeIncidentLightCurve: No stars in FOV.");

	//Open output fits files
	SciCorLightcurve * lightCurve = new SciCorLightcurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+"incidentLightCurve.fits", "CREATE");
	SciCorLightcurve * normalizedLightCurve = new SciCorLightcurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+"incidentNormalizedLightCurve.fits", "CREATE");
	SciCorLightcurve * idealLightCurve = new SciCorLightcurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+"incidentLightCurveWithPhotonNoise.fits", "CREATE");

	//Calculate the mid-time of the image
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double exposureTime = timeConf.getRepetitionPeriod();
	UTC startTime = to_iso_extended_string(timeConf.getStartTime());
	DeltaTime halfExposure = exposureTime/2.;
	UTC lastImageMidTime = startTime + halfExposure;

	//Calculate the integrated throughput*QE
	double effectiveTemperature = 5660.; //Default to effective temperature for G5 star if there are no stars in the FOV
	if (data->getFieldOfView()->getStars().size()>0) effectiveTemperature = (data->getFieldOfView()->getStars()[0])->getEffectiveTemperature();
	double integratedThroughputQE = data->getWavelengthDependence()->getIntegratedThroughputQE(effectiveTemperature,SatelliteData::kDefaultCcdTemperature);

	//Extract the value of the nominal gain
//	RefAppGaincorrection * gainCorrection_file = new RefAppGaincorrection(string(getenv("CHEOPS_SW"))+"/resources/"+data->getGainFilename(),"READONLY");
//	double nominalGain = gainCorrection_file->getKeyGainNom();
//	delete gainCorrection_file;

	//Count the number of pixels within the photometric extraction aperture
	int nPixInAperture = 0;
	for (int ix=412; ix<612; ix++) {
		for (int iy=412; iy<612; iy++) {
			double x = double(ix) - 512.;
			double y = double(iy) - 512.;
			double radius = sqrt(x*x + y*y);
			if (radius < data->getPhotometryParams().m_radius_psf) nPixInAperture++;
		}
	}

	//Write total incident flux target star in electrons (applying global throughput) as a function of time
	//and write ideal light curve, which is the incident light curve with photon noise applied
	double idealLightCurveFlux = 0.;
    const Star * targetStar = data->getFieldOfView()->getStars()[0];
    for (unsigned iTimeStep=0; iTimeStep<timeConf.getNumberOfTimeSteps();  iTimeStep++) {

    	DeltaTime exposureStart = timeConf.getTimeSinceStart(iTimeStep);
    	UTC exposureMidTime = startTime + exposureStart + halfExposure;
    	lastImageMidTime = exposureMidTime;

    	//Calculate the number of incident photons before noise
    	StarData * targetStarTimeData = targetStar->getTimeSeries()[iTimeStep];
    	double fluxFactor = targetStarTimeData->getCombinedFluxFactor();
    	double flux = targetStar->getMeanPhotonFlux() * fluxFactor * PSFGenerator::kTelescopeArea * data->getTimeConfiguration().getExposureTimeAsDouble();
    	//cout << "Target star flux (" << timeConf.getTimeSinceStart(iTimeStep)/3600. << " hours) = " << flux << endl;

    	//Apply global throughput
    	flux *= integratedThroughputQE;

    	//Apply photon noise
    	double fluxWithNoise = flux;
    	//fluxWithNoise *= 0.9644020037405; //Take into account fraction of measured PSF within aperture radius 30 pixels
    	//if (data->getIdealLightCurveParams()->m_photonNoise) {
    		m_poissonNoiseGenerator->distribution().param(boost::poisson_distribution<int,double>::param_type(fluxWithNoise));
    		fluxWithNoise = (*m_poissonNoiseGenerator)();
    	//}

    	//Apply the electronic gain
//    	flux *= nominalGain;
//    	fluxWithNoise *= nominalGain;

    	//Apply bias offset and readout noise to each pixel within the aperture, adding to the flux
//    	if (data->getIdealLightCurveParams()->m_biasOffset>0. || data->getIdealLightCurveParams()->m_readNoise>0.) {
//    		for (int i=0; i< nPixInAperture; i++) fluxWithNoise += (*m_gaussianNoiseGenerator)();
//    	}

    	//Add points to the incident light curve for each unstacked image
    	setLightCurveData(lightCurve,flux,exposureMidTime,0);
    	setLightCurveData(normalizedLightCurve,fluxFactor,exposureMidTime,0);

    	//Add points to light curves with the cadence of the stacking
    	idealLightCurveFlux += fluxWithNoise;
    	if ((iTimeStep+1)%data->getTimeConfiguration().getExposuresPerStack() == 0) {
    		setLightCurveData(idealLightCurve,idealLightCurveFlux,exposureMidTime,0);
    		idealLightCurveFlux = 0.;
    	}

    }

	setLightCurveHeaderKeywords(data,lightCurve,lastImageMidTime);
	setLightCurveHeaderKeywords(data,normalizedLightCurve,lastImageMidTime);
	setLightCurveHeaderKeywords(data,idealLightCurve,lastImageMidTime);

    delete lightCurve;
    delete normalizedLightCurve;
    delete idealLightCurve;

}

void DataReduction::writeOutputLightCurves(Data* data, vector<Flux> & extractedFlux) const {

	//Get the time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double repetitionPeriod = timeConf.getRepetitionPeriod();
	double exposureTimeAsDouble = timeConf.getExposureTimeAsDouble();
    unsigned exposuresPerStack = timeConf.getExposuresPerStack();

	//Obtain flux statistics outside the transit, needed to obtain the normalized light curve.
	boost::math::tools::stats<double> fluxStatsOutsideTransit;
    for (unsigned iTimeStep=0; iTimeStep<extractedFlux.size(); iTimeStep++) {
    	//Work out if current timestep is outside the transit using the Incident light curve, taking into account
    	//the fact that the Incident light curve has one time step per exposure instead of one time step per stacked image
        DeltaTime timeSinceStart = extractedFlux[iTimeStep].m_midTime - UTC(to_iso_extended_string(timeConf.getStartTime()));
    	int inputTimeStep = static_cast<int>(floor(timeSinceStart.getSeconds()/repetitionPeriod));
    	double transitFluxFactor = 1.;
    	if (data->getFieldOfView()->getStars().size()>0) {
    		transitFluxFactor = data->getFieldOfView()->getStars()[0]->getTimeSeries()[inputTimeStep]->getTransitFluxFactor();
    	}
    	if (!(transitFluxFactor<1.)) {
    		fluxStatsOutsideTransit.add(extractedFlux[iTimeStep].m_flux);
    	}
    }

	//Open output fits files
	SciCorLightcurve * lightCurve = new SciCorLightcurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+"extractedLightCurve.fits", "CREATE");
	SciCorLightcurve * normalizedLightCurve;
	SciCorLightcurve * OMinusCLightCurve;
	if (fluxStatsOutsideTransit.count()>0) {
		normalizedLightCurve = new SciCorLightcurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+"extractedNormalizedLightCurve.fits", "CREATE");
		OMinusCLightCurve = new SciCorLightcurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+"extractedOMinusCLightCurve.fits", "CREATE");
	}

	//Set header information
	UTC lastImageMidTime = extractedFlux[extractedFlux.size()-1].m_midTime;
	setLightCurveHeaderKeywords(data,lightCurve,lastImageMidTime,true);
	if (fluxStatsOutsideTransit.count()>0) {
		setLightCurveHeaderKeywords(data,normalizedLightCurve,lastImageMidTime,true);
		setLightCurveHeaderKeywords(data,OMinusCLightCurve,lastImageMidTime,true);
	}

	//Get the intended target location
	Data::SubarrayDimensions subarray = data->getSubarrayDimensions();
	double locationX = subarray.m_targetLocationX - subarray.m_xOffset;
	double locationY = subarray.m_targetLocationY - subarray.m_yOffset;

    //Fill the information into the light curves
	boost::math::tools::stats<double> OMinusCfluxStats;
	vector<double> transitSubtractedNormalizedFlux;
    for (unsigned iTimeStep=0; iTimeStep<extractedFlux.size(); iTimeStep++) {

    	//Raw extracted light curve
    	setLightCurveData(lightCurve,extractedFlux[iTimeStep].m_flux,extractedFlux[iTimeStep].m_midTime,extractedFlux[iTimeStep].m_status,
    			extractedFlux[iTimeStep].m_rollAngle,extractedFlux[iTimeStep].m_centroidX,extractedFlux[iTimeStep].m_centroidY,locationX,locationY);

    	if (fluxStatsOutsideTransit.count()>0) {

    		//Normalized light curve
    		double normalizedFlux = extractedFlux[iTimeStep].m_flux/fluxStatsOutsideTransit.mean();
        	setLightCurveData(normalizedLightCurve,normalizedFlux,extractedFlux[iTimeStep].m_midTime,extractedFlux[iTimeStep].m_status,
        			extractedFlux[iTimeStep].m_rollAngle,extractedFlux[iTimeStep].m_centroidX,extractedFlux[iTimeStep].m_centroidY,locationX,locationY);

        	//Observed minus calculated light curve
    		double transitFluxFactor = 1.;
    		if (data->getFieldOfView()->getStars().size()>0) {
            	int inputTimeStep = iTimeStep*exposuresPerStack + exposuresPerStack/2;
            	transitFluxFactor = data->getFieldOfView()->getStars()[0]->getTimeSeries()[inputTimeStep]->getTransitFluxFactor();
    		}
    		double OMinusCFlux = extractedFlux[iTimeStep].m_flux/fluxStatsOutsideTransit.mean() - transitFluxFactor;
        	setLightCurveData(OMinusCLightCurve,OMinusCFlux,extractedFlux[iTimeStep].m_midTime,extractedFlux[iTimeStep].m_status,
        			extractedFlux[iTimeStep].m_rollAngle,extractedFlux[iTimeStep].m_centroidX,extractedFlux[iTimeStep].m_centroidY,locationX,locationY);
    		OMinusCfluxStats.add(OMinusCFlux);
    		transitSubtractedNormalizedFlux.push_back(OMinusCFlux+1.);

    	}

    }

    //Write light curves to disk
    delete lightCurve;
	if (fluxStatsOutsideTransit.count()>0) {
		delete normalizedLightCurve;
		delete OMinusCLightCurve;
	}

	//Print summary information to the log file, and generate noise curves
    if (fluxStatsOutsideTransit.count()>0) {
    	cout << "Mean flux in extracted light curve: " << setprecision(10) << fluxStatsOutsideTransit.mean() << endl;
    	cout << "Standard deviation in Observed minus Calculated light curve (ppm):   " << pow(10.,6)*sqrt(OMinusCfluxStats.variance1()) << endl;
    	cout << "Standard deviation in extracted light curve excluding transit (ppm): " << pow(10.,6)*sqrt(fluxStatsOutsideTransit.variance1())/fluxStatsOutsideTransit.mean() << endl;
        if (m_generateNoiseCurve) {
        	generateNoiseCurve(data, transitSubtractedNormalizedFlux, extractedFlux[extractedFlux.size()-1].m_midTime, "noiseCurve.fits");
        }
       	double photonFlux = data->getFieldOfView()->getStars().size()>0 ? data->getFieldOfView()->getStars()[0]->getMeanPhotonFlux() : 0.;
        double incidentFlux = photonFlux * PSFGenerator::kTelescopeArea * exposureTimeAsDouble*exposuresPerStack;
        if (incidentFlux>0.) cout << "Mean flux outside transit / incident flux: " << fluxStatsOutsideTransit.mean()/incidentFlux << endl;
    } else {
    	cout << "A normalized light curve, O-C light curve and noise curve have not been generated because no part of the simulation is outside the transit" << endl;
    }

}

void DataReduction::setLightCurveHeaderKeywords(Data * data, SciCorLightcurve * lightCurve, UTC lastImageMidTime, bool photometry) const {

	//Copy fits headers from image cube
	SciRawSubarray * imageCube = nullptr;
	data->getFitsImageCube(&imageCube);
	if (imageCube != nullptr) {
		setLightCurveHeaderKeywordsFromImageCube<SciRawSubarray>(lightCurve,imageCube);
	} else {
		SimRawDoubleprecisionsubarray * imageCube = nullptr;
		data->getFitsImageCube(&imageCube);
		if (imageCube != nullptr) {
			setLightCurveHeaderKeywordsFromImageCube<SimRawDoubleprecisionsubarray>(lightCurve,imageCube);
		} else {
			//No image cube available. This can be the case if only the incident light curve is being written
			setLightCurveHeaderKeywordsWithoutImageCube(data,lightCurve);
		}
	}

	//Get the mid time of the first image
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double imageDuration = (timeConf.getExposuresPerStack()-1)*timeConf.getRepetitionPeriod() + timeConf.getExposureTimeAsDouble();
	UTC startTime = UTC(to_iso_extended_string(timeConf.getStartTime()));
	UTC firstImageMidTime = startTime + DeltaTime(imageDuration/2.);
	UTC endTime = startTime + DeltaTime(timeConf.getDuration());

	lightCurve->setKeyVStrtU(startTime);
	lightCurve->setKeyVStopU(endTime);
	lightCurve->setKeyVStrtM(startTime.getMjd());
	lightCurve->setKeyVStopM(endTime.getMjd());
	lightCurve->setKeyTStrtU(firstImageMidTime);
	lightCurve->setKeyTStopU(lastImageMidTime);
	lightCurve->setKeyTStrtM(firstImageMidTime.getMjd());
	lightCurve->setKeyTStopM(lastImageMidTime.getMjd());
	BarycentricOffset offset(lightCurve->getKeyRaTarg(),lightCurve->getKeyDecTarg());
	lightCurve->setKeyTStrtB(firstImageMidTime.getMjd() + offset);
	lightCurve->setKeyTStopB(lastImageMidTime.getMjd() + offset);
	if (data->hasSatelliteData()) {
		SatelliteData * satData_firstExp = data->getSatelliteData(0,false,false);
		SatelliteData * satData_lastExp = data->getSatelliteData(timeConf.getNumberOfTimeSteps()-1,false,false);
		lightCurve->setKeyBSunA(satData_firstExp->getSunAngle());
		lightCurve->setKeyBMoonA(satData_firstExp->getMoonAngle());
		lightCurve->setKeyBEartA(satData_firstExp->getEarthLimbAngle());
		lightCurve->setKeyESunA(satData_lastExp->getSunAngle());
		lightCurve->setKeyEMoonA(satData_lastExp->getMoonAngle());
		lightCurve->setKeyEEartA(satData_lastExp->getEarthLimbAngle());
	}
	if (photometry) {
		lightCurve->setKeyApRadi(data->getPhotometryParams().m_radius_psf);
		lightCurve->setKeyApType("r33");
		lightCurve->setKeyDataname("Rdef");
	} else {
		lightCurve->setKeyApType("none");
	}

}

template <class T> void DataReduction::setLightCurveHeaderKeywordsFromImageCube(SciCorLightcurve * lightCurve, T * imageCube) const {

	lightCurve->setKeyProcChn(imageCube->getKeyProcChn());
	lightCurve->setKeyProcNum(imageCube->getKeyProcNum());
	lightCurve->setKeyArchRev(imageCube->getKeyArchRev());
	lightCurve->setKeyPiName(imageCube->getKeyPiName());
	lightCurve->setKeyPiUid(imageCube->getKeyPiUid());
	lightCurve->setKeyObsCat(imageCube->getKeyObsCat());
	lightCurve->setKeyProgtype(imageCube->getKeyProgtype());
	lightCurve->setKeyProgId(imageCube->getKeyProgId());
	lightCurve->setKeyReqId(imageCube->getKeyReqId());
	lightCurve->setKeyVisitctr(imageCube->getKeyVisitctr());
	lightCurve->setKeyObsid(imageCube->getKeyObsid());
	lightCurve->setKeyPrpVst1(imageCube->getKeyPrpVst1());
	lightCurve->setKeyPrpVstn(imageCube->getKeyPrpVstn());
	lightCurve->setKeyNexp(imageCube->getKeyNexp());
	lightCurve->setKeyExptime(imageCube->getKeyExptime());
	lightCurve->setKeyTexptime(imageCube->getKeyTexptime());
	lightCurve->setKeyRaTarg(imageCube->getKeyRaTarg());
	lightCurve->setKeyDecTarg(imageCube->getKeyDecTarg());
	lightCurve->setKeyTargname(imageCube->getKeyTargname());
	lightCurve->setKeySpectype(imageCube->getKeySpectype());
	lightCurve->setKeyTEff(imageCube->getKeyTEff());
	lightCurve->setKeyMagG(imageCube->getKeyMagG());
	lightCurve->setKeyMagGerr(imageCube->getKeyMagGerr());
	lightCurve->setKeyMagChps(imageCube->getKeyMagChps());
	lightCurve->setKeyMagCerr(imageCube->getKeyMagCerr());

}

void DataReduction::setLightCurveHeaderKeywordsWithoutImageCube(const Data* data, SciCorLightcurve* lightCurve) const {

	TimeConfiguration timeConf = data->getTimeConfiguration();
	unsigned exposureTimeAsDouble = timeConf.getExposureTimeAsDouble();
	unsigned exposuresPerStack = timeConf.getExposuresPerStack();

	string specType = "G5"; //Default to spectral type to G5 if there are no stars in the FOV
	double Gmag = 9.; //Default to 9th magnitude if there are no stars in the FOV
	double cheopsMag = 9.; //Default to 9th magnitude if there are no stars in the FOV
	double GmagErr = 0.;
	double cheopsMagErr = 0.;
	double TEff = 0;
	if (data->getFieldOfView()->getStars().size()>0) {
		Star * targetStar = data->getFieldOfView()->getStars()[0];
		specType = targetStar->getSpectralTypeString();
		TEff = targetStar->getEffectiveTemperature();
		Gmag = targetStar->getGaiaMagnitude();
		cheopsMag = targetStar->getCheopsMagnitude();
		GmagErr = targetStar->getGaiaMagnitudeError();
		cheopsMagErr = targetStar->getCheopsMagnitudeError();
	}

	Data::Visit visit = data->getVisit();

	lightCurve->setKeyProcChn("CHEOPSim");
	lightCurve->setKeyProcNum(visit.m_versionNum);
	lightCurve->setKeyArchRev(0);
	lightCurve->setKeyPiName(visit.m_piName);
	lightCurve->setKeyPiUid(visit.m_piUid);
	lightCurve->setKeyObsCat(visit.m_obsCat);
	lightCurve->setKeyProgtype(visit.m_progType);
	lightCurve->setKeyProgId(visit.m_progId);
	lightCurve->setKeyReqId(visit.m_reqId);
	lightCurve->setKeyVisitctr(visit.m_visitCtr);
	lightCurve->setKeyObsid(visit.m_obsId);
	lightCurve->setKeyPrpVst1(visit.m_prpFirst);
	lightCurve->setKeyPrpVstn(visit.m_prpLast);
	lightCurve->setKeyNexp(exposuresPerStack);
	lightCurve->setKeyExptime(exposureTimeAsDouble);
	lightCurve->setKeyTexptime(exposureTimeAsDouble*exposuresPerStack);
	lightCurve->setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	lightCurve->setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
	lightCurve->setKeyTargname(visit.m_targName);
	lightCurve->setKeySpectype(specType);
	lightCurve->setKeyTEff(TEff);
	lightCurve->setKeyMagG(Gmag);
	lightCurve->setKeyMagChps(cheopsMag);
	lightCurve->setKeyMagGerr(GmagErr);
	lightCurve->setKeyMagCerr(cheopsMagErr);

}

void DataReduction::setLightCurveData(SciCorLightcurve * lightCurve, double flux, UTC midTime, int status, double rollAngle,
		double centroidX, double centroidY, double locationX, double locationY) const {

	lightCurve->setCellFlux(flux);
	lightCurve->setCellFluxerr(0.);
	lightCurve->setCellUtcTime(midTime);
	lightCurve->setCellMjdTime(midTime.getMjd());
	BarycentricOffset offset(lightCurve->getKeyRaTarg(),lightCurve->getKeyDecTarg());
	lightCurve->setCellBjdTime(midTime.getMjd() + offset);
	lightCurve->setCellStatus(status);
	if (centroidX>0. && centroidY>0. && locationX>0. && locationY>0.) {
		lightCurve->setCellCentroidX(centroidX);
		lightCurve->setCellCentroidY(centroidY);
		lightCurve->setCellLocationX(locationX);
		lightCurve->setCellLocationY(locationY);
	}
	if (rollAngle>0.) lightCurve->setCellRollAngle(rollAngle);
	lightCurve->WriteRow();

}

void DataReduction::generateNoiseCurve(const Data * data, const vector<double> flux, UTC lastImageMidTime, string filename) const {

	SimAnaNoisecurve noiseCurve(data->getOutputDirectory() + "/" + m_outputDirectory+"/"+filename, "CREATE");

	//Copy fits headers from image cube
	SciRawSubarray * imageCube = nullptr;
	data->getFitsImageCube(&imageCube);
	if (imageCube != nullptr) {
		setNoiseCurveHeaderKeywordsFromImageCube<SciRawSubarray>(noiseCurve,imageCube);
	} else {
		SimRawDoubleprecisionsubarray * imageCube = nullptr;
		data->getFitsImageCube(&imageCube);
		if (imageCube != nullptr) {
			setNoiseCurveHeaderKeywordsFromImageCube<SimRawDoubleprecisionsubarray>(noiseCurve,imageCube);
		} else {
			throw runtime_error("Error in DataReduction::generateNoiseCurve: image cube pointer is null");
		}
	}

	//Get firstimage mid time
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double imageDuration = (timeConf.getExposuresPerStack()-1)*timeConf.getRepetitionPeriod() + timeConf.getExposureTimeAsDouble();
	UTC startTime = to_iso_extended_string(timeConf.getStartTime());
	UTC firstImageMidTime = startTime + DeltaTime(imageDuration/2.);
	UTC endTime = startTime + DeltaTime(timeConf.getDuration());

	noiseCurve.setKeyVStrtU(startTime);
	noiseCurve.setKeyVStopU(endTime);
	noiseCurve.setKeyVStrtM(startTime.getMjd());
	noiseCurve.setKeyVStopM(endTime.getMjd());
	noiseCurve.setKeyTStrtU(firstImageMidTime);
	noiseCurve.setKeyTStopU(lastImageMidTime);
	noiseCurve.setKeyTStrtM(firstImageMidTime.getMjd());
	noiseCurve.setKeyTStopM(lastImageMidTime.getMjd());
	BarycentricOffset offset(noiseCurve.getKeyRaTarg(),noiseCurve.getKeyDecTarg());
	noiseCurve.setKeyTStrtB(firstImageMidTime.getMjd() + offset);
	noiseCurve.setKeyTStopB(lastImageMidTime.getMjd() + offset);

	for (unsigned nGroup=1; nGroup<flux.size()/2;  nGroup++) {

		vector<double> groupedFlux;
		double meanFluxInGroup = 0;
		unsigned indexWithinGroup = 0;
		for (unsigned i=0; i<flux.size();  i++) {

			meanFluxInGroup += flux[i];
			indexWithinGroup++;
			if (indexWithinGroup == nGroup) {
				groupedFlux.push_back(meanFluxInGroup/nGroup);
				meanFluxInGroup = 0.;
				indexWithinGroup = 0;
			}
		}

		boost::math::tools::stats<double> fluxStats;
		for (unsigned i=0; i<groupedFlux.size();  i++) fluxStats.add(groupedFlux[i]);
		double noise = 0.;
		if (fluxStats.mean()>0.) noise = pow(10.,6)*sqrt(fluxStats.variance1())/fluxStats.mean();

    	noiseCurve.setCellNoise(noise);
    	noiseCurve.setCellTimeBin(double(nGroup)*noiseCurve.getKeyTexptime()/3600.);
    	noiseCurve.WriteRow();

		//cout << nGroup << " " << groupedFlux.size() << " " << noise << endl;
	}

}

template <class T> void DataReduction::setNoiseCurveHeaderKeywordsFromImageCube(SimAnaNoisecurve & noiseCurve, T * imageCube) const {

	noiseCurve.setKeyProcNum(imageCube->getKeyProcNum());
	noiseCurve.setKeyArchRev(imageCube->getKeyArchRev());
	noiseCurve.setKeyPiName(imageCube->getKeyPiName());
	noiseCurve.setKeyPiUid(imageCube->getKeyPiUid());
	noiseCurve.setKeyObsCat(imageCube->getKeyObsCat());
	noiseCurve.setKeyProgtype(imageCube->getKeyProgtype());
	noiseCurve.setKeyProgId(imageCube->getKeyProgId());
	noiseCurve.setKeyReqId(imageCube->getKeyReqId());
	noiseCurve.setKeyVisitctr(imageCube->getKeyVisitctr());
	noiseCurve.setKeyObsid(imageCube->getKeyObsid());
	noiseCurve.setKeyPrpVst1(imageCube->getKeyPrpVst1());
	noiseCurve.setKeyPrpVstn(imageCube->getKeyPrpVstn());
	noiseCurve.setKeyNexp(imageCube->getKeyNexp());
	noiseCurve.setKeyExptime(imageCube->getKeyExptime());
	noiseCurve.setKeyTexptime(imageCube->getKeyTexptime());
	noiseCurve.setKeyRaTarg(imageCube->getKeyRaTarg());
	noiseCurve.setKeyDecTarg(imageCube->getKeyDecTarg());
	noiseCurve.setKeyTargname(imageCube->getKeyTargname());
	noiseCurve.setKeySpectype(imageCube->getKeySpectype());
	noiseCurve.setKeyTEff(imageCube->getKeyTEff());
	noiseCurve.setKeyMagG(imageCube->getKeyMagG());
	noiseCurve.setKeyMagGerr(imageCube->getKeyMagGerr());
	noiseCurve.setKeyMagChps(imageCube->getKeyMagChps());
	noiseCurve.setKeyMagCerr(imageCube->getKeyMagCerr());

}

void DataReduction::writeExtractedFlux(Data * data, double flux) const {

	unsigned exposureTimeAsDouble = data->getTimeConfiguration().getExposureTimeAsDouble();
	flux/=exposureTimeAsDouble;
	cout << "extracted flux: " << flux << " ADU/s" << endl;

	double magnitude = data->getFieldOfView()->getStars()[0]->getVbandMagnitude();
	Star::spectralType spectralType = data->getFieldOfView()->getStars()[0]->getSpectralType();

	ofstream extractedFlux_file("extractedFlux.txt", ios::app);
	if (spectralType==Star::G0) extractedFlux_file << magnitude;
	extractedFlux_file << "	" << setprecision(10) << flux;
	if (spectralType==Star::M9) extractedFlux_file << endl;
	extractedFlux_file.close();

}
