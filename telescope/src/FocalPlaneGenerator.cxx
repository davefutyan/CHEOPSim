/*
 * FocalPlaneGenerator.cxx
 *
 *  Created on: 5 Feb 2014
 *      Author: futyand
 */

#include <fstream>

#include "boost/algorithm/string.hpp"
#include "fitsio.h"

#include "FocalPlaneGenerator.hxx"

void FocalPlaneGenerator::initialize(const ModuleParams& params) {

	m_subArrayXDim = params.GetAsInt("subArrayXDim");
	m_subArrayYDim = params.GetAsInt("subArrayYDim");
	m_subArrayXOffset = params.GetAsInt("subArrayXOffset");
	m_subArrayYOffset = params.GetAsInt("subArrayYOffset");
	m_targetLocationX = params.GetAsDouble("targetLocationX");
	m_targetLocationY = params.GetAsDouble("targetLocationY");
	string targetLocationListX_str = params.GetAsString("targetLocationListX");
	string targetLocationListY_str = params.GetAsString("targetLocationListY");

	if (targetLocationListX_str!="null" && targetLocationListX_str!="") {

		vector<string> substrings_x,substrings_y;
		boost::split(substrings_x,targetLocationListX_str,boost::is_any_of(","));
		boost::split(substrings_y,targetLocationListY_str,boost::is_any_of(","));

		if (substrings_x.size() != substrings_y.size()) {
			throw runtime_error("Error in FocalPlaneGenerator::initialize: inconsistent array sizes for targetLocationListX and targetLocationListY");
		}

		for (unsigned i=0; i<substrings_x.size(); i++) {
			m_targetLocationListX.push_back(boost::lexical_cast<double>(substrings_x[i]));
			m_targetLocationListY.push_back(boost::lexical_cast<double>(substrings_y[i]));
		}

	}

	m_plateScale = 1.;

}

void FocalPlaneGenerator::doBegin(Data* data, bool fullFrame) {

	data->setSubarrayDimensions(m_subArrayXDim,m_subArrayYDim,m_subArrayXOffset,m_subArrayYOffset,m_targetLocationX,m_targetLocationY);

	TimeConfiguration timeConf = data->getTimeConfiguration();

	if (m_subArrayXDim == Image::kXDim && m_subArrayYDim==Image::kYDim) {
		data->getMpsPreVisits()->setCellNFullframeExp(timeConf.getNumberOfStackedImages() * timeConf.getExposuresPerStack() + (timeConf.doFullFrame()?1:0));
		data->getMpsPreVisits()->setCellNWindowframeExp(0);
	}

	bool omitStackedImage = false;
	unsigned exposuresPerStack = timeConf.getExposuresPerStack();
	unsigned numberOfValidStackedImages = timeConf.getNumberOfStackedImages();

	for (unsigned timeStep=0; timeStep<timeConf.getNumberOfTimeSteps();  timeStep++) {

		SatelliteData * satData = data->getSatelliteData(timeStep,(!data->moduleIsActive("JitterProducer")&&!data->moduleIsActive("OrbitSimulator")&&!data->moduleIsActive("StrayLightGenerator")));

		//Check for Earth occultation, stray light or SAA in any exposure for the current stacked image
		if ((data->omitSAA() && satData->getSAAFlag()) ||
			(data->omitEarthOccultation() && satData->getEarthOccultationFlag()) ||
			(data->omitStrayLight() && satData->getStrayLightFlag())) omitStackedImage = true;

		if ((timeStep+1)%exposuresPerStack == 0) { //this timestep corresponds to the start of a new stacked image
			//Identify stacked images to be omitted due to option to omit first N images, or option to only write out every Nth image
			if (((timeStep+1)/exposuresPerStack <= data->getNumberOfCEsToOmitAtStart()) ||
				((timeStep+1)/exposuresPerStack%data->getOnlyWriteEveryNthCE() != 0)) omitStackedImage = true;
			if (omitStackedImage) {
				numberOfValidStackedImages--; //remove this image from the number of valid stacked images
				//Flag all time steps corresponding to the stacked image as discarded
				for (unsigned i=timeStep+1-exposuresPerStack; i<timeStep+1; i++) {
					data->getSatelliteData(i,false,false)->setDiscardFlag(true);
				}
			}
			omitStackedImage = false; //reset ready for the next stacked image
		}

	}

	data->setNumberOfStackedImagesPITL(numberOfValidStackedImages);

	if (m_targetLocationListX.size() != 0 && m_targetLocationListX.size() <  timeConf.getNumberOfTimeSteps()) {
		throw runtime_error("Error in FocalPlaneGenerator::doBegin: target location list is non-zero but shorter than the number of exposures");
	}

	string pixelScaleFilename = "CH_TU2016-01-01T00-00-00_REF_APP_PixelScale_V0101.fits";
	RefAppPixelscale * pixelScale = new RefAppPixelscale(string(getenv("CHEOPS_SW"))+"/resources/"+pixelScaleFilename,"READONLY");
	pixelScale->ReadRow();
	m_plateScale = pixelScale->getCellPixelScale();
	cout << "Plate scale = " << m_plateScale << endl;

	ofstream referenceFilesList("reference_files.txt", ios::app);
	referenceFilesList << pixelScaleFilename << endl;
	referenceFilesList.close();

}

void FocalPlaneGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	//If a list of different target locations has been provided for each time step (e.g. for M and C data), set the target location for the current time step
	if (m_targetLocationListX.size() != 0) {
		data->setTargetLocation(m_targetLocationListX[timeStep],m_targetLocationListY[timeStep]);
		cout << "Target location: (" << data->getSubarrayDimensions().m_targetLocationX << "," << data->getSubarrayDimensions().m_targetLocationY << ")" << endl;
	}

	//Initialize the images
	Image * image;
	if (fullFrame) {
		//Full frame image
		image = new Image(Image::kXDim,Image::kYDim,0,0);
	} else {
		//Sub-frame image with default dimensions
		image = new Image(m_subArrayXDim,m_subArrayYDim,m_subArrayXOffset,m_subArrayYOffset);
	}
	data->addImage(image);

	unsigned nJitterPerExposure = data->getTimeConfiguration().getJitterSamplingsPerExposure();

	//Get the vector of roll angles and the vector of pointing directions (filled if OrbitSimulator has been run)
	vector<double> rollAngles = data->getSatelliteData(timeStep)->getRollAngles();
	vector<SkyPosition> pointingDirections = data->getSatelliteData(timeStep)->getPointingDirections();
	if (rollAngles.size() != pointingDirections.size()) {
		throw runtime_error("Error in FocalPlaneGenerator::process: size of roll angle vector does not match size of point directions vector");
	} else if (rollAngles.size() != 0 && rollAngles.size() != nJitterPerExposure) {
		throw runtime_error("Error in FocalPlaneGenerator::process: roll angle vector must have size equal to exposure time or be empty");
	}

	//Get the jitter APEs if available (i.e. if JitterProducer has been run)
	vector<SatelliteData::APE> jitterAPEs = data->getSatelliteData(timeStep)->getJitterAPEs();
	//Size of APE vector is 2 seconds longer than the exposure duration for the last exposure to allow extra rows to be written to SCI_RAW_Attitude (needed for DFS)
	if (jitterAPEs.size() != 0 && jitterAPEs.size() != nJitterPerExposure && jitterAPEs.size() != nJitterPerExposure+2) {
		throw runtime_error("Error in FocalPlaneGenerator::process: jitter APE vector must have size equal to exposure time or be empty");
	}

	const SkyFieldOfView * fov = data->getFieldOfView();

	for (unsigned istar=0; istar<fov->getStars().size(); istar++) {

		Star * star = fov->getStars()[istar];
		SkyPosition starPosition = star->getSkyPosition();

		vector<double> xpositions,ypositions;
		int status = 0;
		char projection[5] = "-TAN";

		for (unsigned i=0; i<nJitterPerExposure; i++) {

			double rollAngle = 0.;
			SkyPosition pointingDirection = fov->getPointingDirection();
			if (rollAngles.size() > 0) { //OrbitSimulator has been run (with or without JitterProducer)
				rollAngle = rollAngles[i];
				pointingDirection = pointingDirections[i];
			} else if (jitterAPEs.size() > 0) {
				throw runtime_error("Error in configuration: Cannot run JitterProducer without also running OrbitSimulator");
			}

			double xpos = 0.;
			double ypos = 0.;
			//Replace RA values with 360-RA to put the movement on the CCD plane in the right direction and multiply the roll angle by -1
			//since the rotation is about the negative of the line of sight vector when facing the CCD (instructions from Andrew Hyslop 19th Jul)
			fits_world_to_pix(360.-starPosition.getRightAscension()/3600., starPosition.getDeclination()/3600.,
					360.-pointingDirection.getRightAscension()/3600., pointingDirection.getDeclination()/3600.,
					data->getSubarrayDimensions().m_targetLocationX, data->getSubarrayDimensions().m_targetLocationY, m_plateScale/3600., m_plateScale/3600., -rollAngle, projection, &xpos, &ypos, &status);
			//Truncate the CCD plane position values to 6 decimal places so that a star located at the centre stays exactly centred when there is roll but no jitter.
			//Otherwise (xpixel,ypixel) for the PSF position in PSFGenerator is unstable and can flip between (511,511), (511,512), (512, 511) and (512,512)
			xpositions.push_back(round(xpos*1.e6)/1.e6);
			ypositions.push_back(round(ypos*1.e6)/1.e6);

//			if (istar==0) {
//				cout << "FocalPlaneGenerator: star RA = " << 360.-starPosition.getRightAscension()/3600. << ", star Dec = " << starPosition.getDeclination()/3600.
//				     << ", pointing RA = " << 360.-pointingDirection.getRightAscension()/3600. << ", pointing Dec = " << pointingDirection.getDeclination()/3600.
//				     << ", roll angle = " << rollAngle << ", xpos = " << xpositions[i]-512 << ", ypos = " << ypositions[i]-512 << ", status=" << status;
//				if (jitterAPEs.size() > 0) cout << ", yaw = " << jitterAPEs[i].getYawOffset() << ", pitch = " << jitterAPEs[i].getPitchOffset();
//				cout << endl;
//			}

		}

		star->setFocalPlaneX(xpositions);
		star->setFocalPlaneY(ypositions);

	}

}
