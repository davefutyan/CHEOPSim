/*
 * Data.cxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

#include "REF_APP_GainCorrection.hxx"
#include "CreateFitsFile.hxx"

#include "Data.hxx"

const int Data::kProgramType = 90;
constexpr double Data::kRawHKCadence;
constexpr double Data::kAveragedHKCadence;

Data::Data(const TimeConfiguration& timeConf, string modulesToRun) : m_visit(Visit()) {

	m_modulesToRun = modulesToRun;
	m_timeConf = timeConf;
	m_fov = nullptr;
	m_wavelengthDependence = nullptr;
	m_idealLightCurveParams = new IdealLightCurveParams();
	m_breathingTemperatureRange = make_pair(0.,0.);
	m_MPSVisits = nullptr;
	m_mpsVisitConstraints = nullptr;
	m_stackedImageCount = 0;
	m_discardedImages = 0;
	m_stackedImageCube = nullptr;
	m_stackedMetaData = nullptr;
	m_stackedMetaData_unstacked = nullptr;
	m_unstackedImageCube = nullptr;
	m_unstackedMetaData = nullptr;
	m_imagetteCube = nullptr;
	m_imagetteMetaData = nullptr;
	m_stackedImageCube_DP = nullptr;
	m_stackedImageCube_smear = nullptr;
	m_stackedImageCube_cosmics = nullptr;
	m_centroid = nullptr;
	m_stackedTruthMetaData = nullptr;
	m_unstackedTruthMetaData = nullptr;
	m_imagetteTruthMetaData = nullptr;
	m_cosmicRayImage = nullptr;
	m_flatField = nullptr;
	m_flatFieldToSubtract = nullptr;
	m_badPixels = nullptr;
	m_badPixelsLeft = nullptr;
	m_badPixelsRight = nullptr;
	m_badPixelsTop = nullptr;
	m_photPixels = nullptr;
	m_photPixelsLeft = nullptr;
	m_photPixelsRight = nullptr;
	m_photPixelsTop = nullptr;
	m_telegraphicPixels = new vector<TelegraphicPixel>();
	m_frameTransferSmear = false;
	m_chargeTransferEol = false;
	m_doFullFrame = timeConf.doFullFrame();
	m_readoutTime = 0.;
	m_biasOffset = 0.;
	m_gainFilename = "";
	m_nominalVoltage[0] = -999.;
	m_omitSAA = false;
	m_omitEarthOccultation = false;
	m_omitStrayLight = false;
	m_numberOfCEsToOmitAtStart = 0;
	m_onlyWriteEveryNthCE = 1;
	m_redundantHardware = false;
	m_gainNonLinearity = false;
	m_serialReadRate = 0.;

    initializeStackedCcdMargins();
    initializeUnstackedCcdMargins();

}

bool Data::moduleIsActive(string moduleName) {

	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(",");
	tokenizer tokens(m_modulesToRun, sep);
	for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
		if (string(*tok_iter) == moduleName) return true;
	}
	return false;

}
void Data::initializeStackedCcdMargins() {

	m_stackedBlankLeftMarginCube = nullptr;
	m_stackedBlankRightMarginCube = nullptr;
	m_stackedDarkLeftMarginCube = nullptr;
	m_stackedDarkRightMarginCube = nullptr;
	m_stackedDarkTopMarginCube = nullptr;
	m_stackedOverscanLeftMarginCube = nullptr;
	m_stackedOverscanRightMarginCube = nullptr;
	m_stackedOverscanTopMarginCube = nullptr;

}

void Data::initializeUnstackedCcdMargins() {

	m_unstackedBlankLeftImageCube = nullptr;
	m_unstackedBlankRightImageCube = nullptr;
	m_unstackedDarkLeftImageCube = nullptr;
	m_unstackedDarkRightImageCube = nullptr;
	m_unstackedDarkTopImageCube = nullptr;
	m_unstackedOverscanLeftImageCube = nullptr;
	m_unstackedOverscanRightImageCube = nullptr;
	m_unstackedOverscanTopImageCube = nullptr;

}

Data::~Data() {

	delete m_fov;
	if (m_wavelengthDependence != nullptr) delete m_wavelengthDependence;
	if (m_idealLightCurveParams != nullptr) delete m_idealLightCurveParams;
	for (vector<SatelliteData*>::const_iterator it = m_satelliteData.begin(); it!=m_satelliteData.end(); ++it) delete (*it);

	if (m_MPSVisits != nullptr) delete m_MPSVisits;
	if (m_mpsVisitConstraints != nullptr) delete m_mpsVisitConstraints;
	if (m_unstackedImageCube != nullptr) delete m_unstackedImageCube;
	if (m_unstackedMetaData != nullptr) delete m_unstackedMetaData;
	if (m_stackedImageCube != nullptr) delete m_stackedImageCube;
	if (m_stackedMetaData != nullptr) delete m_stackedMetaData;
	if (m_stackedMetaData_unstacked != nullptr) delete m_stackedMetaData_unstacked;
	if (m_imagetteCube != nullptr) delete m_imagetteCube;
	if (m_imagetteMetaData != nullptr) delete m_imagetteMetaData;
	if (m_stackedImageCube_DP != nullptr) delete m_stackedImageCube_DP;
	if (m_stackedImageCube_smear != nullptr) delete m_stackedImageCube_smear;
	if (m_stackedImageCube_cosmics != nullptr) delete m_stackedImageCube_cosmics;
	if (m_centroid != nullptr) delete m_centroid;
	if (m_stackedTruthMetaData != nullptr) delete m_stackedTruthMetaData;
	if (m_unstackedTruthMetaData != nullptr) delete m_unstackedTruthMetaData;
	if (m_imagetteTruthMetaData != nullptr) delete m_imagetteTruthMetaData;
	if (m_cosmicRayImage != nullptr) delete m_cosmicRayImage;
	if (m_flatField != nullptr) delete m_flatField;
	if (m_flatFieldToSubtract != nullptr) delete m_flatFieldToSubtract;
	if (m_badPixelsLeft != nullptr) delete m_badPixelsLeft;
	if (m_badPixelsRight != nullptr) delete m_badPixelsRight;
	if (m_badPixelsTop != nullptr) delete m_badPixelsTop;
	if (m_photPixels != nullptr) delete m_photPixels;
	if (m_photPixelsLeft != nullptr) delete m_photPixelsLeft;
	if (m_photPixelsRight != nullptr) delete m_photPixelsRight;
	if (m_photPixelsTop != nullptr) delete m_photPixelsTop;
	if (m_badPixels != nullptr) delete m_badPixels;
	delete m_telegraphicPixels;

	deleteStackedCcdMargins();
	deleteUnstackedCcdMargins();

}

void Data::clearImages() {
	for (vector<Image*>::const_iterator it = m_images.begin(); it!=m_images.end(); ++it) delete (*it);
	for (vector<Image*>::const_iterator it = m_stackedImages.begin(); it!=m_stackedImages.end(); ++it) delete (*it);
	for (vector<Image*>::const_iterator it = m_imagettes.begin(); it!=m_imagettes.end(); ++it) delete (*it);
	for (vector<Image*>::const_iterator it = m_stackedImagettes.begin(); it!=m_stackedImagettes.end(); ++it) delete (*it);
	for (vector<Image*>::const_iterator it = m_imagesToSmearUp.begin(); it!=m_imagesToSmearUp.end(); ++it) delete (*it);
	for (vector<Image*>::const_iterator it = m_imagesToSmearDown.begin(); it!=m_imagesToSmearDown.end(); ++it) delete (*it);
	m_images.clear();
	m_stackedImages.clear();
	m_imagettes.clear();
	m_stackedImagettes.clear();
	m_imagesToSmearUp.clear();
	m_imagesToSmearDown.clear();
	if (m_cosmicRayImage != nullptr) {
		delete m_cosmicRayImage;
		m_cosmicRayImage = nullptr;
	}
}

string Data::getOutputDirectory() const {

	string progType = boost::lexical_cast<string>(m_visit.m_progType);

	string progId = boost::lexical_cast<string>(m_visit.m_progId);
	for (int i=3;i>0;i--) {
		if (m_visit.m_progId<pow(10,i)) {
			progId = "0"+progId;
		} else {
			break;
		}
	}
	string reqId = boost::lexical_cast<string>(m_visit.m_reqId);
	for (int i=3;i>0;i--) {
		if (m_visit.m_reqId<pow(10,i)) {
			reqId = "0"+reqId;
		} else {
			break;
		}
	}
	string visitCtr = boost::lexical_cast<string>(m_visit.m_visitCtr);
	for (int i=1;i>0;i--) {
		if (m_visit.m_visitCtr<pow(10,i)) {
			visitCtr = "0"+visitCtr;
		} else {
			break;
		}
	}
	string versionNum = boost::lexical_cast<string>(m_visit.m_versionNum);
	for (int i=2;i>0;i--) {
		if (m_visit.m_versionNum<pow(10,i)) {
			versionNum = "0"+versionNum;
		} else {
			break;
		}
	}

	return "CH_PR"+progType+progId+"_TG"+reqId+visitCtr;

}

SkyFieldOfView* Data::getFieldOfView() const {
	if (m_fov) {
		return m_fov;
	} else {
		throw runtime_error("Error in Data::getFieldOfView() :Field of view undefined");
	}
}

void Data::setSubarrayDimensions(int xDim, int yDim, int xOffset, int yOffset, double targetLocationX, double targetLocationY) {
	m_subarrayDimensions.m_xDim = xDim;
	m_subarrayDimensions.m_yDim = yDim;
	m_subarrayDimensions.m_xOffset = xOffset;
	m_subarrayDimensions.m_yOffset = yOffset;
	setTargetLocation(targetLocationX,targetLocationY);
}

void Data::setTargetLocation(double targetLocationX, double targetLocationY) {
	m_subarrayDimensions.m_targetLocationX = targetLocationX;
	m_subarrayDimensions.m_targetLocationY = targetLocationY;
}

void Data::setPhotometryParams(double radius_barycentre, double radius_bkgInner, double radius_bkgOuter, double radius_psf, bool subtractBackground) {
	m_photometryParams.m_radius_barycentre = radius_barycentre;
	m_photometryParams.m_radius_bkgInner = radius_bkgInner;
	m_photometryParams.m_radius_bkgOuter = radius_bkgOuter;
	m_photometryParams.m_radius_psf = radius_psf;
	m_photometryParams.m_subtractBackground = subtractBackground;
}

void Data::setFitsImageCube(SciRawImagette* imageCube, SciRawImagettemetadata* metaData) {
	m_imagetteCube = imageCube;
	m_imagetteMetaData = metaData;
}

void Data::setFitsImageCube(SciRawSubarray* imageCube, SciRawImagemetadata* metaData, SciRawUnstackedimagemetadata* metaData_unstacked) {
	m_stackedImageCube = imageCube;
	m_stackedMetaData = metaData;
	m_stackedMetaData_unstacked = metaData_unstacked;
}

void Data::setFitsImageCube(SimRawUnstackedsubarray* imageCube, SciRawImagemetadata* metaData) {
	m_unstackedImageCube = imageCube;
	m_unstackedMetaData = metaData;
}

void Data::setFitsImageCube(SimRawDoubleprecisionsubarray* imageCube, SciRawImagemetadata* metaData, SciRawUnstackedimagemetadata* metaData_unstacked) {
	m_stackedImageCube_DP = imageCube;
	m_stackedMetaData = metaData;
	m_stackedMetaData_unstacked = metaData_unstacked;
}

void Data::setFitsImageCube_smear(SimRawDoubleprecisionsubarray* imageCube) {
	m_stackedImageCube_smear = imageCube;
}

void Data::setFitsImageCube_cosmics(SimRawDoubleprecisionsubarray* imageCube) {
	m_stackedImageCube_cosmics = imageCube;
}

void Data::setTruthMetaData(SimTruSubarray* truthMetaData, IMAGETYPE type) {
	switch (type) {
	case stacked :
		m_stackedTruthMetaData = truthMetaData;
		break;
	case unstacked :
		m_unstackedTruthMetaData = truthMetaData;
		break;
	case imagette :
		m_imagetteTruthMetaData = truthMetaData;
		break;
	default :
		cout << "ERROR in Data::setTruthMetaData: invalid image type"<< endl;
	}
}

SatelliteData* Data::getSatelliteData(int timeStep, bool create, bool checkTimeStep) {

	unsigned absoluteTimeStep = timeStep + m_timeConf.getNumberOfPreSubArrayTimeSteps();
	//cout << m_satelliteData.size() << " " << absoluteTimeStep << " " << m_timeConf.getNumberOfTimeSteps() << " " << m_timeConf.getNumberOfPreSubArrayTimeSteps() << endl;

	if (create) {
		if (m_satelliteData.size()==absoluteTimeStep) {
			SatelliteData * satData = new SatelliteData();
			m_satelliteData.push_back(satData);
			return satData;
		} else {
			cout << m_satelliteData.size() << " " << timeStep << " " << absoluteTimeStep << " " << m_timeConf.getNumberOfTimeSteps()+m_timeConf.getNumberOfPreSubArrayTimeSteps() << endl;
			throw runtime_error("Error in Data::getSatelliteData (create mode) : Satellite data array size does not match the current time step - abort");
		}
	} else {
		if (m_satelliteData.size()==absoluteTimeStep+1 ||
			m_satelliteData.size()==m_timeConf.getNumberOfTimeSteps()+m_timeConf.getNumberOfPreSubArrayTimeSteps() || (!checkTimeStep && m_satelliteData.size()>absoluteTimeStep)) {
			return m_satelliteData[absoluteTimeStep];
		} else {
			cout << m_satelliteData.size() << " " << timeStep << " " << absoluteTimeStep << " " << m_timeConf.getNumberOfTimeSteps()+m_timeConf.getNumberOfPreSubArrayTimeSteps() << endl;
			throw runtime_error("Error in Data::getSatelliteData : Satellite data array size not consistent with the number of time steps - abort");
		}
	}

}

SciRawImagemetadata* Data::getFitsImageMetaData(IMAGETYPE type) {
	switch (type) {
	case stacked : return m_stackedMetaData; break;
	case unstacked : return m_unstackedMetaData; break;
	default :
		cout << "ERROR in Data::getFitsImageMetaData: invalid image type" << endl;
		return nullptr;
	}
}

SciRawUnstackedimagemetadata* Data::getFitsUnstackedImageMetaData(IMAGETYPE type) {
	switch (type) {
	case stacked : return m_stackedMetaData_unstacked; break;
	default :
		cout << "ERROR in Data::getFitsUnstackedImageMetaData: invalid image type" << endl;
		return nullptr;
	}
}

SimTruSubarray* Data::getFitsTruthMetaData(IMAGETYPE type) {
	switch (type) {
	case stacked : return m_stackedTruthMetaData; break;
	case unstacked : return m_unstackedTruthMetaData; break;
	case imagette : return m_imagetteTruthMetaData; break;
	default :
		cout << "ERROR in Data::getFitsTruthMetaData: invalid image type" << endl;
		return nullptr;
	}
}

void Data::resetFitsImageCube(IMAGETYPE type) {
	if (type==unstacked) {
	    if (m_unstackedImageCube != nullptr) delete m_unstackedImageCube;
	    if (m_unstackedMetaData != nullptr) delete m_unstackedMetaData;
	    if (m_unstackedTruthMetaData != nullptr) delete m_unstackedTruthMetaData;
	    m_unstackedImageCube = nullptr;
	    m_unstackedMetaData = nullptr;
	    m_unstackedTruthMetaData = nullptr;
	    deleteUnstackedCcdMargins();
	    initializeUnstackedCcdMargins();
	} else if (type==stacked) {
		if (m_stackedImageCube != nullptr) delete m_stackedImageCube;
		if (m_stackedImageCube_DP != nullptr) delete m_stackedImageCube_DP;
		if (m_stackedMetaData != nullptr) delete m_stackedMetaData;
		if (m_stackedMetaData_unstacked != nullptr) delete m_stackedMetaData_unstacked;
		if (m_stackedTruthMetaData != nullptr) delete m_stackedTruthMetaData;
		if (m_stackedImageCube_smear != nullptr) delete m_stackedImageCube_smear;
		if (m_stackedImageCube_cosmics != nullptr) delete m_stackedImageCube_cosmics;
	    m_stackedImageCube = nullptr;
	    m_stackedImageCube_DP = nullptr;
	    m_stackedMetaData = nullptr;
	    m_stackedMetaData_unstacked = nullptr;
	    m_stackedTruthMetaData = nullptr;
	    m_stackedImageCube_smear = nullptr;
	    m_stackedImageCube_cosmics = nullptr;
	    deleteStackedCcdMargins();
	    initializeStackedCcdMargins();
	} else if (type==imagette) {
		if (m_imagetteCube != nullptr) delete m_imagetteCube;
		if (m_imagetteMetaData != nullptr) delete m_imagetteMetaData;
		if (m_imagetteTruthMetaData != nullptr) delete m_imagetteTruthMetaData;
	    m_imagetteCube = nullptr;
	    m_imagetteMetaData = nullptr;
	    m_imagetteTruthMetaData = nullptr;
	}
}

unsigned Data::getNumberOfImagesPerCube(IMAGETYPE type) const {
	unsigned imagesPerCube;
	if (type==stacked || type==imagette) {
		unsigned numberOfStackedImages = m_timeConf.getNumberOfStackedImagesPITL();
//		imagesPerCube = numberOfStackedImages/m_timeConf.getNumberOfStackedImageCubes();
//		int remainder = numberOfStackedImages%m_timeConf.getNumberOfStackedImageCubes();
//		if (remainder != 0) imagesPerCube += 1;
		imagesPerCube = numberOfStackedImages; // For now, always use a single image cube for stacked images and imagettes, since DRP cannot currently work with multiple input files
		if (type==imagette) {
			imagesPerCube *= m_timeConf.getExposuresPerStack()/m_timeConf.getImagetteStackingNumber();
			//If there is not a whole number of imagette stacks per sub-array stack, there is an additional partial imagette stack for each sub-array stack.
			if (m_timeConf.getExposuresPerStack()%m_timeConf.getImagetteStackingNumber() != 0) {
				imagesPerCube += numberOfStackedImages;
				//This shouldn't happen in the current implementation of TimeConfiguration::setImagetteStackingNumber() which sets the imagette stacking number
				//according to CHEOPSim_web/LUT_image_stacking.txt, which should guarantee that exposuresPerStack/imagetteStackingNumber is a whole number,
				//so throw a runtime error while this is the implementation
				throw runtime_error("ERROR in Data::getNumberOfImagesPerCube(imagette): exposuresPerStack/imagetteStackingNumber is not a whole number");
			}
		}
	} else {
		unsigned numberOfTimeSteps = m_timeConf.getNumberOfTimeStepsPITL();
		imagesPerCube = numberOfTimeSteps/m_timeConf.getNumberOfUnstackedImageCubes();
		int remainder = numberOfTimeSteps%m_timeConf.getNumberOfUnstackedImageCubes();
		if (remainder != 0) imagesPerCube += 1;
		//Force images per cube to be a multiple of the number of images per stack
		//so that writing out of unstacked images can be done once per stacked image
		if (imagesPerCube%m_timeConf.getExposuresPerStack() != 0) imagesPerCube = ((imagesPerCube/m_timeConf.getExposuresPerStack())+1)*m_timeConf.getExposuresPerStack();
	}
	return imagesPerCube;
}

void Data::deleteStackedCcdMargins() {

	if (m_stackedBlankLeftMarginCube != nullptr) delete m_stackedBlankLeftMarginCube;
	if (m_stackedBlankRightMarginCube != nullptr) delete m_stackedBlankRightMarginCube;
	if (m_stackedDarkLeftMarginCube != nullptr) delete m_stackedDarkLeftMarginCube;
	if (m_stackedDarkRightMarginCube != nullptr) delete m_stackedDarkRightMarginCube;
	if (m_stackedDarkTopMarginCube != nullptr) delete m_stackedDarkTopMarginCube;
	if (m_stackedOverscanLeftMarginCube != nullptr) delete m_stackedOverscanLeftMarginCube;
	if (m_stackedOverscanRightMarginCube != nullptr) delete m_stackedOverscanRightMarginCube;
	if (m_stackedOverscanTopMarginCube != nullptr) delete m_stackedOverscanTopMarginCube;

}

void Data::deleteUnstackedCcdMargins() {

	if (m_unstackedBlankLeftImageCube != nullptr) delete m_unstackedBlankLeftImageCube;
	if (m_unstackedBlankRightImageCube != nullptr) delete m_unstackedBlankRightImageCube;
	if (m_unstackedDarkLeftImageCube != nullptr) delete m_unstackedDarkLeftImageCube;
	if (m_unstackedDarkRightImageCube != nullptr) delete m_unstackedDarkRightImageCube;
	if (m_unstackedDarkTopImageCube != nullptr) delete m_unstackedDarkTopImageCube;
	if (m_unstackedOverscanLeftImageCube != nullptr) delete m_unstackedOverscanLeftImageCube;
	if (m_unstackedOverscanRightImageCube != nullptr) delete m_unstackedOverscanRightImageCube;
	if (m_unstackedOverscanTopImageCube != nullptr) delete m_unstackedOverscanTopImageCube;

}

void Data::setFitsCcdMargins(SciRawBlankleft * blankLeftMarginCube,
					   	     SciRawBlankright * blankRightMarginCube,
							 SciRawDarkleft * darkLeftMarginCube,
							 SciRawDarkright * darkRightMarginCube,
							 SciRawDarktop * darkTopMarginCube,
							 SciRawOverscanleft * overscanLeftMarginCube,
							 SciRawOverscanright * overscanRightMarginCube,
							 SciRawOverscantop * overscanTopMarginCube) {

	m_stackedBlankLeftMarginCube = blankLeftMarginCube;
	m_stackedBlankRightMarginCube = blankRightMarginCube;
	m_stackedDarkLeftMarginCube = darkLeftMarginCube;
	m_stackedDarkRightMarginCube = darkRightMarginCube;
	m_stackedDarkTopMarginCube = darkTopMarginCube;
	m_stackedOverscanLeftMarginCube = overscanLeftMarginCube;
	m_stackedOverscanRightMarginCube = overscanRightMarginCube;
	m_stackedOverscanTopMarginCube = overscanTopMarginCube;

}

void Data::setFitsCcdMarginImages(SimRawUnstackedblankleftimage* blankLeftImageCube,
								  SimRawUnstackedblankrightimage* blankRightImageCube,
								  SimRawUnstackeddarkleftimage* darkLeftImageCube,
								  SimRawUnstackeddarkrightimage* darkRightImageCube,
								  SimRawUnstackeddarktopimage* darkTopImageCube,
								  SimRawUnstackedoverscanleftimage* overscanLeftImageCube,
								  SimRawUnstackedoverscanrightimage* overscanRightImageCube,
								  SimRawUnstackedoverscantopimage* overscanTopImageCube) {

	m_unstackedBlankLeftImageCube = blankLeftImageCube;
	m_unstackedBlankRightImageCube = blankRightImageCube;
	m_unstackedDarkLeftImageCube = darkLeftImageCube;
	m_unstackedDarkRightImageCube = darkRightImageCube;
	m_unstackedDarkTopImageCube = darkTopImageCube;
	m_unstackedOverscanLeftImageCube = overscanLeftImageCube;
	m_unstackedOverscanRightImageCube = overscanRightImageCube;
	m_unstackedOverscanTopImageCube = overscanTopImageCube;

}

void Data::setNominalVoltages(string gainFilename) {

	m_gainFilename = gainFilename;

	//In CHEOPS-UGE-SYS-PR-019 Issue 3 Revision 1, nominal voltage values are defined relative VSS, so add nominal VSS to get the absolute value
	RefAppGaincorrection * gainCorrection_file = new RefAppGaincorrection(string(getenv("CHEOPS_SW"))+"/resources/"+gainFilename,"READONLY");
	m_nominalVoltage[SatelliteData::VSS] = gainCorrection_file->getKeyVssOff();
	m_nominalVoltage[SatelliteData::VOD] = gainCorrection_file->getKeyVodOff() + m_nominalVoltage[SatelliteData::VSS];
	m_nominalVoltage[SatelliteData::VRD] = gainCorrection_file->getKeyVrdOff() + m_nominalVoltage[SatelliteData::VSS];
	m_nominalVoltage[SatelliteData::VOG] = gainCorrection_file->getKeyVogOff() + m_nominalVoltage[SatelliteData::VSS];
	m_nominalVoltage[SatelliteData::TEMP] = gainCorrection_file->getKeyTempOff();

}

double Data::getTemperatureDependentVoltage(double temperature, SatelliteData::VOLTAGE_TYPE type) const {

	if (m_nominalVoltage[type] == -999.) {
		throw runtime_error("ERROR in Data::getTemperatureDependentVoltage: nominal values uninitialized. Data::setNominalVoltages must be called first.");
	}

	switch(type) {
	case SatelliteData::VOD:
		return m_nominalVoltage[SatelliteData::VOD] + (temperature-SatelliteData::kDefaultFeeBiasTemperature)*(kVODTempCoeff/1.E6);
	case SatelliteData::VRD:
		return m_nominalVoltage[SatelliteData::VRD] + (temperature-SatelliteData::kDefaultFeeBiasTemperature)*(kVRDTempCoeff/1.E6);
	case SatelliteData::VOG:
		return m_nominalVoltage[SatelliteData::VOG] + (temperature-SatelliteData::kDefaultFeeBiasTemperature)*(kVOGTempCoeff/1.E6);
	case SatelliteData::VSS:
		return m_nominalVoltage[SatelliteData::VSS] + (temperature-SatelliteData::kDefaultFeeBiasTemperature)*(kVSSTempCoeff/1.E6);
	default:
		throw runtime_error("ERROR in Data::getTemperatureDependentVoltage: Invalid voltage type");
	}
}

double Data::getVoltageWithDrift(double timeSinceVisitStart, SatelliteData::VOLTAGE_TYPE type) const {

	if (type != SatelliteData::VOD && type != SatelliteData::VRD &&
		type != SatelliteData::VOG && type != SatelliteData::VSS) {
		throw runtime_error("ERROR in Data::getVoltageWithDrift: invalid voltage type: "+type);
	}

	if (m_nominalVoltage[type] == -999.) {
		throw runtime_error("ERROR in Data::getVoltageWithDrift: nominal values uninitialized. Data::setNominalVoltages must be called first.");
	}

	//The voltage drifts are defined in the payload calibration document for voltages relative to VSS, so need to
	//take into account VSS nominal and drifted values to calculate absolute drifted values for VOD, VOG and VRD
	double vssNominal = 0.;
	double vssWithDrift = 0.;
	if (type!=SatelliteData::VSS) {
		vssNominal = m_nominalVoltage[SatelliteData::VSS];
		vssWithDrift = vssNominal + (timeSinceVisitStart/(24.*3600.))*(SatelliteData::kVoltageDrift[m_redundantHardware][SatelliteData::VSS]/1.E6);
	}

	double voltageWithDrift_relative = (m_nominalVoltage[type] - vssNominal) + (timeSinceVisitStart/(24.*3600.))*(SatelliteData::kVoltageDrift[m_redundantHardware][type]/1.E6);
	double voltageWithDrift = voltageWithDrift_relative + vssWithDrift;
	//cout << m_nominalVoltage[type] << " " <<  m_timeConf.getTimeSinceVisitStart(timeStep)
	//		<< " " << SatelliteData::kVoltageDrift[m_redundantHardware][type] << " " << voltageWithDrift << endl;
	return voltageWithDrift;

}

void Data::initializeBadPixelMap(string reference_file) {

	if (m_badPixels == nullptr) {

		string dirname = getOutputDirectory()+"/data";
		if (!boost::filesystem::exists(dirname)) boost::filesystem::create_directory(dirname);
		boost::posix_time::ptime startTime = getTimeConfiguration().getStartTime();
		double duration = getTimeConfiguration().getDuration();
		long seconds = static_cast<long>(lround((floor(duration))));
		long milliseconds = static_cast<long>(lround((duration - floor(duration))*1000.));
		boost::posix_time::ptime endTime = startTime + boost::posix_time::seconds(seconds) + boost::posix_time::milliseconds(milliseconds);

		//startTime = boost::posix_time::ptime(boost::gregorian::date(2018,1,1));
		//endTime = boost::posix_time::ptime(boost::gregorian::date(2050,1,1));

		m_badPixels = new RefAppBadpixelmap(buildFileName(dirname, to_iso_extended_string(startTime),
				VisitId(),PassId(), std::string(), RefAppBadpixelmap::getExtName()), "CREATE");
		m_badPixelsLeft = Append_RefAppBadpixelmapleft(m_badPixels, {Image::kNDarkCols,Image::kYDim});
		m_badPixelsRight = Append_RefAppBadpixelmapright(m_badPixels, {Image::kNDarkCols,Image::kYDim});
		m_badPixelsTop = Append_RefAppBadpixelmaptop(m_badPixels, {Image::kXDim,Image::kNDarkRows});
		m_photPixels = Append_RefAppPhotpixelmap(m_badPixels, {Image::kXDim,Image::kYDim});
		m_photPixelsLeft = Append_RefAppPhotpixelmapleft(m_badPixels, {Image::kNDarkCols,Image::kYDim});
		m_photPixelsRight = Append_RefAppPhotpixelmapright(m_badPixels, {Image::kNDarkCols,Image::kYDim});
		m_photPixelsTop = Append_RefAppPhotpixelmaptop(m_badPixels, {Image::kXDim,Image::kNDarkRows});

		m_badPixels->setKeyArchRev(0);
		m_badPixels->setKeyProcNum(m_visit.m_versionNum);
		m_badPixels->setKeyVStrtU(to_iso_extended_string(startTime));
		m_badPixels->setKeyVStopU(to_iso_extended_string(endTime));
		m_badPixels->setKeyProvider("CHEOPSim");
		m_badPixels->setKeyDescrip("Bad pixels assigned according to CHEOPSim configuration");

		m_badPixelsLeft->setKeyArchRev(0);
		m_badPixelsLeft->setKeyProcNum(m_visit.m_versionNum);
		m_badPixelsLeft->setKeyVStrtU(to_iso_extended_string(startTime));
		m_badPixelsLeft->setKeyVStopU(to_iso_extended_string(endTime));

		m_badPixelsRight->setKeyArchRev(0);
		m_badPixelsRight->setKeyProcNum(m_visit.m_versionNum);
		m_badPixelsRight->setKeyVStrtU(to_iso_extended_string(startTime));
		m_badPixelsRight->setKeyVStopU(to_iso_extended_string(endTime));

		m_badPixelsTop->setKeyArchRev(0);
		m_badPixelsTop->setKeyProcNum(m_visit.m_versionNum);
		m_badPixelsTop->setKeyVStrtU(to_iso_extended_string(startTime));
		m_badPixelsTop->setKeyVStopU(to_iso_extended_string(endTime));

		m_photPixels->setKeyArchRev(0);
		m_photPixels->setKeyProcNum(m_visit.m_versionNum);
		m_photPixels->setKeyVStrtU(to_iso_extended_string(startTime));
		m_photPixels->setKeyVStopU(to_iso_extended_string(endTime));

		m_photPixelsLeft->setKeyArchRev(0);
		m_photPixelsLeft->setKeyProcNum(m_visit.m_versionNum);
		m_photPixelsLeft->setKeyVStrtU(to_iso_extended_string(startTime));
		m_photPixelsLeft->setKeyVStopU(to_iso_extended_string(endTime));

		m_photPixelsRight->setKeyArchRev(0);
		m_photPixelsRight->setKeyProcNum(m_visit.m_versionNum);
		m_photPixelsRight->setKeyVStrtU(to_iso_extended_string(startTime));
		m_photPixelsRight->setKeyVStopU(to_iso_extended_string(endTime));

		m_photPixelsTop->setKeyArchRev(0);
		m_photPixelsTop->setKeyProcNum(m_visit.m_versionNum);
		m_photPixelsTop->setKeyVStrtU(to_iso_extended_string(startTime));
		m_photPixelsTop->setKeyVStopU(to_iso_extended_string(endTime));

	}

	if (reference_file.find("Dark") != std::string::npos) {
		m_badPixels->setKeyDarkRf(reference_file);
	} else if (reference_file.find("Flat") != std::string::npos) {
		m_badPixels->setKeyFfRf(reference_file);
	} else if (reference_file.find("GainCorrection") != std::string::npos) {
		m_badPixels->setKeyGainRf(reference_file);
	}

}

void Data::setBadPixel(int i, int j, int value) {

	if (m_badPixels == nullptr) throw runtime_error("Error in Data::setBadPixel: bad pixel map data structure was not initialized");

	if (i>=Image::kNBlankCols+Image::kNOverscanCols && i<Image::kLeftMargin && j<Image::kYDim) {
		(*m_badPixelsLeft)[j][i-Image::kNBlankCols-Image::kNOverscanCols] = value;
		(*m_photPixelsLeft)[j][i-Image::kNBlankCols-Image::kNOverscanCols] = (value==-1 || value==0 || value==1) ? 0 : 1;
	} else if (i>=Image::kXDim+Image::kLeftMargin && i<Image::kXDim+Image::kLeftMargin+Image::kNDarkCols && j<Image::kYDim) {
		(*m_badPixelsRight)[j][i-Image::kXDim-Image::kLeftMargin] = value;
		(*m_photPixelsRight)[j][i-Image::kXDim-Image::kLeftMargin] = (value==-1 || value==0 || value==1) ? 0 : 1;
	} else if (i>=Image::kLeftMargin && i<Image::kXDim+Image::kLeftMargin && j>=Image::kYDim && j<Image::kYDim+Image::kNDarkRows) {
		(*m_badPixelsTop)[j-Image::kYDim][i-Image::kLeftMargin] = value;
		(*m_photPixelsTop)[j-Image::kYDim][i-Image::kLeftMargin] = (value==-1 || value==0 || value==1) ? 0 : 1;
	} else  if (i>=Image::kLeftMargin && i<Image::kXDim+Image::kLeftMargin && j<Image::kYDim) {
		(*m_badPixels)[j][i-Image::kLeftMargin] = value;
		(*m_photPixels)[j][i-Image::kLeftMargin] = (value==-1 || value==0 || value==1) ? 0 : 1;
	}

}

void Data::calculateHKAverages() {

	//Loop over the total duration in steps of 20s (the cadence for HK data sent to ground)
	double startTimeOffset = (m_timeConf.getStartTimeUTC() - m_timeConf.getVisitStartTimeUTC()).getSeconds();
	for (double time = 0.; time <= startTimeOffset+m_timeConf.getDuration()+kAveragedHKCadence; time+=kAveragedHKCadence) {

		//Identify the closest entry in the raw HK vector (1.2s cadence)
		int iHKRaw = static_cast<int>(lround((time)/kRawHKCadence));
		if (iHKRaw>int(m_rawHKData.size())-1) {
			throw runtime_error("Error in Data::calculateHKAverage: no entry in raw HK data vector corresponding to requested time");
		}
		//cout << endl << "iKHRaw " << time << " " << iHKRaw << endl;

		//For each voltage and temperature, calculate the mean over the previous 16 raw measurements
		//(HK to ground cadence / onboard measurement cadence = 20/1.2 = 16.67)
		unsigned nHK = 0;
		double ccdTemp = 0.;
		double biasTemp = 0.;
		double adcTemp = 0.;
		double vod = 0.;
		double vrd = 0.;
		double vog = 0.;
		double vss = 0.;
		double voltages[4];
		for (int i=iHKRaw; i>=0 && i>iHKRaw-16; i--) {
			//cout << i << endl;
			nHK++;
			HKData hkData = m_rawHKData[i];
			ccdTemp += hkData.m_ccdTemp;
			biasTemp += hkData.m_biasTemp;
			adcTemp += hkData.m_adcTemp;
			vod += hkData.m_vod;
			vrd += hkData.m_vrd;
			vog += hkData.m_vog;
			vss += hkData.m_vss;
		}
		ccdTemp /= nHK;
		biasTemp /= nHK;
		adcTemp /= nHK;
		voltages[SatelliteData::VOD] = vod/nHK;
		voltages[SatelliteData::VRD] = vrd/nHK;
		voltages[SatelliteData::VOG] = vog/nHK;
		voltages[SatelliteData::VSS] = vss/nHK;
		//cout << nHK << " " << setprecision(10) << ccdTemp << endl;

		m_averagedHKData.push_back(HKData(ccdTemp,biasTemp,adcTemp,voltages));

	}

}

Data::HKData Data::getClosestRawHKData(double timeSinceVisitStart) const {

	if (timeSinceVisitStart < 0) throw runtime_error("Error in Data::getClosestRawHKData: timeSinceVisitStart argument cannot be negative");

	//Identify the closest entry in the raw HK vector (1.2s cadence)
	unsigned index = lround((timeSinceVisitStart)/kRawHKCadence);
	//cout << index << " " << m_rawHKData.size() << endl;
	if (index > m_rawHKData.size()-1) {
		throw runtime_error("Error in Data::getClosestRawHKData: no entry in raw HK data vector corresponding to requested time");
	}

	return m_rawHKData[index];

}

Data::HKData Data::getClosestAveragedHKData(double timeSinceVisitStart) const {

	if (timeSinceVisitStart < 0) throw runtime_error("Error in Data::getClosestAveragedHKData: timeSinceVisitStart argument cannot be negative");

	//Identify the closest entry in the averaged HK vector (20s cadence)
	unsigned index = lround((timeSinceVisitStart)/kAveragedHKCadence);
	//cout << "getClosestAveragedHKData: " << index << " " << m_averagedHKData.size() << " " << m_averagedHKData[index].m_ccdTemp << endl;
	if (index > m_averagedHKData.size()-1) {
		throw runtime_error("Error in Data::getClosestAveragedHKData: no entry in averaged HK data vector corresponding to requested time");
	}

	return m_averagedHKData[index];

}
