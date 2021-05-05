/*
 * ImageWriter.cxx
 *
 *  Created on: 17 Nov 2014
 *      Author: futyand
 */

#include <sstream>

#include "boost/date_time/gregorian/gregorian_types.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"
#include "boost/accumulators/accumulators.hpp"
#include <boost/accumulators/statistics/stats.hpp>
#include "boost/accumulators/statistics/mean.hpp"
#include "boost/accumulators/statistics/variance.hpp"
using namespace boost::accumulators;

#include "SIM_TRU_FullArray.hxx"
#include "SIM_RAW_DoublePrecisionSubArray.hxx"
#include "SCI_RAW_Centroid.hxx"
#include "REF_APP_GainCorrection.hxx"
#include "CreateFitsFile.hxx"

#include "CommonTools.hxx"

#include "ImageWriter.hxx"

void ImageWriter::initialize(const ModuleParams & params) {

	m_writeUnstackedImages = params.GetAsBool("writeUnstackedImages");
	m_writeStackedImages = params.GetAsBool("writeStackedImages");
	m_writeImagettes = params.GetAsBool("writeImagettes");
	m_dynamicImagettes = params.GetAsBool("dynamicImagettes");
	m_stackImagettes = params.GetAsBool("stackImagettes");
	m_imagetteShape = params.GetAsString("imagetteShape");
	m_subarrayShape = params.GetAsString("subarrayShape");
	m_doublePrecisionStackedImages = params.GetAsBool("doublePrecisionStackedImages");
	m_writeCentroid = params.GetAsBool("writeCentroid");
	m_writeTruthData = params.GetAsBool("writeTruthData");
	m_stackingMethod = params.GetAsString("stackingMethod");
	m_marginStackingMethod = params.GetAsString("marginStackingMethod");
	m_imagetteSize = params.GetAsInt("imagetteSize");
	m_truthBarycentre = params.GetAsBool("truthBarycentre");
	m_targetLocationFromJitter = params.GetAsBool("targetLocationFromJitter");
	m_ceCounterOffset = params.GetAsInt("ceCounterOffset");

}

void ImageWriter::doBegin(Data* data, bool fullFrame) {

	//Initialize the visit data
	Data::Visit visit = data->getVisit();
	m_visit = visit;
	unsigned exposuresPerStack = data->getTimeConfiguration().getExposuresPerStack();

	cout << "  Margin mode: " << visit.m_marginMode << endl;
	if (visit.m_marginMode == "image") {
		m_marginMode = image;
	} else if (visit.m_marginMode == "reduced") {
		m_marginMode = reduced;
	} else {
		m_marginMode = total_collapsed;
	}

	//Define the spectral type to be used for header keywords
	m_specType = "G5"; //Default to spectral type to G5 if there are no stars in the FOV
	m_TEff = 5660.; //Default to effective temperature 5660K if there are no stars in the FOV
	m_Gmag = 9.; //Default to 9th magnitude if there are no stars in the FOV
	m_cheopsMag = 9.; //Default to 9th magnitude if there are no stars in the FOV
	m_GmagErr = 0.;
	m_cheopsMagErr = 0.;
	if (data->getFieldOfView()->getStars().size()>0) {
		Star * targetStar = data->getFieldOfView()->getStars()[0];
		m_specType = targetStar->getSpectralTypeString();
		m_TEff = targetStar->getEffectiveTemperature();
		m_Gmag = targetStar->getGaiaMagnitude();
		m_cheopsMag = targetStar->getCheopsMagnitude();
		m_GmagErr = targetStar->getGaiaMagnitudeError();
		m_cheopsMagErr = targetStar->getCheopsMagnitudeError();
	}

	//Get the norminal gain to be used for SCI_RAW_UnstackedImageMetadata
	RefAppGaincorrection * gainCorrection_file = new RefAppGaincorrection(string(getenv("CHEOPS_SW"))+"/resources/"+data->getGainFilename(),"READONLY");
	m_nominalGain = gainCorrection_file->getKeyGainNom();
	delete gainCorrection_file;

	//Generate the output directory structure
	m_dirname = data->getOutputDirectory();
	if (m_writeStackedImages || (m_writeImagettes && exposuresPerStack>1)) boost::filesystem::create_directory(m_dirname+"/data"); //Only write imagettes if subarrays are stacked
	if (m_writeTruthData && !boost::filesystem::exists(m_dirname+"/truth")) boost::filesystem::create_directory(m_dirname+"/truth");
	if (m_writeUnstackedImages) boost::filesystem::create_directory(m_dirname+"/dfs");

	//Initialize the Photometry class for barycentre determination
	m_photometry = new Photometry(data->getSubarrayDimensions(),data->getPhotometryParams());

	//Determine the number of imagettes corresponding to a stacked sub-array.
	//If there is not a whole number of imagette stacks per sub-array stack, there is an additional partial imagette stack.
	if (m_stackImagettes) {
		data->setImagetteStackingNumber();
		if (data->getTimeConfiguration().getImagetteStackingNumber()==1) m_stackImagettes=false;
	}
	data->getMpsPreVisits()->setCellNexpImagettes(exposuresPerStack>1 ? data->getTimeConfiguration().getImagetteStackingNumber() : 0);
	cout << "  Imagette stacking number: " << (exposuresPerStack>1 ? data->getTimeConfiguration().getImagetteStackingNumber() : 0) << endl;
	m_imagettesPerStackedSubarray = exposuresPerStack/data->getTimeConfiguration().getImagetteStackingNumber();
	if (exposuresPerStack%data->getTimeConfiguration().getImagetteStackingNumber() != 0) m_imagettesPerStackedSubarray++;

	m_redundantHardware = data->redundantHardware();

}

void ImageWriter::process(Data * data, int timeStep, bool fullFrame) const {

	unsigned exposuresPerStack = data->getTimeConfiguration().getExposuresPerStack();

	if (fullFrame) {

		writeFullFrameImage(data,timeStep,true);
		data->clearImages();

	} else if ((timeStep+1)%exposuresPerStack == 0) { //do writing out each time there is a new stacked image

		//Require no stacking if dimensions correspond to full frame
		unsigned xdim = (*data->getImages().begin())->getXDim();
		unsigned ydim = (*data->getImages().begin())->getYDim();
		if (xdim==Image::kXDim && ydim==Image::kYDim && exposuresPerStack != 1) {
			throw runtime_error("Error in ImageWriter::process: Full frame images cannot be stacked. Either set the number of exposures per stacked image to 1 or reduce the sub-array dimensions.");
		}

		//Find out if any of the exposures making up the current stacked image are invalid
		//(target Earth occulted and/or spacecraft within SAA and/or stray light exceeding threshold, according to configuration)
		bool earthOccultation = false;
		bool saa = false;
		bool strayLight = false;
		for (int exposure = timeStep+1-exposuresPerStack; exposure<timeStep+1; exposure++) {
			SatelliteData * satData = data->getSatelliteData(exposure);
			if (data->omitEarthOccultation() && satData->getEarthOccultationFlag()) earthOccultation = true;
			if (data->omitSAA() && satData->getSAAFlag()) saa = true;
			if (data->omitStrayLight() && satData->getStrayLightFlag()) strayLight = true;
		}

		//Only write out images for which the stacked image is valid
		if ((timeStep+1)/exposuresPerStack <= data->getNumberOfCEsToOmitAtStart()) {

			cout << "omitting first " << data->getNumberOfCEsToOmitAtStart() << " stacked images" << endl;
			data->incrementDiscardedImages();

		} else if ((timeStep+1)/exposuresPerStack%data->getOnlyWriteEveryNthCE() != 0) {

			cout << "omitting stacked image because only writing of every " << data->getOnlyWriteEveryNthCE() << " stacked image has been requested" << endl;
			data->incrementDiscardedImages();

		} else if (earthOccultation || saa || strayLight) {

			cout << "omitting stacked image (" << (earthOccultation?" Earth occultation ":"") << (saa?" SAA ":"") << (strayLight?" Stray light ":"") << ")" << endl;
			data->incrementDiscardedImages();

		} else {

			if (data->getImages().size() != exposuresPerStack) {
				throw runtime_error("Error in ImageWriter::process: size of image vector does not match the number of images per stack");
			}

			if (m_writeStackedImages) {
				if (m_doublePrecisionStackedImages) {
					writeImages<SimRawDoubleprecisionsubarray>(data,timeStep,Data::stacked);
				} else  if (xdim==Image::kXDim && ydim==Image::kYDim) {
					writeFullFrameImage(data,timeStep,false);
				} else{
					writeImages<SciRawSubarray>(data,timeStep,Data::stacked);
				}
			}
			if (m_writeUnstackedImages) writeImages<SimRawUnstackedsubarray>(data,timeStep,Data::unstacked);
			if (m_writeImagettes && exposuresPerStack>1) writeImages<SciRawImagette>(data,timeStep,Data::imagette); //Only write imagettes if subarrays are stacked

			if (m_writeCentroid) {
				if (data->getFitsCentroid() == nullptr) initializeCentroid(data);
				writeCentroid(data,(timeStep+1)-exposuresPerStack);
			}

			data->incrementStackedImageCount();

		}

		data->clearImages();

	}

}


void ImageWriter::doEnd(Data* data) const {

	cout << "In total " << data->stackedImageCount() << " stacked were written out, and " << data->getNumberOfDiscardedImages() << " stacked images were discarded." << endl;

}

template <class T> void ImageWriter::writeImages(Data * data, int timeStep, Data::IMAGETYPE type) const {

	TimeConfiguration timeConf = data->getTimeConfiguration();
	unsigned exposuresPerStack = timeConf.getExposuresPerStack();
	unsigned stackedImageCount = data->stackedImageCount();
	unsigned unstackedImageCount = stackedImageCount*exposuresPerStack;
	unsigned imagetteCount = stackedImageCount*m_imagettesPerStackedSubarray;
	//cout << (timeStep+1)/exposuresPerStack  << " " << data->getNumberOfImagesPerCube(Data::stacked) << " " << ((timeStep+1)/exposuresPerStack)%data->getNumberOfImagesPerCube(Data::stacked) << endl;
	//cout << (timeStep+1) << " " << data->getNumberOfImagesPerCube(Data::unstacked) << " " << (timeStep+1)%data->getNumberOfImagesPerCube(Data::unstacked) << endl;

	string type_str;
	switch (type) {
	case Data::stacked:
		type_str = "stacked images";
		break;
	case Data::unstacked:
		type_str = "unstacked images";
		break;
	case Data::imagette:
		type_str = "imagettes";
		break;
	default:
		throw runtime_error("ERROR: Invalid data type in ImageWriter::writeImages");
	}

	//Start a new image cube if this is the first time writing out, or if the existing cube is full
	if ((type == Data::unstacked && unstackedImageCount%data->getNumberOfImagesPerCube(Data::unstacked) == 0) ||
		(type != Data::unstacked && stackedImageCount%data->getNumberOfImagesPerCube(Data::stacked) == 0)) {

		//Get the start time and number of layers for the new image cube

		double timeSinceStart = timeConf.getTimeSinceStart((timeStep+1)-exposuresPerStack);
		long seconds = static_cast<long>(lround((floor(timeSinceStart))));
		long milliseconds = static_cast<long>(lround((timeSinceStart - floor(timeSinceStart))*1000.));
		boost::posix_time::ptime startTime = timeConf.getStartTime() + boost::posix_time::seconds(seconds) + boost::posix_time::milliseconds(milliseconds);
		unsigned numberOfImages = data->getNumberOfImagesPerCube(type);

		//If the number of images remaining to be processed is fewer than the number of cube layers,
		//set the number of cube layers to the number of remaining images
		unsigned numberOfStackedImages = timeConf.getNumberOfStackedImagesPITL();
		unsigned numberOfTimeSteps = timeConf.getNumberOfTimeStepsPITL();
		if (type == Data::unstacked && (numberOfTimeSteps - unstackedImageCount) <
									  data->getNumberOfImagesPerCube(Data::unstacked)) {
			numberOfImages = numberOfTimeSteps - unstackedImageCount;
		} else if (type != Data::unstacked && (numberOfStackedImages - stackedImageCount) <
									  data->getNumberOfImagesPerCube(Data::stacked)) {
			unsigned remainingStackedImages = numberOfStackedImages - stackedImageCount;
			if (type==Data::imagette) {
				numberOfImages = remainingStackedImages * (exposuresPerStack/data->getTimeConfiguration().getImagetteStackingNumber());
				//If there is not a whole number of imagette stacks per sub-array stack, there is an additional partial imagette stack for each sub-array stack.
				if (exposuresPerStack % data->getTimeConfiguration().getImagetteStackingNumber() != 0) numberOfImages += remainingStackedImages;
			} else {
				numberOfImages = remainingStackedImages;
			}
		}

		//Initialize the new image cube
		cout << "Starting new image cube for " << type_str << " to contain " << numberOfImages << " images" << endl;
		data->resetFitsImageCube(type);
		initializeFitsOutput<T>(data,type,startTime,numberOfImages);
	}

	//Apply ADC saturation (max 2^16 ADU) to all images before stacking
	saturateImages(data,16);
	//Write the images to the current image cube
	unsigned numberOfWrittenStackedSubarrays = ((timeStep+1)/exposuresPerStack)-1;
	if (type==Data::stacked) {
		stackImages(data);
		writeImageData<T>(data, numberOfWrittenStackedSubarrays, stackedImageCount%data->getNumberOfImagesPerCube(type), type);
		if (m_writeTruthData) {
			if (data->doFrameTransferSmearing()) writeTruthSmearTrails(data, ((timeStep+1)/exposuresPerStack)-1, stackedImageCount%data->getNumberOfImagesPerCube(type));
			if (data->getCosmicRayImage() != nullptr) writeTruthCosmicImage(data, ((timeStep+1)/exposuresPerStack)-1, stackedImageCount%data->getNumberOfImagesPerCube(type));
		}
	} else if (type==Data::imagette) {
		extractImagettes(data);
		if (m_stackImagettes) stackImagettes(data);
		writeImageData<T>(data, numberOfWrittenStackedSubarrays*m_imagettesPerStackedSubarray, imagetteCount%data->getNumberOfImagesPerCube(type), type);
	} else {
		writeImageData<T>(data, numberOfWrittenStackedSubarrays*exposuresPerStack, unstackedImageCount%data->getNumberOfImagesPerCube(type), type);
	}
	//update the end time keywords
	setEndTimeKeywords<T>(data,timeStep,type);

}

template<class T> void ImageWriter::initializeFitsOutput(Data * data, Data::IMAGETYPE type, boost::posix_time::ptime startTime, unsigned numberOfImages) const {

	//Get mid time of first image
	UTC startTime_utc = UTC(to_iso_extended_string(startTime));
	UTC firstImageMidTime_utc = startTime_utc + DeltaTime(getImageDuration(data->getTimeConfiguration(),type)/2.);

	//Open image output file
	int xdim = type==Data::imagette ? m_imagetteSize : (*data->getImages().begin())->getXDim();
	int ydim = type==Data::imagette ? m_imagetteSize : (*data->getImages().begin())->getYDim();
	string subdir = type==Data::unstacked ? "/dfs/" : "/data/";
	T * imageCube = createFitsFile<T>(m_dirname+subdir, firstImageMidTime_utc,
			VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), {xdim,ydim,numberOfImages}, PassId());
	setImageKeywords<T>(imageCube,data,startTime,type);

	if (type==Data::stacked || type==Data::unstacked) {
		initializeCcdMargins(imageCube,data,xdim,ydim,startTime,numberOfImages,type);
	}

	if (m_writeTruthData && type==Data::stacked) {

		if (data->doFrameTransferSmearing()) {
			SimRawDoubleprecisionsubarray * imageCube_smear = createFitsFile<SimRawDoubleprecisionsubarray>(m_dirname+"/truth/", firstImageMidTime_utc,
					VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), {xdim,ydim,numberOfImages}, PassId(), "smear");
			setImageKeywords<SimRawDoubleprecisionsubarray>(imageCube_smear,data,startTime,type);
			data->setFitsImageCube_smear(imageCube_smear);
		}

		if (data->getCosmicRayImage() != nullptr) {
			SimRawDoubleprecisionsubarray * imageCube_cosmics = createFitsFile<SimRawDoubleprecisionsubarray>(m_dirname+"/truth/", firstImageMidTime_utc,
					VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), {xdim,ydim,numberOfImages}, PassId(), "cosmics");
			setImageKeywords<SimRawDoubleprecisionsubarray>(imageCube_cosmics,data,startTime,type);
			data->setFitsImageCube_cosmics(imageCube_cosmics);
		}

	}

	setFitsImageCubeWithMetaData<T>(data,type,imageCube);

	//Open the truth metadata fits file and define fits header keywords
	if ((type==Data::stacked || type==Data::unstacked || type==Data::imagette) && m_writeTruthData) {

		string type_str = "";
		if (type==Data::unstacked) {
			type_str = "Unstacked";
		}
		else if (type==Data::imagette) {
			type_str = "Imagette";
		}

		SimTruSubarray * truthMetaData = createFitsFile<SimTruSubarray>(m_dirname+"/truth/", firstImageMidTime_utc,
				                                                        VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), PassId(), type_str);
		setVisitKeywords(truthMetaData);
		setTargetKeywords(truthMetaData);
		truthMetaData->setKeyVStrtU(imageCube->getKeyVStrtU());
		truthMetaData->setKeyVStrtM(imageCube->getKeyVStrtM());
		data->setTruthMetaData(truthMetaData,type);

	}

}

template <class T> void ImageWriter::initializeCcdMargins(T * imageCube, Data * data, int xdim, int ydim, boost::posix_time::ptime startTime, unsigned numberOfImages, Data::IMAGETYPE type) const {

	SciRawBlankleft * blankLeftMarginCube = nullptr;
	SciRawBlankright * blankRightMarginCube = nullptr;
	SciRawDarkleft * darkLeftMarginCube = nullptr;
	SciRawDarkright * darkRightMarginCube = nullptr;
	SciRawDarktop * darkTopMarginCube = nullptr;
	SciRawOverscanleft * overscanLeftMarginCube = nullptr;
	SciRawOverscanright * overscanRightMarginCube = nullptr;
	SciRawOverscantop * overscanTopMarginCube = nullptr;

	bool fullFrame = (type==Data::fullframe);

	if (m_marginMode == image || fullFrame) {
		blankLeftMarginCube = Append_SciRawBlankleft(imageCube, {Image::kNBlankCols,ydim,numberOfImages});
		blankRightMarginCube = Append_SciRawBlankright(imageCube, {Image::kNBlankCols,ydim,numberOfImages});
		darkLeftMarginCube = Append_SciRawDarkleft(imageCube, {Image::kNDarkCols,ydim,numberOfImages});
		darkRightMarginCube = Append_SciRawDarkright(imageCube, {Image::kNDarkCols,ydim,numberOfImages});
		darkTopMarginCube = Append_SciRawDarktop(imageCube, {xdim,Image::kNDarkRows,numberOfImages});
		if (m_redundantHardware) {
			overscanRightMarginCube = Append_SciRawOverscanright(imageCube, {Image::kNOverscanCols,ydim,numberOfImages});
		} else {
			overscanLeftMarginCube = Append_SciRawOverscanleft(imageCube, {Image::kNOverscanCols,ydim,numberOfImages});
		}
		overscanTopMarginCube = Append_SciRawOverscantop(imageCube, {xdim,Image::kNOverscanRows,numberOfImages});
		setCcdMarginDescriptionKeywords<SciRawBlankleft>(blankLeftMarginCube,data,fullFrame);
		setCcdMarginDescriptionKeywords<SciRawBlankright>(blankRightMarginCube,data,fullFrame);
		setCcdMarginDescriptionKeywords<SciRawDarkleft>(darkLeftMarginCube,data,fullFrame);
		setCcdMarginDescriptionKeywords<SciRawDarkright>(darkRightMarginCube,data,fullFrame);
		setCcdMarginDescriptionKeywords<SciRawDarktop>(darkTopMarginCube,data,fullFrame);
		if (m_redundantHardware) {
			setCcdMarginDescriptionKeywords<SciRawOverscanright>(overscanRightMarginCube,data,fullFrame);
		} else {
			setCcdMarginDescriptionKeywords<SciRawOverscanleft>(overscanLeftMarginCube,data,fullFrame);
		}
		setCcdMarginDescriptionKeywords<SciRawOverscantop>(overscanTopMarginCube,data,fullFrame);
	} else if (m_marginMode == reduced) {
		blankLeftMarginCube = Append_SciRawBlankleft(imageCube, {4,1,numberOfImages});
		blankRightMarginCube = Append_SciRawBlankright(imageCube, {4,1,numberOfImages});
		darkLeftMarginCube = Append_SciRawDarkleft(imageCube, {3,ydim,numberOfImages});
		darkRightMarginCube = Append_SciRawDarkright(imageCube, {3,ydim,numberOfImages});
		darkTopMarginCube = Append_SciRawDarktop(imageCube, {xdim,Image::kNDarkRows,numberOfImages});
		if (m_redundantHardware) {
			overscanRightMarginCube = Append_SciRawOverscanright(imageCube, {4,1,numberOfImages});
		} else {
			overscanLeftMarginCube = Append_SciRawOverscanleft(imageCube, {4,1,numberOfImages});
		}
		overscanTopMarginCube = Append_SciRawOverscantop(imageCube, {xdim,Image::kNOverscanRows,numberOfImages});
		setCcdMarginDescriptionKeywords<SciRawBlankleft>(blankLeftMarginCube,data,fullFrame,4);
		setCcdMarginDescriptionKeywords<SciRawBlankright>(blankRightMarginCube,data,fullFrame,4);
		setCcdMarginDescriptionKeywords<SciRawDarkleft>(darkLeftMarginCube,data,fullFrame,3);
		setCcdMarginDescriptionKeywords<SciRawDarkright>(darkRightMarginCube,data,fullFrame,3);
		setCcdMarginDescriptionKeywords<SciRawDarktop>(darkTopMarginCube,data,fullFrame);
		if (m_redundantHardware) {
			setCcdMarginDescriptionKeywords<SciRawOverscanright>(overscanRightMarginCube,data,fullFrame,4);
		} else {
			setCcdMarginDescriptionKeywords<SciRawOverscanleft>(overscanLeftMarginCube,data,fullFrame,4);
		}
		setCcdMarginDescriptionKeywords<SciRawOverscantop>(overscanTopMarginCube,data,fullFrame);
	} else if (m_marginMode == total_collapsed) {
		blankLeftMarginCube = Append_SciRawBlankleft(imageCube, {4,1,numberOfImages});
		blankRightMarginCube = Append_SciRawBlankright(imageCube, {4,1,numberOfImages});
		darkLeftMarginCube = Append_SciRawDarkleft(imageCube, {4,1,numberOfImages});
		darkRightMarginCube = Append_SciRawDarkright(imageCube, {4,1,numberOfImages});
		darkTopMarginCube = Append_SciRawDarktop(imageCube, {1,4,numberOfImages});
		if (m_redundantHardware) {
			overscanRightMarginCube = Append_SciRawOverscanright(imageCube, {4,1,numberOfImages});
		} else {
			overscanLeftMarginCube = Append_SciRawOverscanleft(imageCube, {4,1,numberOfImages});
		}
		overscanTopMarginCube = Append_SciRawOverscantop(imageCube, {1,4,numberOfImages});
		setCcdMarginDescriptionKeywords<SciRawBlankleft>(blankLeftMarginCube,data,fullFrame,4);
		setCcdMarginDescriptionKeywords<SciRawBlankright>(blankRightMarginCube,data,fullFrame,4);
		setCcdMarginDescriptionKeywords<SciRawDarkleft>(darkLeftMarginCube,data,fullFrame,4);
		setCcdMarginDescriptionKeywords<SciRawDarkright>(darkRightMarginCube,data,fullFrame,4);
		setCcdMarginDescriptionKeywords<SciRawDarktop>(darkTopMarginCube,data,fullFrame,4);
		if (m_redundantHardware) {
			setCcdMarginDescriptionKeywords<SciRawOverscanright>(overscanRightMarginCube,data,fullFrame,4);
		} else {
			setCcdMarginDescriptionKeywords<SciRawOverscanleft>(overscanLeftMarginCube,data,fullFrame,4);
		}
		setCcdMarginDescriptionKeywords<SciRawOverscantop>(overscanTopMarginCube,data,fullFrame,4);
	}

	setCcdMarginColumnKeywords<SciRawBlankleft>(blankLeftMarginCube,data,startTime,type);
	setCcdMarginColumnKeywords<SciRawBlankright>(blankRightMarginCube,data,startTime,type);
	setCcdMarginColumnKeywords<SciRawDarkleft>(darkLeftMarginCube,data,startTime,type);
	setCcdMarginColumnKeywords<SciRawDarkright>(darkRightMarginCube,data,startTime,type);
	setCcdMarginRowKeywords<SciRawDarktop>(darkTopMarginCube,data,startTime,type);
	if (m_redundantHardware) {
		setCcdMarginColumnKeywords<SciRawOverscanright>(overscanRightMarginCube,data,startTime,type);
	} else {
		setCcdMarginColumnKeywords<SciRawOverscanleft>(overscanLeftMarginCube,data,startTime,type);
	}
	setCcdMarginRowKeywords<SciRawOverscantop>(overscanTopMarginCube,data,startTime,type);

	data->setFitsCcdMargins(blankLeftMarginCube,blankRightMarginCube,
			darkLeftMarginCube,darkRightMarginCube,darkTopMarginCube,
			overscanLeftMarginCube,overscanRightMarginCube,overscanTopMarginCube);

}

void ImageWriter::initializeCcdMargins(SimRawUnstackedsubarray * imageCube, Data * data, int xdim, int ydim, boost::posix_time::ptime startTime, unsigned numberOfImages, Data::IMAGETYPE type) const {

	SimRawUnstackedoverscanleftimage * overscanLeftImageCube = nullptr;
	SimRawUnstackedoverscanrightimage * overscanRightImageCube = nullptr;

	SimRawUnstackedblankleftimage * blankLeftImageCube = Append_SimRawUnstackedblankleftimage(imageCube, {xdim,ydim,numberOfImages});
	SimRawUnstackedblankrightimage * blankRightImageCube = Append_SimRawUnstackedblankrightimage(imageCube, {xdim,ydim,numberOfImages});
	SimRawUnstackeddarkleftimage * darkLeftImageCube = Append_SimRawUnstackeddarkleftimage(imageCube, {xdim,ydim,numberOfImages});
	SimRawUnstackeddarkrightimage * darkRightImageCube = Append_SimRawUnstackeddarkrightimage(imageCube, {xdim,ydim,numberOfImages});
	SimRawUnstackeddarktopimage * darkTopImageCube = Append_SimRawUnstackeddarktopimage(imageCube, {xdim,ydim,numberOfImages});
	if (m_redundantHardware) {
		overscanRightImageCube = Append_SimRawUnstackedoverscanrightimage(imageCube, {xdim,ydim,numberOfImages});
	} else {
		overscanLeftImageCube = Append_SimRawUnstackedoverscanleftimage(imageCube, {xdim,ydim,numberOfImages});
	}
	SimRawUnstackedoverscantopimage * overscanTopImageCube = Append_SimRawUnstackedoverscantopimage(imageCube, {xdim,ydim,numberOfImages});

	setCcdMarginColumnKeywords<SimRawUnstackedblankleftimage>(blankLeftImageCube,data,startTime,type);
	setCcdMarginColumnKeywords<SimRawUnstackedblankrightimage>(blankRightImageCube,data,startTime,type);
	setCcdMarginColumnKeywords<SimRawUnstackeddarkleftimage>(darkLeftImageCube,data,startTime,type);
	setCcdMarginColumnKeywords<SimRawUnstackeddarkrightimage>(darkRightImageCube,data,startTime,type);
	setCcdMarginRowKeywords<SimRawUnstackeddarktopimage>(darkTopImageCube,data,startTime,type);
	if (m_redundantHardware) {
		setCcdMarginColumnKeywords<SimRawUnstackedoverscanrightimage>(overscanRightImageCube,data,startTime,type);
	} else {
		setCcdMarginColumnKeywords<SimRawUnstackedoverscanleftimage>(overscanLeftImageCube,data,startTime,type);
	}
	setCcdMarginRowKeywords<SimRawUnstackedoverscantopimage>(overscanTopImageCube,data,startTime,type);

	data->setFitsCcdMarginImages(blankLeftImageCube,blankRightImageCube,
			darkLeftImageCube,darkRightImageCube,darkTopImageCube,
			overscanLeftImageCube,overscanRightImageCube,overscanTopImageCube);

}

void ImageWriter::initializeCentroid(Data *data) const {

	//Get the time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();

	//Get start and end times
	UTC startTime_utc = timeConf.getVisitStartTimeUTC() + DeltaTime(30.); //30s delay between start of visit and first image
	UTC endTime_utc = timeConf.getStartTimeUTC() + DeltaTime(timeConf.getDuration());

	//Open the centroid fits file and define fits header keywords
	SciRawCentroid * centroid = createFitsFile<SciRawCentroid>(m_dirname+"/data/", startTime_utc,
			VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), PassId());

	//Define fits header keywords for image
	setVisitKeywords(centroid);
	setTargetKeywords(centroid);
	centroid->setKeyProcChn("CHEOPSim");
	centroid->setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	centroid->setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
	centroid->setKeyVStrtU(startTime_utc);
	centroid->setKeyVStrtM(startTime_utc.getMjd());
	centroid->setKeyVStopU(endTime_utc);
	centroid->setKeyVStopM(endTime_utc.getMjd());

	//Add the centroid to the data
	data->setFitsCentroid(centroid);

}

template<class T> void ImageWriter::setImageKeywords(T * imageCube, Data * data, boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {

	//Define fits header keywords for image
	setCommonKeywords(imageCube,data,startTime,type);
	setTargetKeywords(imageCube);
	setImageOriginKeywords(imageCube,data);

}

/** *************************************************************************
 *  @brief Template specialization of setImageKeywords<class T>
 *  for T = SciRawFullarray
 */
template<> void ImageWriter::setImageKeywords(SciRawFullarray * imageCube, Data * data, boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {

	//Define fits header keywords for image
	setCommonKeywords(imageCube,data,startTime,type);
	setTargetKeywords(imageCube);
	setImageOriginKeywords(imageCube,data);
	setReadoutKeywords(imageCube,data,type);
	imageCube->setKeyMrgMode("image");

}

/** *************************************************************************
 *  @brief Template specialization of setImageKeywords<class T>
 *  for T = SciRawSubarray
 */
template<> void ImageWriter::setImageKeywords(SciRawSubarray * imageCube, Data * data, boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {

	//Define fits header keywords for image
	setCommonKeywords(imageCube,data,startTime,type);
	setTargetKeywords(imageCube);
	setImageOriginKeywords(imageCube,data);
	setReadoutKeywords(imageCube,data,type);
	imageCube->setKeyMrgMode(data->getVisit().m_marginMode);
	imageCube->setKeyStacking(data->getTimeConfiguration().getExposuresPerStack()>1 ? m_stackingMethod : "none");
	imageCube->setKeyShape(m_subarrayShape);

}

/** *************************************************************************
 *  @brief Template specialization of setImageKeywords<class T>
 *  for T = SimRawDoubleprecisionsubarray
 */
template<> void ImageWriter::setImageKeywords(SimRawDoubleprecisionsubarray * imageCube, Data * data, boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {

	//Define fits header keywords for image
	setCommonKeywords(imageCube,data,startTime,type);
	setTargetKeywords(imageCube);
	setImageOriginKeywords(imageCube,data);
	imageCube->setKeyShape(m_subarrayShape);

}

/** *************************************************************************
 *  @brief Template specialization of setImageKeywords<class T>
 *  for T = SciRawImagette
 */
template<> void ImageWriter::setImageKeywords(SciRawImagette * imageCube, Data * data, boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {

	TimeConfiguration timeConf = data->getTimeConfiguration();

	//Define fits header keywords for image
	setCommonKeywords(imageCube,data,startTime,type);
	setTargetKeywords(imageCube);
	imageCube->setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	imageCube->setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
	imageCube->setKeyStacking(timeConf.getImagetteStackingNumber()>1 ? m_stackingMethod : "none");
	imageCube->setKeyCropping((m_dynamicImagettes && !m_stackImagettes) ? "moving window" : "static window");
	imageCube->setKeyShape(m_imagetteShape);
	imageCube->setKeyRounding(0);
	imageCube->setKeyNlinCor(!data->gainNonLinearity()); // If non-linearity is simulated, set the non-linearity correction flag to false and vice versa
	imageCube->setKeyNexp(timeConf.getImagetteStackingNumber());
	imageCube->setKeyTexptime(timeConf.getExposureTimeAsDouble()*timeConf.getImagetteStackingNumber());

}


template<class T> void ImageWriter::setReadoutKeywords(T * imageCube, Data * data, Data::IMAGETYPE type) const {

	//Define fits header keywords for image
	imageCube->setKeyRounding(0);
	imageCube->setKeyNlinCor(!data->gainNonLinearity()); // If non-linearity is simulated, set the non-linearity correction flag to false and vice versa
	imageCube->setKeyRoScrpt(readoutScript(data->getVisit().m_readMode,data->getTimeConfiguration().getExposureTimeAsDouble(),type));
	imageCube->setKeyRoHw(m_redundantHardware ? "redundant" : "main");
	imageCube->setKeyRoFrequ(type==Data::fullframe ? 230000. : lround(data->getSerialReadRate()*1000.));

}

int ImageWriter::readoutScript(string readMode, double exposureTime, Data::IMAGETYPE type) const {

	// Define CCD readout timing script ID according to https://redmine.isdc.unige.ch/projects/cheops/wiki/CCD_Read-Out_Timing_Scripts
	int script_id;
	if (type==Data::fullframe) {
		if (exposureTime < 5.) {
			script_id = m_redundantHardware ? 1 : 2;
		} else {
			script_id = m_redundantHardware ? 3 : 4;
		}
	} else if (readMode == "bright") {
		script_id = m_redundantHardware ? 5 : 6;
	} else if (readMode == "ultrabright") {
		script_id = m_redundantHardware ? 7 : 8;
	} else if (readMode == "faint") {
		script_id = m_redundantHardware ? 32 : 31;
	} else if (readMode == "faint fast") {
		script_id = m_redundantHardware ? 34 : 33;
	}

	return script_id;

}

template<class T> void ImageWriter::setEndTimeKeywords(Data * data, int timeStep, Data::IMAGETYPE type) const {

	//Get end time and mid time of last image
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double timeSinceStart = timeConf.getTimeSinceStart(timeStep+1);
	long seconds = static_cast<long>(lround((floor(timeSinceStart))));
	long milliseconds = static_cast<long>(lround((timeSinceStart - floor(timeSinceStart))*1000.));
	boost::posix_time::ptime endTime = timeConf.getStartTime() + boost::posix_time::seconds(seconds) + boost::posix_time::milliseconds(milliseconds);
	UTC endTime_utc = UTC(to_iso_extended_string(endTime));
	UTC lastImageMidTime_utc;
	unsigned imagesPerStack = 1;
	if (type==Data::stacked) {
		imagesPerStack = timeConf.getExposuresPerStack();;
	} else if (type==Data::imagette && timeConf.getImagetteStackingNumber()>1) {
		imagesPerStack = timeConf.getImagetteStackingNumber();
	}
	lastImageMidTime_utc = endTime_utc - DeltaTime(imagesPerStack*timeConf.getRepetitionPeriod()) + DeltaTime(getImageDuration(timeConf,type)/2.);

	T * imageCube = nullptr;
	data->getFitsImageCube(&imageCube);
	if (imageCube == nullptr) {
		throw runtime_error("Error in ImageWriter::setEndTimeKeywords: image cube has not been initialized");
	}

	imageCube->setKeyVStopU(endTime_utc);
	imageCube->setKeyVStopM(endTime_utc.getMjd());
	imageCube->setKeyTStopU(lastImageMidTime_utc);
	imageCube->setKeyTStopM(lastImageMidTime_utc.getMjd());
	imageCube->setKeyTStopO(lastImageMidTime_utc.getObt());

	if (type==Data::imagette) {
		SciRawImagettemetadata * metaData = data->getFitsImagetteMetaData();
		metaData->setKeyVStopU(imageCube->getKeyVStopU());
		metaData->setKeyVStopM(imageCube->getKeyVStopM());
	} else {
		SciRawImagemetadata * metaData = data->getFitsImageMetaData(type);
		metaData->setKeyVStopU(imageCube->getKeyVStopU());
		metaData->setKeyVStopM(imageCube->getKeyVStopM());
	}
	if (type==Data::stacked && !m_doublePrecisionStackedImages) {
		SciRawUnstackedimagemetadata * metaData_unstacked = data->getFitsUnstackedImageMetaData(type);
		metaData_unstacked->setKeyVStopU(imageCube->getKeyVStopU());
		metaData_unstacked->setKeyVStopM(imageCube->getKeyVStopM());
	}
	if ((type==Data::stacked || type==Data::unstacked || type==Data::imagette) && m_writeTruthData) {
		SimTruSubarray * truthMetaData = data->getFitsTruthMetaData(type);
		truthMetaData->setKeyVStopU(imageCube->getKeyVStopU());
		truthMetaData->setKeyVStopM(imageCube->getKeyVStopM());
	}

	if (m_writeTruthData && type==Data::stacked) {

		if (data->doFrameTransferSmearing()) {
			SimRawDoubleprecisionsubarray * imageCube_smear = nullptr;
			data->getFitsImageCube_smear(&imageCube_smear);
			imageCube_smear->setKeyVStopU(endTime_utc);
			imageCube_smear->setKeyVStopM(endTime_utc.getMjd());
			imageCube_smear->setKeyTStopU(lastImageMidTime_utc);
			imageCube_smear->setKeyTStopM(lastImageMidTime_utc.getMjd());
			imageCube_smear->setKeyTStopO(lastImageMidTime_utc.getObt());
		}

		if (data->getCosmicRayImage() != nullptr) {
			SimRawDoubleprecisionsubarray * imageCube_cosmics = nullptr;
			data->getFitsImageCube_cosmics(&imageCube_cosmics);
			imageCube_cosmics->setKeyVStopU(endTime_utc);
			imageCube_cosmics->setKeyVStopM(endTime_utc.getMjd());
			imageCube_cosmics->setKeyTStopU(lastImageMidTime_utc);
			imageCube_cosmics->setKeyTStopM(lastImageMidTime_utc.getMjd());
			imageCube_cosmics->setKeyTStopO(lastImageMidTime_utc.getObt());
		}

	}

	if (type == Data::unstacked) {
		setCcdMarginImageEndTime<SimRawUnstackedblankleftimage>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SimRawUnstackedblankrightimage>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SimRawUnstackeddarkleftimage>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SimRawUnstackeddarkrightimage>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SimRawUnstackeddarktopimage>(data,endTime_utc,lastImageMidTime_utc);
		if (m_redundantHardware) {
			setCcdMarginImageEndTime<SimRawUnstackedoverscanrightimage>(data,endTime_utc,lastImageMidTime_utc);
		} else {
			setCcdMarginImageEndTime<SimRawUnstackedoverscanleftimage>(data,endTime_utc,lastImageMidTime_utc);
		}
		setCcdMarginImageEndTime<SimRawUnstackedoverscantopimage>(data,endTime_utc,lastImageMidTime_utc);
	} else if (type == Data::stacked) {
		setCcdMarginImageEndTime<SciRawBlankleft>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SciRawBlankright>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SciRawDarkleft>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SciRawDarkright>(data,endTime_utc,lastImageMidTime_utc);
		setCcdMarginImageEndTime<SciRawDarktop>(data,endTime_utc,lastImageMidTime_utc);
		if (m_redundantHardware) {
			setCcdMarginImageEndTime<SciRawOverscanright>(data,endTime_utc,lastImageMidTime_utc);
		} else {
			setCcdMarginImageEndTime<SciRawOverscanleft>(data,endTime_utc,lastImageMidTime_utc);
		}
		setCcdMarginImageEndTime<SciRawOverscantop>(data,endTime_utc,lastImageMidTime_utc);
	}

}

template<class T> void ImageWriter::setCcdMarginImageEndTime(Data *data, UTC endTime_utc, UTC lastImageMidTime_utc) const {
	T * imageCube = nullptr;
	data->getFitsCcdMargin(&imageCube);
	imageCube->setKeyVStopU(endTime_utc);
	imageCube->setKeyVStopM(endTime_utc.getMjd());
	imageCube->setKeyTStopU(lastImageMidTime_utc);
	imageCube->setKeyTStopM(lastImageMidTime_utc.getMjd());
	imageCube->setKeyTStopO(lastImageMidTime_utc.getObt());
}

template<class T> void ImageWriter::setImageOriginKeywords(T * imageCube, Data * data) const {
	imageCube->setKeyXWinoff((*data->getImages().begin())->getXOffset());
	imageCube->setKeyYWinoff((*data->getImages().begin())->getYOffset());
	imageCube->setKeyRaTarg(data->getFieldOfView()->getPointingDirection().getRightAscension()/3600.);
	imageCube->setKeyDecTarg(data->getFieldOfView()->getPointingDirection().getDeclination()/3600.);
}

template<class T> void ImageWriter::setCcdMarginRowKeywords(T * imageCube, Data * data,
		boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {
	setCommonKeywords(imageCube,data,startTime,type);
	imageCube->setKeyXWinoff((*data->getImages().begin())->getXOffset());
}

template<class T> void ImageWriter::setCcdMarginColumnKeywords(T * imageCube, Data * data,
		boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {
	setCommonKeywords(imageCube,data,startTime,type);
	imageCube->setKeyYWinoff((*data->getImages().begin())->getYOffset());
}

template<class T> void ImageWriter::setCcdMarginDescriptionKeywords(T * imageCube, Data *data, bool fullFrame, int nAverages) const {

	imageCube->setKeyStacking((!fullFrame && data->getTimeConfiguration().getExposuresPerStack()>1) ? m_marginStackingMethod : "none");
	imageCube->setKeyRounding(0);
	imageCube->setKeyNlinCor(!data->gainNonLinearity()); // If non-linearity is simulated, set the non-linearity correction flag to false and vice versa

	if (nAverages == 3) {
		long cubeSize[3];
		imageCube->GetSize(cubeSize);
		if (cubeSize[0] > cubeSize[1]) {
			imageCube->setKeyMrgProc("col collapsed");
		} else {
			imageCube->setKeyMrgProc("row collapsed");
		}
	} else if (nAverages == 4) {
		imageCube->setKeyMrgProc("total collapsed");
	} else {
		imageCube->setKeyMrgProc("image");
	}

	if (nAverages == 3 || nAverages == 4) {
		imageCube->setKeyMrgDty1("mean");
		imageCube->setKeyMrgDty2("stdev");
		imageCube->setKeyMrgDty3("median");
	}

	if (nAverages == 4) {
		imageCube->setKeyMrgDty4("mad");
	}

}

template<class T> void ImageWriter::setCommonKeywords(T * imageCube, Data * data,
		boost::posix_time::ptime startTime, Data::IMAGETYPE type) const {

	//Get time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double exposureTimeAsDouble = timeConf.getExposureTimeAsDouble();
	unsigned exposuresPerStack = type==Data::fullframe ? 1 : timeConf.getExposuresPerStack();

	//Get start time and mid time of first image
	UTC startTime_utc = UTC(to_iso_extended_string(startTime));
	UTC firstImageMidTime_utc = startTime_utc + DeltaTime(getImageDuration(timeConf,type)/2.);

	//Define fits header keywords for image
	setVisitKeywords(imageCube);
	imageCube->setKeyProcChn("CHEOPSim");
	imageCube->setKeyNexp(exposuresPerStack);
	imageCube->setKeyExptime(exposureTimeAsDouble);
	imageCube->setKeyTexptime(exposureTimeAsDouble*exposuresPerStack);
	imageCube->setKeyVStrtU(startTime_utc);
	imageCube->setKeyVStrtM(startTime_utc.getMjd());
	imageCube->setKeyTStrtU(firstImageMidTime_utc);
	imageCube->setKeyTStrtM(firstImageMidTime_utc.getMjd());
	imageCube->setKeyTStrtO(firstImageMidTime_utc.getObt());
}

template <class T> void ImageWriter::setVisitKeywords(T * fitsFile) const {
	//Define first set of fixed fits header keywords
	fitsFile->setKeyProcNum(m_visit.m_versionNum);
	fitsFile->setKeyArchRev(0);
	fitsFile->setKeyPiName(m_visit.m_piName);
	fitsFile->setKeyPiUid(m_visit.m_piUid);
	fitsFile->setKeyObsCat(m_visit.m_obsCat);
	fitsFile->setKeyProgtype(m_visit.m_progType);
	fitsFile->setKeyProgId(m_visit.m_progId);
	fitsFile->setKeyReqId(m_visit.m_reqId);
	fitsFile->setKeyVisitctr(m_visit.m_visitCtr);
	fitsFile->setKeyObsid(m_visit.m_obsId);
	fitsFile->setKeyPrpVst1(m_visit.m_prpFirst);
	fitsFile->setKeyPrpVstn(m_visit.m_prpLast);
}

template <class T> void ImageWriter::setTargetKeywords(T * fitsFile) const {
	//Define second set of fixed fits header keywords
	fitsFile->setKeyTargname(m_visit.m_targName);
	fitsFile->setKeySpectype(m_specType);
	fitsFile->setKeyTEff(m_TEff);
	fitsFile->setKeyMagG(m_Gmag);
	fitsFile->setKeyMagChps(m_cheopsMag);
	fitsFile->setKeyMagGerr(m_GmagErr);
	fitsFile->setKeyMagCerr(m_cheopsMagErr);
}

/** *************************************************************************
 *  @brief Template specialization of setFitsImageCubeWithMetaData<class T>
 *  for T = SciRawSubarray
 */
template<> void ImageWriter::setFitsImageCubeWithMetaData(Data * data, Data::IMAGETYPE type, SciRawSubarray* imageCube) const{

	SciRawImagemetadata * metaData = Append_SciRawImagemetadata(imageCube);
	setVisitKeywords(metaData);
	setTargetKeywords(metaData);
	metaData->setKeyVStrtU(imageCube->getKeyVStrtU());
	metaData->setKeyVStrtM(imageCube->getKeyVStrtM());
	metaData->setKeyProcChn(imageCube->getKeyProcChn()); //kept outside setVisitKeywords because not used for truth metadata
	metaData->setKeyRdMode(type==Data::fullframe ? "full frame" : data->getVisit().m_readMode);

	//Compression Entity keywords
	double repetitionPeriod = data->getTimeConfiguration().getRepetitionPeriod();
	metaData->setKeyAcqMode(1); //1: DUMP (sub-array) 2: DIGIT 3: FULL
	metaData->setKeyOversamp(true); //ToDo: Check this is correct. Documentation says: if true then averaging of several exposures is done
	metaData->setKeyFSource(0); //0: CCD 1: PATTERN 2:SIMULATION
	metaData->setKeyRepetit(repetitionPeriod);

	SciRawUnstackedimagemetadata * metaData_unstacked = Append_SciRawUnstackedimagemetadata(imageCube);
	setVisitKeywords(metaData_unstacked);
	setTargetKeywords(metaData_unstacked);
	metaData_unstacked->setKeyVStrtU(imageCube->getKeyVStrtU());
	metaData_unstacked->setKeyVStrtM(imageCube->getKeyVStrtM());
	metaData_unstacked->setKeyProcChn(imageCube->getKeyProcChn()); //kept outside setVisitKeywords because not used for truth metadata

	data->setFitsImageCube(imageCube, metaData, metaData_unstacked);

}

/** *************************************************************************
 *  @brief Template specialization of setFitsImageCubeWithMetaData<class T>
 *  for T = SimRawDoubleprecisionsubarray
 */
template<> void ImageWriter::setFitsImageCubeWithMetaData(Data * data, Data::IMAGETYPE type, SimRawDoubleprecisionsubarray* imageCube) const
{

	SciRawImagemetadata * metaData = Append_SciRawImagemetadata(imageCube);
	setVisitKeywords(metaData);
	setTargetKeywords(metaData);
	metaData->setKeyVStrtU(imageCube->getKeyVStrtU());
	metaData->setKeyVStrtM(imageCube->getKeyVStrtM());
	metaData->setKeyDataLvl("SIM");
	metaData->setKeyProcChn(imageCube->getKeyProcChn()); //kept outside setVisitKeywords because not used for truth metadata
	metaData->setKeyRdMode(type==Data::fullframe ? "full frame" : data->getVisit().m_readMode);

	//Compression Entity keywords
	double repetitionPeriod = data->getTimeConfiguration().getRepetitionPeriod();
	metaData->setKeyAcqMode(1); //1: DUMP (sub-array) 2: DIGIT 3: FULL
	metaData->setKeyOversamp(true); //ToDo: Check this is correct. Documentation says: if true then averaging of several exposures is done
	metaData->setKeyFSource(0); //0: CCD 1: PATTERN 2:SIMULATION
	metaData->setKeyRepetit(repetitionPeriod);

	data->setFitsImageCube(imageCube, metaData, nullptr);

}

/** *************************************************************************
 *  @brief Template specialization of setFitsImageCubeWithMetaData<class T>
 *  for T = SimRawUnstackedsubarray
 */
template<> void ImageWriter::setFitsImageCubeWithMetaData(Data * data, Data::IMAGETYPE type, SimRawUnstackedsubarray* imageCube) const
{

	SciRawImagemetadata * metaData = Append_SciRawImagemetadata(imageCube);
	setVisitKeywords(metaData);
	setTargetKeywords(metaData);
	metaData->setKeyVStrtU(imageCube->getKeyVStrtU());
	metaData->setKeyVStrtM(imageCube->getKeyVStrtM());
	metaData->setKeyDataLvl("SIM");
	metaData->setKeyProcChn(imageCube->getKeyProcChn()); //kept outside setVisitKeywords because not used for truth metadata
	metaData->setKeyRdMode(type==Data::fullframe ? "full frame" : data->getVisit().m_readMode);

	//Compression Entity keywords
	double repetitionPeriod = data->getTimeConfiguration().getRepetitionPeriod();
	metaData->setKeyAcqMode(1); //1: DUMP (sub-array) 2: DIGIT 3: FULL
	metaData->setKeyOversamp(true); //ToDo: Check this is correct. Documentation says: if true then averaging of several exposures is done
	metaData->setKeyFSource(0); //0: CCD 1: PATTERN 2:SIMULATION
	metaData->setKeyRepetit(repetitionPeriod);

	data->setFitsImageCube(imageCube, metaData);

}

/** *************************************************************************
 *  @brief Template specialization of setFitsImageCubeWithMetaData<class T>
 *  for T = SciRawImagette
 */
template<> void ImageWriter::setFitsImageCubeWithMetaData(Data * data, Data::IMAGETYPE type, SciRawImagette* imageCube) const {

	SciRawImagettemetadata * metaData = Append_SciRawImagettemetadata(imageCube);
	setVisitKeywords(metaData);
	setTargetKeywords(metaData);
	metaData->setKeyVStrtU(imageCube->getKeyVStrtU());
	metaData->setKeyVStrtM(imageCube->getKeyVStrtM());
	metaData->setKeyProcChn(imageCube->getKeyProcChn()); //kept outside setVisitKeywords because not used for truth metadata

	data->setFitsImageCube(imageCube, metaData);

}

void ImageWriter::saturateImages(Data* data, unsigned nBits) const {

	vector<Image*> images = data->getImages();
	for (vector<Image*>::const_iterator image = images.begin(); image!=images.end(); ++image) {
		(*image)->saturate(nBits);
	}

}

void ImageWriter::stackImages(Data* data) const {

	if (m_stackingMethod !="coadd" && m_stackingMethod !="mean") {
		throw runtime_error("Error in ImageWriter::stackImages: image stacking method must be coadd or mean");
	}

	if (m_marginStackingMethod !="coadd" && m_marginStackingMethod !="mean") {
		throw runtime_error("Error in ImageWriter::stackImages: margin stacking method must be coadd or mean");
	}

	vector<Image*> images = data->getImages();
	Image * firstImage = *(images.begin());
	Image * stackedImage = new Image(firstImage->getXDim(),firstImage->getYDim(),
									 firstImage->getXOffset(),firstImage->getYOffset());

	//Perform image coaddition
	for (vector<Image*>::const_iterator image = images.begin(); image!=images.end(); ++image) {
		if (!m_doublePrecisionStackedImages) (*image)->roundToIntegers();
		*stackedImage += **image;
	}

	//For mean rather than coadd (main image), divide each pixel of the main image by the number of coadded images
	if (m_stackingMethod=="mean") {
		for (int ix=0; ix<Image::kXDim; ix++) {
			for (int iy=0; iy<Image::kYDim; iy++) {
				stackedImage->setPixelValue(ix,iy,stackedImage->getPixelValue(ix,iy)/static_cast<double>(images.size()));
			}
		}
		if (!m_doublePrecisionStackedImages) stackedImage->roundToIntegers();
	}

	//For mean rather than coadd (margins), divide each margin pixel by the number of coadded images
	if (m_marginStackingMethod=="mean") {
		for (int ix=-Image::kLeftMargin; ix<Image::kXTotal-Image::kLeftMargin; ix++) {
			for (int iy=0; iy<Image::kYTotal; iy++) {
				if (ix<0 || ix>=Image::kYDim || iy>=Image::kYDim) {
					stackedImage->setPixelValue(ix,iy,stackedImage->getPixelValue(ix,iy)/static_cast<double>(images.size()));
				}
			}
		}
		if (!m_doublePrecisionStackedImages) stackedImage->roundToIntegers();

	}
	data->addStackedImage(stackedImage);

	//printTruthData(stackedImage->getTruthData());
}

void ImageWriter::extractImagettes(Data* data) const {

	vector<Image*> images = data->getImages();
	Data::SubarrayDimensions subarray = data->getSubarrayDimensions();

	for (vector<Image*>::const_iterator image = images.begin(); image!=images.end(); ++image) {

		pair<double,double> barycentre(0.,0.);
		//If dynamic imagettes is requested, get the PSF barycentre relative to the centre of the image, provided there is no imagette stacking
		if (m_dynamicImagettes && !m_stackImagettes) barycentre = getBarycentre(*image,subarray);

		//Define offset of edge of imagette w.r.t. full frame. The -0.5 is included for consistency with the calculation made by IFSW
		//(EngineeringAlgorithms getStartCoord and CompressionEntity::Compress)
		unsigned imageXOffset = lround(subarray.m_targetLocationX - 0.5 - 0.5*m_imagetteSize) + lround(barycentre.first);
		unsigned imageYOffset = lround(subarray.m_targetLocationY - 0.5 - 0.5*m_imagetteSize) + lround(barycentre.second);
		//cout << imageXOffset << " " << imageYOffset << " " << barycentre.first << " " << barycentre.second << endl;
		//cout << (512 + barycentre.first - imageXOffset - 17) << " " << (512 + barycentre.second - imageYOffset - 17) << endl;

		Image * imagette = new Image(m_imagetteSize,m_imagetteSize,imageXOffset,imageYOffset);

		//Fill fits image array
		for (unsigned j = 0; j<m_imagetteSize; j++) {
			for (unsigned i = 0; i<m_imagetteSize; i++) {
				imagette->setPixelValue(i,j,(*image)->getPixelValue(imageXOffset+i,imageYOffset+j));
			}
		}

		imagette->setTruthData((*image)->getTruthData());

		data->addImagette(imagette);

	}

}

void ImageWriter::stackImagettes(Data* data) const {

	vector<Image*> imagettes = data->getImagettes();
	Image * stackedImagette = nullptr;
	vector<unsigned> numberPerStack;
	unsigned n = 0;
	for (unsigned i=0; i<imagettes.size(); i++) {
		if (i%data->getTimeConfiguration().getImagetteStackingNumber() == 0) {
			if (i != 0) {
				data->addStackedImagette(stackedImagette);
				numberPerStack.push_back(n);
			}
			stackedImagette = new Image(m_imagetteSize,m_imagetteSize,imagettes.front()->getXOffset(),imagettes.front()->getYOffset());
			n = 0;
		}
		imagettes[i]->roundToIntegers();
		*stackedImagette += *(imagettes[i]);
		n++;
	}
	data->addStackedImagette(stackedImagette);
	numberPerStack.push_back(n);

//	cout << "numberPerStack: ";
//	for (unsigned i=0; i<numberPerStack.size(); i++) cout << numberPerStack[i] << " ";
//	cout << endl;

	vector<Image*> stackedImagettes = data->getStackedImagettes();
	for (unsigned i=0; i<stackedImagettes.size(); i++) {
		if (m_stackingMethod=="mean") {
			*stackedImagettes[i] *= (1./static_cast<double>(numberPerStack[i]));
			stackedImagettes[i]->roundToIntegers();
		}
	}


}

void ImageWriter::printTruthData(TruthData* truthData) const {

	for (unsigned i=0; i<truthData->getPSFs().size();  i++) {
		cout << "PSF: " << truthData->getPSFs()[i].getMeanXPosition() << " "
						<< truthData->getPSFs()[i].getMeanYPosition() << " "
						<< truthData->getPSFs()[i].getFlux() << endl;
		vector<double> xPositions = truthData->getPSFs()[i].getXPositions();
		vector<double> yPositions = truthData->getPSFs()[i].getYPositions();
		cout << "PSF positions for each exposure:" << endl;
		for (unsigned i=0; i<xPositions.size();  i++) {
			cout<< "(" << xPositions[i] << "," << yPositions[i] << ")" << endl;
		}
	}
	for (unsigned i=0; i<truthData->getDeadPixels().size();  i++) {
		cout << "dead pixel: " << truthData->getDeadPixels()[i].m_xPixel << " "
							  << truthData->getDeadPixels()[i].m_yPixel << " "
							  << truthData->getDeadPixels()[i].m_quantumEfficiency << endl;
	}
	for (unsigned i=0; i<truthData->getHotPixels().size();  i++) {
		cout << "hot pixel: " << truthData->getHotPixels()[i].m_xPixel << " "
							  << truthData->getHotPixels()[i].m_yPixel << " "
							  << truthData->getHotPixels()[i].m_nElectrons << " "
		  	  	  	  	  	  << truthData->getHotPixels()[i].m_type << endl;
	}

}

template <class T> void ImageWriter::writeImageData(Data * data, int imageCount, int cubeLayer, Data::IMAGETYPE type) const {

	//cout << "type: " << type << endl;

	//Get the images
	vector<Image*> images;
	switch (type) {
	case Data::stacked:
		images = data->getStackedImages();
		break;
	case Data::imagette:
		if (m_stackImagettes) {
			images = data->getStackedImagettes();
		} else {
			images = data->getImagettes();
		}
		break;
	default:
		images = data->getImages();
		break;
	}

	//Get the image cube
	T * imageCube = nullptr;
	data->getFitsImageCube(&imageCube);
	if (imageCube == nullptr) throw runtime_error("Error in ImageWriter::writeImageData: image cube has not been initialized");

	//Loop over images
	for (vector<Image*>::const_iterator image = images.begin(); image!=images.end(); ++image) {

		//cout << "cubeLayer = " << cubeLayer << endl;

		Averages blankLeftAverages;
		Averages blankRightAverages;
		Averages darkLeftAverages;
		Averages darkRightAverages;
		Averages darkTopAverages;
		Averages overscanLeftAverages;
		Averages overscanTopAverages;

		int imageXOffset = 0;
		int imageYOffset = 0;
		unsigned imageXSize = 0;
		unsigned imageYSize = 0;

		if (type==Data::imagette) {

			imageXOffset = (*image)->getXOffset();
			imageYOffset = (*image)->getYOffset();
			imageXSize = (*image)->getXDim();
			imageYSize = (*image)->getYDim();

			SciRawImagettemetadata * metaData = data->getFitsImagetteMetaData();
			if (metaData == nullptr) throw runtime_error("Error in ImageWriter::writeImageData: metadata has not been initialized");
			metaData->setCellXOffFullArray(imageXOffset);
			metaData->setCellYOffFullArray(imageYOffset);
			metaData->setCellXOffSubArray(imageXOffset-data->getSubarrayDimensions().m_xOffset);
			metaData->setCellYOffSubArray(imageYOffset-data->getSubarrayDimensions().m_yOffset);

			//Fill fits image array
			for (unsigned j = 0; j<m_imagetteSize; j++) {
				for (unsigned i = 0; i<m_imagetteSize; i++) {
					double pixelValue = (*image)->getPixelValue(i,j);
					if (m_imagetteShape=="circular") {
						int x = i - m_imagetteSize/2;
						int y = j - m_imagetteSize/2;
						if (sqrt(x*x + y*y) > m_imagetteSize/2) pixelValue = 0.;
					}
					((*imageCube)[cubeLayer])[j][i] = pixelValue>0. ? static_cast<uint32_t>(lround(pixelValue)+0.5) : 0;
				}
			}

		} else {

			imageXOffset = (*image)->getXOffset();
			imageYOffset = (*image)->getYOffset();
			imageXSize = (*image)->getXDim();
			imageYSize = (*image)->getYDim();

			array2D * overscanLeftImage = new array2D(boost::extents[Image::kNOverscanCols][(*image)->getYDim()]);
			array2D * blankLeftImage = new array2D(boost::extents[Image::kNBlankCols][(*image)->getYDim()]);
			array2D * darkLeftImage = new array2D(boost::extents[Image::kNDarkCols][(*image)->getYDim()]);
			array2D * blankRightImage = new array2D(boost::extents[Image::kNBlankCols][(*image)->getYDim()]);
			array2D * darkRightImage = new array2D(boost::extents[Image::kNDarkCols][(*image)->getYDim()]);
			array2D * darkTopImage = new array2D(boost::extents[(*image)->getXDim()][Image::kNDarkRows]);
			array2D * overscanTopImage = new array2D(boost::extents[(*image)->getXDim()][Image::kNOverscanRows]);

			//Fill fits image array

			for (int j = 0; j<Image::kYTotal; j++) {
				for (int i = 0; i<Image::kXTotal; i++) {

					int ix = i-Image::kLeftMargin;
					int iy = j;
					double pixelValue = (*image)->getPixelValue(ix,iy);

					if (iy>=imageYOffset && iy<imageYOffset+(*image)->getYDim() &&
						ix>=imageXOffset && ix<imageXOffset+(*image)->getXDim()) {

						if (m_subarrayShape=="circular" && type != Data::unstacked) {
							double radius = min((*image)->getXDim(),(*image)->getYDim())/2.;
							int x = ix - (imageXOffset + (*image)->getXDim()/2);
							int y = iy - (imageYOffset + (*image)->getYDim()/2);
							if (sqrt(x*x + y*y) > radius) pixelValue = 0.;
						}
						writePixelData<T>(imageCube,cubeLayer,ix-imageXOffset,iy-imageYOffset,type,pixelValue);

					} else {

						if (iy>=imageYOffset && iy<imageYOffset+(*image)->getYDim()) {

							if (i<Image::kNOverscanCols) {
								(*overscanLeftImage)[i][j-imageYOffset] = pixelValue;
							} else if (i < Image::kNOverscanCols+Image::kNBlankCols) {
								(*blankLeftImage)[i-Image::kNOverscanCols][j-imageYOffset] = pixelValue;
							} else if (i < Image::kLeftMargin) {
								(*darkLeftImage)[i-Image::kNOverscanCols-Image::kNBlankCols][j-imageYOffset] = pixelValue;
							} else if (ix>=Image::kXDim && ix<Image::kXDim+Image::kNDarkCols) {
								(*darkRightImage)[ix-Image::kXDim][j-imageYOffset] = pixelValue;
							} else if (ix>=Image::kXDim+Image::kNDarkCols) {
								(*blankRightImage)[ix-Image::kXDim-Image::kNDarkCols][j-imageYOffset] = pixelValue;
							}

						} else if (ix>=imageXOffset && ix<imageXOffset+(*image)->getXDim()) {

							if (iy>=Image::kYDim && iy<Image::kYDim+Image::kNDarkRows) {
								(*darkTopImage)[ix-imageXOffset][j-Image::kYDim] = pixelValue;
							} else if (iy>=Image::kYDim+Image::kNDarkRows) {
								(*overscanTopImage)[ix-imageXOffset][j-Image::kYDim-Image::kNDarkRows] = pixelValue;
							}

						}

					}

				}
			}

			blankLeftAverages = ccdMarginAverages(blankLeftImage);
			blankRightAverages = ccdMarginAverages(blankRightImage);
			darkLeftAverages = ccdMarginAverages(darkLeftImage);
			darkRightAverages = ccdMarginAverages(darkRightImage);
			darkTopAverages = ccdMarginAverages(darkTopImage);
			overscanLeftAverages = ccdMarginAverages(overscanLeftImage);
			overscanTopAverages = ccdMarginAverages(overscanTopImage);

			if (type == Data::unstacked) {
				writeCcdMarginImage<SimRawUnstackedblankleftimage>(data,blankLeftImage,cubeLayer,type);
				writeCcdMarginImage<SimRawUnstackedblankrightimage>(data,blankRightImage,cubeLayer,type);
				writeCcdMarginImage<SimRawUnstackeddarkleftimage>(data,darkLeftImage,cubeLayer,type);
				writeCcdMarginImage<SimRawUnstackeddarkrightimage>(data,darkRightImage,cubeLayer,type);
				writeCcdMarginImage<SimRawUnstackeddarktopimage>(data,darkTopImage,cubeLayer,type);
				if (m_redundantHardware) {
					//For the redundant hardware case, write the overscan data to SCI_RAW_OverscanRight (the overscan data is in the left-most
					//column of the Image array, hence overscanLeftImage, but we write it to SCI_RAW_OverscanRight)
					writeCcdMarginImage<SimRawUnstackedoverscanrightimage>(data,overscanLeftImage,cubeLayer,type);
				} else {
					writeCcdMarginImage<SimRawUnstackedoverscanleftimage>(data,overscanLeftImage,cubeLayer,type);
				}
				writeCcdMarginImage<SimRawUnstackedoverscantopimage>(data,overscanTopImage,cubeLayer,type);
			} else {
				if (m_marginMode == reduced) {
					writeCcdMargin<SciRawBlankleft>(data,blankLeftImage,blankLeftAverages,cubeLayer,type,total_collapsed);
					writeCcdMargin<SciRawBlankright>(data,blankRightImage,blankRightAverages,cubeLayer,type,total_collapsed);
					writeCcdMargin<SciRawDarkleft>(data,darkLeftImage,darkLeftAverages,cubeLayer,type,reduced);
					writeCcdMargin<SciRawDarkright>(data,darkRightImage,darkRightAverages,cubeLayer,type,reduced);
					writeCcdMargin<SciRawDarktop>(data,darkTopImage,darkTopAverages,cubeLayer,type,ImageWriter::image);
					if (m_redundantHardware) {
						//For the redundant hardware case, write the overscan data to SCI_RAW_OverscanRight (the overscan data is in the left-most
						//column of the Image array, hence overscanLeftImage,overscanLeftAverages, but we write it to SCI_RAW_OverscanRight)
						writeCcdMargin<SciRawOverscanright>(data,overscanLeftImage,overscanLeftAverages,cubeLayer,type,total_collapsed);
					} else {
						writeCcdMargin<SciRawOverscanleft>(data,overscanLeftImage,overscanLeftAverages,cubeLayer,type,total_collapsed);
					}
					writeCcdMargin<SciRawOverscantop>(data,overscanTopImage,overscanTopAverages,cubeLayer,type,ImageWriter::image);
				} else {
					writeCcdMargin<SciRawBlankleft>(data,blankLeftImage,blankLeftAverages,cubeLayer,type,m_marginMode);
					writeCcdMargin<SciRawBlankright>(data,blankRightImage,blankRightAverages,cubeLayer,type,m_marginMode);
					writeCcdMargin<SciRawDarkleft>(data,darkLeftImage,darkLeftAverages,cubeLayer,type,m_marginMode);
					writeCcdMargin<SciRawDarkright>(data,darkRightImage,darkRightAverages,cubeLayer,type,m_marginMode);
					writeCcdMargin<SciRawDarktop>(data,darkTopImage,darkTopAverages,cubeLayer,type,m_marginMode);
					if (m_redundantHardware) {
						//For the redundant hardware case, write the overscan data to SCI_RAW_OverscanRight (the overscan data is in the left-most
						//column of the Image array, hence overscanLeftImage,overscanLeftAverages, but we write it to SCI_RAW_OverscanRight)
						writeCcdMargin<SciRawOverscanright>(data,overscanLeftImage,overscanLeftAverages,cubeLayer,type,m_marginMode);
					} else {
						writeCcdMargin<SciRawOverscanleft>(data,overscanLeftImage,overscanLeftAverages,cubeLayer,type,m_marginMode);
					}
					writeCcdMargin<SciRawOverscantop>(data,overscanTopImage,overscanTopAverages,cubeLayer,type,m_marginMode);
				}
			}

			delete overscanLeftImage;
			delete blankLeftImage;
			delete darkLeftImage;
			delete blankRightImage;
			delete darkRightImage;
			delete darkTopImage;
			delete overscanTopImage;

		}

		//Fill the image meta data
		writeMetaData(data,type,imageCount);

		//Fill truth meta data
		if ((type==Data::stacked || type==Data::unstacked || type==Data::imagette) && m_writeTruthData) {
			writeTruthMetaData(data,(*image)->getTruthData(),imageXOffset,imageYOffset,imageXSize,imageYSize,type,imageCount);
		}

		imageCount += 1;
		cubeLayer += 1;

	}

}

pair<double,double> ImageWriter::getBarycentre(const Image * image, const Data::SubarrayDimensions & subarray) const {

	//Determine PSF barycentre using either truth information or photometry
	pair<double,double> barycentre;
	if (m_truthBarycentre) {
		vector<PSF> truthPSFs = image->getTruthData()->getPSFs();
		if (truthPSFs.size()>0) {
			barycentre = make_pair(truthPSFs[0].getMeanXPosition() - subarray.m_targetLocationX,
								   truthPSFs[0].getMeanYPosition() - subarray.m_targetLocationY);
		} else {
			throw runtime_error("Error in ImageWriter::getBarycentre: no PSF truth data available for centroid determination.");
		}
	} else {
		barycentre = m_photometry->getBarycentre(image);
	}

	return barycentre;

}

void ImageWriter::writeMetaData(Data * data, Data::IMAGETYPE type, int imageCount) const {

	//Get time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();
	unsigned exposuresPerImage,remainder;
	double startTimeOffset;
	switch(type) {
	case Data::stacked:
		exposuresPerImage = timeConf.getExposuresPerStack();
		startTimeOffset = imageCount*exposuresPerImage*timeConf.getRepetitionPeriod();
		break;
	case Data::imagette:
		exposuresPerImage = timeConf.getImagetteStackingNumber();
		//For imagettes the start time offset is determined according to the number of exposures in previously processed stacked sub-arrays
		//plus the number of imagette stacks corresponding to the current stacked sub-array
		startTimeOffset = ((imageCount/m_imagettesPerStackedSubarray)*timeConf.getExposuresPerStack() +  (imageCount%m_imagettesPerStackedSubarray)*exposuresPerImage)*timeConf.getRepetitionPeriod();
		//If there is not a whole number of imagette stacks per sub-array stack,
		//the last imagette stack corresponding to the current stacked sub-array is a partial imagette stack containing the remainder
		remainder = timeConf.getExposuresPerStack()%timeConf.getImagetteStackingNumber();
		if (remainder != 0 && imageCount%m_imagettesPerStackedSubarray==m_imagettesPerStackedSubarray-1) exposuresPerImage = remainder;
		//cout << imageCount/m_imagettesPerStackedSubarray << " " << timeConf.getExposuresPerStack() << " " << imageCount%m_imagettesPerStackedSubarray << " " << exposuresPerImage << " " << timeConf.getRepetitionPeriod() << endl;
		break;
	default:
		exposuresPerImage = 1;
		startTimeOffset = imageCount*exposuresPerImage*timeConf.getRepetitionPeriod();
		break;
	}

	//Get time information for meta data
	UTC imageMidTime_utc = UTC(to_iso_extended_string(timeConf.getStartTime())) + DeltaTime(startTimeOffset + getImageDuration(timeConf,type)/2.);
	//cout << type << " " << to_iso_extended_string(timeConf.getStartTime()) << " " << startTimeOffset << " " << getImageDuration(timeConf,type)/2. << " " << imageMidTime_utc.getUtc() << endl;

	//Concatenate the satellite data corresponding to the image duration
	vector<SatelliteData*> satData;
	if (data->hasSatelliteData()) {
		if (type==Data::stacked) {
			for (unsigned i=imageCount*exposuresPerImage; i<(imageCount+1)*exposuresPerImage; i++) {
				satData.push_back(data->getSatelliteData(i));
			}
		} else {
			satData.push_back(data->getSatelliteData(imageCount));
		}
	}

	//Calculate the mean telescope and DPU temperatures, sun, moon and earth angles, and spacecraft latitude and longitude,
	//averaged over all time steps making up the image
	double telescopeTemperature = 0.;
	double dpuTemperature = 0.;
	double sunAngle = 0.;
	double moonAngle = 0.;
	double earthAngle = 0.;
	double latitude = 0.;
	double longitude = 0.;
	if (satData.size() != 0) {
		int timeStep = imageCount*exposuresPerImage;
		for (vector<SatelliteData*>::const_iterator it = satData.begin(); it!=satData.end(); ++it) {
			telescopeTemperature += (*it)->getTelescopeTemperature();
			dpuTemperature += (*it)->getDpuTemperature();
			sunAngle += (*it)->getSunAngle();
			moonAngle += (*it)->getMoonAngle();
			earthAngle += (*it)->getEarthLimbAngle();
			latitude += (*it)->getLatitude();
			longitude += (*it)->getLongitude();
			timeStep++;
		}
		telescopeTemperature/=exposuresPerImage;
		dpuTemperature/=exposuresPerImage;
		sunAngle/=exposuresPerImage;
		moonAngle/=exposuresPerImage;
		earthAngle/=exposuresPerImage;
		latitude/=exposuresPerImage;
		longitude/=exposuresPerImage;
	} else {
		telescopeTemperature = SatelliteData::kDefaultTelescopeTemperature;
		dpuTemperature = SatelliteData::kDefaultDpuTemperature;
	}

	//Get the metadata
	if (type==Data::imagette) {

		SciRawImagettemetadata * metaData = data->getFitsImagetteMetaData();
		if (metaData == nullptr) throw runtime_error("Error in ImageWriter::writeMetaData: metadata has not been initialized");

		//Fill metadata table
		metaData->setCellUtcTime(imageMidTime_utc);
		metaData->setCellMjdTime(imageMidTime_utc.getMjd());
		metaData->setCellObtTime(imageMidTime_utc.getObt());
		metaData->setCellNexp(exposuresPerImage);
		//CE counter for sub-array images starts from 1
		//unless counter=1 is already taken by a full frame image in which case the counter starts from 2
		int indexOffset = m_ceCounterOffset + (data->doFullFrame() ? 2 : 1);
		metaData->setCellCeCounter(data->stackedImageCount()+data->getNumberOfDiscardedImages()+indexOffset);

		metaData->WriteRow();

	} else {

		//Get HK data
		Data::HKData hkData;
		if (type==Data::unstacked) {
			//For unstacked images, use the raw HK data closest in time to the midpoint of the exposure
			hkData = data->getClosestRawHKData((imageMidTime_utc-timeConf.getVisitStartTimeUTC()).getSeconds());
		} else {
			//For stacked images and full frame, use the averaged HK data closest in time to 10s after the midpoint of the combined exposure (since the HK value at a given time is the average over the previous 20s)
			hkData = data->getClosestAveragedHKData((imageMidTime_utc-timeConf.getVisitStartTimeUTC()).getSeconds() + 10.);
		}

		SciRawImagemetadata * metaData = data->getFitsImageMetaData(type);
		if (metaData == nullptr) throw runtime_error("Error in ImageWriter::writeMetaData: metadata has not been initialized");

		//Fill metadata table
		metaData->setCellUtcTime(imageMidTime_utc);
		metaData->setCellMjdTime(imageMidTime_utc.getMjd());
		metaData->setCellObtTime(imageMidTime_utc.getObt());
		if (satData.size() != 0) {
			metaData->setCellLosToSunAngle(sunAngle);
			metaData->setCellLosToMoonAngle(moonAngle);
			metaData->setCellLosToEarthAngle(earthAngle);
			metaData->setCellLatitude(latitude);
			metaData->setCellLongitude(longitude);
		} else {
			metaData->setNullLosToSunAngle();
			metaData->setNullLosToMoonAngle();
			metaData->setNullLosToEarthAngle();
			metaData->setNullLatitude();
			metaData->setNullLongitude();
		}
		//CE counter for sub-array images starts from 1
		//unless counter=1 is already taken by a full frame image in which case the counter starts from 2
		int indexOffset = m_ceCounterOffset + (data->doFullFrame() ? 2 : 1);
		metaData->setCellCeCounter(data->stackedImageCount()+data->getNumberOfDiscardedImages()+indexOffset);
		metaData->setCellCeIntegrity(0);
		metaData->setCellPixDataOffset(static_cast<uint16_t>(lround(data->getBiasOffset())+0.5));
		metaData->setCellCcdTimingScript(readoutScript(data->getVisit().m_readMode,data->getTimeConfiguration().getExposureTimeAsDouble(),type));
		metaData->setCellHkSource(type==Data::unstacked ? "ce" : "hk tm");
		metaData->setCellHkTempFeeCcd(hkData.m_ccdTemp-273.15);
		metaData->setCellHkTempFeeBias(hkData.m_biasTemp-273.15);
		metaData->setCellHkTempFeeAdc(hkData.m_adcTemp-273.15);
		metaData->setCellAdcTemp1(dpuTemperature-273.15);
		metaData->setCellAdcN5v(SatelliteData::kDefaultDpuVoltage);
		metaData->setCellThermaft4(telescopeTemperature-273.15+0.4+0.5);
		metaData->setCellThermaft3(telescopeTemperature-273.15+0.3+0.5);
		metaData->setCellThermaft2(telescopeTemperature-273.15+0.2+0.5);
		metaData->setCellThermaft1(telescopeTemperature-273.15+0.1+0.5);
		metaData->setCellThermfront1(telescopeTemperature-273.15+0.5);
		metaData->setCellThermfront2(telescopeTemperature-273.15-0.1+0.5);
		metaData->setCellThermfront3(telescopeTemperature-273.15-0.2+0.5);
		metaData->setCellThermfront4(telescopeTemperature-273.15-0.3+0.5);
		metaData->setCellHkVoltFeeVod(hkData.m_vod);
		metaData->setCellHkVoltFeeVrd(hkData.m_vrd);
		metaData->setCellHkVoltFeeVog(hkData.m_vog);
		metaData->setCellHkVoltFeeVss(hkData.m_vss);
		metaData->setCellLeftDarkColMask(65535);
		metaData->setCellRightDarkColMask(65535);

		metaData->WriteRow();

		//Stacked images have an additional extension SCI_RAW_UnstackedImageMetadata to store some information with the cadence of unstacked images
		if (type==Data::stacked && !m_doublePrecisionStackedImages) {

			SciRawUnstackedimagemetadata * metaData_unstacked = data->getFitsUnstackedImageMetaData(type);
			if (metaData_unstacked == nullptr) throw runtime_error("Error in ImageWriter::writeMetaData: unstacked metadata for stacked images has not been initialized");

			//Loop over the unstacked images
			unsigned iexp=0;
			vector<Image*> images = data->getImages();
			for (vector<Image*>::const_iterator image = images.begin(); image!=images.end(); ++image) {

				//get the mid-time of the current unstacked image
				imageMidTime_utc = UTC(to_iso_extended_string(timeConf.getStartTime())) + DeltaTime(startTimeOffset + timeConf.getRepetitionPeriod()*iexp + getImageDuration(timeConf,Data::unstacked)*0.5);
				metaData_unstacked->setCellUtcTime(imageMidTime_utc);
				metaData_unstacked->setCellMjdTime(imageMidTime_utc.getMjd());
				metaData_unstacked->setCellObtTime(imageMidTime_utc.getObt());

				//CE counter for sub-array images starts from 1
				//unless counter=1 is already taken by a full frame image in which case the counter starts from 2
				int indexOffset = m_ceCounterOffset + (data->doFullFrame() ? 2 : 1);
				metaData_unstacked->setCellCeCounter(data->stackedImageCount()+data->getNumberOfDiscardedImages()+indexOffset);

				metaData_unstacked->setCellGain0(m_nominalGain);
				metaData_unstacked->setCellBias0(data->getBiasOffset());

				//Calculate the bias as the mean of the left overscan column for the current unstacked image
				double bias=0;
				for (int j=(*image)->getYOffset(); j<(*image)->getYOffset()+(*image)->getYDim(); j++) {
					for (int i = 0; i<Image::kNOverscanCols; i++) {
						bias += (*image)->getPixelValue(i-Image::kLeftMargin,j);
					}
				}
				bias /= ((*image)->getYDim()*Image::kNOverscanCols);
				metaData_unstacked->setCellBias(bias);

				//For unstacked images, use the raw HK data closest in time to the midpoint of the exposure
				hkData = data->getClosestRawHKData((imageMidTime_utc-timeConf.getVisitStartTimeUTC()).getSeconds());
				metaData_unstacked->setCellCeTempFeeCcd(hkData.m_ccdTemp-273.15);
				metaData_unstacked->setCellCeVoltFeeVod(hkData.m_vod);
				metaData_unstacked->setCellCeVoltFeeVrd(hkData.m_vrd);
				metaData_unstacked->setCellCeVoltFeeVog(hkData.m_vog);
				metaData_unstacked->setCellCeVoltFeeVss(hkData.m_vss);

				metaData_unstacked->WriteRow();

				iexp++;

			}

		}

	}

}

void ImageWriter::writeCentroid(Data * data, int timeStep) const {

	//Get time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double repetitionPeriod = timeConf.getRepetitionPeriod();

	UTC startTime_utc = UTC(to_iso_extended_string(timeConf.getStartTime()));

	SciRawCentroid * centroid = data->getFitsCentroid();
	if (centroid == nullptr) {
		throw runtime_error("Error in ImageWriter::writeCentroid: centroid fits file has not been initialized");
	}

	//Loop over images
	vector<Image*> images = data->getImages();
	for (vector<Image*>::const_iterator image = images.begin(); image!=images.end(); ++image) {

		//Get the subarray dimensions (includes intended target location)
		Data::SubarrayDimensions subarray = data->getSubarrayDimensions();

		//Get the PSF barycentre, relative to the centre of the image
		pair<double,double> barycentre = getBarycentre(*image,subarray);

		//Get integration start time and end time
		UTC imageStartTime_utc = startTime_utc + DeltaTime(timeStep*repetitionPeriod);
		UTC imageEndTime_utc = imageStartTime_utc + DeltaTime(timeConf.getExposureTimeAsDouble());

		//Fill centroid table
		centroid->setCellUtcStart(imageStartTime_utc);
		centroid->setCellMjdStart(imageStartTime_utc.getMjd());
		centroid->setCellObtStart(imageStartTime_utc.getObt());
		centroid->setCellUtcStop(imageEndTime_utc);
		centroid->setCellMjdStop(imageEndTime_utc.getMjd());
		centroid->setCellObtStop(imageEndTime_utc.getObt());
		centroid->setCellFullFrame(false);
		//CE counter for sub-array images starts from 1
		//unless counter=1 is already taken by a full frame image in which case the counter starts from 2
		int indexOffset = m_ceCounterOffset + (data->doFullFrame() ? 2 : 1);
		centroid->setCellCeCounter(data->stackedImageCount()+data->getNumberOfDiscardedImages()+indexOffset);
		if (m_targetLocationFromJitter) {
			centroid->setCellOffsetX(0.);
			centroid->setCellOffsetY(0.);
			centroid->setCellLocationX(subarray.m_targetLocationX + barycentre.first);
			centroid->setCellLocationY(subarray.m_targetLocationY + barycentre.second);
		} else {
			centroid->setCellOffsetX(barycentre.first);
			centroid->setCellOffsetY(barycentre.second);
			centroid->setCellLocationX(subarray.m_targetLocationX);
			centroid->setCellLocationY(subarray.m_targetLocationY);
		}
		centroid->setCellDataCadence(repetitionPeriod);
		centroid->setCellValidity(0);

		centroid->WriteRow();

		timeStep += 1;

	}

}

void ImageWriter::writeTruthMetaData(Data * data, TruthData * truthData, int imageXOffset, int imageYOffset,
									 int imageXSize, int imageYSize, Data::IMAGETYPE type, int imageCount) const {

	//Get time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();
	unsigned exposuresPerImage,remainder;
	double startTimeOffset;
	switch(type) {
	case Data::stacked:
		exposuresPerImage = timeConf.getExposuresPerStack();
		startTimeOffset = imageCount*exposuresPerImage*timeConf.getRepetitionPeriod();
		break;
	case Data::imagette:
		exposuresPerImage = timeConf.getImagetteStackingNumber();
		startTimeOffset = ((imageCount/m_imagettesPerStackedSubarray)*timeConf.getExposuresPerStack() +  (imageCount%m_imagettesPerStackedSubarray)*exposuresPerImage)*timeConf.getRepetitionPeriod();
		remainder = timeConf.getExposuresPerStack()%timeConf.getImagetteStackingNumber();
		if (remainder != 0 && imageCount%m_imagettesPerStackedSubarray==m_imagettesPerStackedSubarray-1) exposuresPerImage = remainder;
		break;
	default:
		exposuresPerImage = 1;
		startTimeOffset = imageCount*exposuresPerImage*timeConf.getRepetitionPeriod();
		break;
	}

	//Get time information for meta data
	UTC imageMidTime_utc = UTC(to_iso_extended_string(timeConf.getStartTime())) + DeltaTime(startTimeOffset + getImageDuration(timeConf,type)/2.);

	//Concatenate the satellite data corresponding to the image duration
	vector<SatelliteData*> satData;
	if (data->hasSatelliteData()) {
		if (type==Data::stacked) {
			for (unsigned i=imageCount*exposuresPerImage; i<(imageCount+1)*exposuresPerImage; i++) {
				satData.push_back(data->getSatelliteData(i));
			}
		} else {
			satData.push_back(data->getSatelliteData(imageCount));
		}
	}

	//Set the jitter AOCS and science validity flags
	bool validAocs = true;
	bool validScience = true;
	if (satData.size() != 0) {
		for (vector<SatelliteData*>::const_iterator it = satData.begin(); it!=satData.end(); ++it) {
			if (!(*it)->validAocs()) validAocs = false;
			if (!(*it)->validScience()) validScience = false;
		}
	}

	SimTruSubarray * truthMetaData = data->getFitsTruthMetaData(type);

	truthMetaData->setCellUtcTime(imageMidTime_utc);
	truthMetaData->setCellMjdTime(imageMidTime_utc.getMjd());
	truthMetaData->setCellObtTime(imageMidTime_utc.getObt());
	truthMetaData->setCellValidAocs(validAocs);
	truthMetaData->setCellValidScience(validScience);
	truthMetaData->setCellFullWellSaturated(truthData->isFullWellSaturated());
	truthMetaData->setCellAdcSaturated(truthData->isAdcSaturated());
	truthMetaData->setCellGlobalThroughput(truthData->getGlobalThroughput());
	truthMetaData->setCellGain(truthData->getGain());
	truthMetaData->setCellStrayLight(truthData->getStrayLight());
	truthMetaData->setCellZodiacalLight(truthData->getZodiacalLight());

	fillTruthData<SimTruSubarray>(truthMetaData,truthData,satData,imageXOffset,imageYOffset,imageXSize,imageYSize,false,type==Data::imagette);

}

template <class T> void ImageWriter::writeCcdMargin(Data * data, array2D * image, ImageWriter::Averages averages, int cubeLayer, Data::IMAGETYPE type, MARGIN_MODE marginMode) const {

	T * imageCube = nullptr;
	data->getFitsCcdMargin(&imageCube);
	if (imageCube == nullptr) {
		throw runtime_error("Error in ImageWriter::writeCcdMargin: image cube has not been initialized");
	}

	if (marginMode == ImageWriter::image) {

		for (unsigned j=0; j<image->shape()[1]; j++) {
			for (unsigned i=0; i<image->shape()[0]; i++) {
				double pixelValue = (*image)[i][j];
				writePixelData<T>(imageCube,cubeLayer,i,j,type,pixelValue,false);
			}
		}

	} else if (marginMode == ImageWriter::reduced) {

		int length = image->shape()[0];//row
		if (image->shape()[1] > image->shape()[0]) {//column
			length = image->shape()[1];
		}

		for (int i=0; i<length; i++) {
			averages = ccdMarginAverages(image,i);
			if (image->shape()[1] > image->shape()[0]) { //column
				writePixelData<T>(imageCube,cubeLayer,0,i,type,averages.mean,false);
				writePixelData<T>(imageCube,cubeLayer,1,i,type,averages.stddev,false);
				writePixelData<T>(imageCube,cubeLayer,2,i,type,averages.median,false);
			} else { //row
				writePixelData<T>(imageCube,cubeLayer,i,0,type,averages.mean,false);
				writePixelData<T>(imageCube,cubeLayer,i,1,type,averages.stddev,false);
				writePixelData<T>(imageCube,cubeLayer,i,2,type,averages.median,false);
			}
		}

	} else {

		if (image->shape()[1] > image->shape()[0]) { //column
			writePixelData<T>(imageCube,cubeLayer,0,0,type,averages.mean,false);
			writePixelData<T>(imageCube,cubeLayer,1,0,type,averages.stddev,false);
			writePixelData<T>(imageCube,cubeLayer,2,0,type,averages.median,false);
			writePixelData<T>(imageCube,cubeLayer,3,0,type,averages.mad,false);
		} else { //row
			writePixelData<T>(imageCube,cubeLayer,0,0,type,averages.mean,false);
			writePixelData<T>(imageCube,cubeLayer,0,1,type,averages.stddev,false);
			writePixelData<T>(imageCube,cubeLayer,0,2,type,averages.median,false);
			writePixelData<T>(imageCube,cubeLayer,0,3,type,averages.mad,false);
		}

	}

}

template <class T> void ImageWriter::writeCcdMarginImage(Data * data, array2D * image, int cubeLayer, Data::IMAGETYPE type) const {

	T * imageCube = nullptr;
	data->getFitsCcdMargin(&imageCube);
	if (imageCube == nullptr) {
		throw runtime_error("Error in ImageWriter::writeCcdMarginImage: image cube has not been initialized");
	}

	for (unsigned j=0; j<image->shape()[1]; j++) {
		for (unsigned i=0; i<image->shape()[0]; i++) {
			double pixelValue = (*image)[i][j];
			writePixelData<T>(imageCube,cubeLayer,i,j,type,pixelValue);
		}
	}

}

ImageWriter::Averages ImageWriter::ccdMarginAverages(array2D * image, int rowcol_index) const {

	accumulator_set<double, stats<tag::mean, tag::variance, tag::count> > acc;
	vector<double> pixelValues, absoluteDeviations;

	bool isColumn = true;
	int length = image->shape()[1];
	int width = image->shape()[0];
	if (width > length) {
		isColumn = false;
		length = image->shape()[0];
		width = image->shape()[1];
	}

	int jmin = rowcol_index;
	int jmax = rowcol_index+1;
	if (rowcol_index == -1) { //total collapsed mode: calculate average over length of the margin
		jmin = 0;
		jmax = length;
	}

	//accumulate statistics
	for (int j=jmin; j<jmax; j++) {
		for (int i=0; i<width; i++) {
			double pixelValue = isColumn ? (*image)[i][j] : (*image)[j][i];
			//cout << i << " " << j << " " << pixelValue << endl;
			acc(pixelValue);
			pixelValues.push_back(pixelValue);
		}
	}

	//calculate mean, median, stddev and mad collapsed to a single value
	double mean = boost::accumulators::mean(acc);
	int n = boost::accumulators::count(acc);
	double stddev = sqrt(variance(acc)*n/(n-1));
	double median = CommonTools::median(pixelValues);
	//cout << "median: " << median << " " << n << endl << endl << endl;

	//calculate mean absolute deviation (MAD) collapsed to a single column/row
	for (int j=jmin; j<jmax; j++) {
		for (int i=0; i<width; i++) {
			double pixelValue = isColumn ? (*image)[i][j] : (*image)[j][i];
			absoluteDeviations.push_back(abs(pixelValue-median));
			//cout << "absoluteDeviations: " << abs(pixelValue-median) << endl;
		}
	}
	double mad = CommonTools::median(absoluteDeviations);
	//cout << "mad: " << mad << endl;

	return Averages(mean,median,stddev,mad);

}

template <class T> void ImageWriter::writePixelData(T *imageCube, int cubeLayer, int i, int j, Data::IMAGETYPE type, double pixelValue, bool round) const {

	if (type==Data::unstacked) {
		((*imageCube)[cubeLayer])[j][i] = pixelValue>0. ? static_cast<uint16_t>(lround(pixelValue)+0.5) : 0;
	} else {
		if (m_doublePrecisionStackedImages || !round) {
			((*imageCube)[cubeLayer])[j][i] = pixelValue;
		} else {
			((*imageCube)[cubeLayer])[j][i] = pixelValue>0. ? static_cast<uint32_t>(lround(pixelValue)+0.5) : 0;
		}
	}

}

template <class T> void ImageWriter::fillTruthData(T * truthMetaData, TruthData * truthData, vector<SatelliteData*> satData,
												   int imageXOffset, int imageYOffset, unsigned imageXSize, unsigned imageYSize,
												   bool fullFrame, bool imagette) const {

	//Fill the mean roll angle over the time interval corresponding to the image
	vector<double> rollAngles;
	for (vector<SatelliteData*>::const_iterator it = satData.begin(); it!=satData.end(); ++it) {
		vector<double> rollAngles_temp = (*it)->getRollAngles();
		rollAngles.insert(rollAngles.end(),rollAngles_temp.begin(),rollAngles_temp.end());
	}
	//Avoid 0-360 wrap around effect (assumes roll angle value should always increase with time, except for 0-360 wrap around)
	for (unsigned i=1; i< rollAngles.size(); i++) {
		if (rollAngles[i]<rollAngles[i-1]) rollAngles[i]+=360;
	}
	float meanRollAngle = accumulate(rollAngles.begin(), rollAngles.end(), 0.)/rollAngles.size();
	if (meanRollAngle > 360.) meanRollAngle -= 360.;
	truthMetaData->setCellRollAngle(meanRollAngle);

	//Fill the truth PSF positions for the target star
	vector<PSF> psfs = truthData->getPSFs();
	float targetPsfX[kJitterTruthSize], targetPsfY[kJitterTruthSize];
	vector<double> targetPsfXpositions, targetPsfYpositions;
	if (psfs.size() > 0) {
		targetPsfXpositions = psfs[0].getXPositions();
		targetPsfYpositions = psfs[0].getYPositions();
	}
	for (unsigned i=0; i<kJitterTruthSize; i++) {
		if (i<targetPsfXpositions.size()) {
			targetPsfX[i] = targetPsfXpositions[i] - imageXOffset;
			targetPsfY[i] = targetPsfYpositions[i] - imageYOffset;
		} else {
			targetPsfX[i] = kNullFloat;
			targetPsfY[i] = kNullFloat;
		}
	}
	truthMetaData->setCellTargetPsfX(targetPsfX);
	truthMetaData->setCellTargetPsfY(targetPsfY);
	for (unsigned i=0; i<kJitterTruthSize; i++) {
		if (targetPsfX[i] == kNullFloat) {
			truthMetaData->setNullTargetPsfX(i);
			truthMetaData->setNullTargetPsfY(i);
		}
	}

	//Fill the PSF truth information for each star
	float psfMeanX[kPsfTruthSize], psfMeanY[kPsfTruthSize], psfFlux[kPsfTruthSize];
	for (unsigned i=0; i<kPsfTruthSize; i++) {
		if (i<psfs.size() && psfs[i].getFlux()>0.) { //out of frame PSFs have flux set to zero in PSFGenerator
			psfMeanX[i] = psfs[i].getMeanXPosition() - imageXOffset;
			psfMeanY[i] = psfs[i].getMeanYPosition() - imageYOffset;
			psfFlux[i] = psfs[i].getFlux();
		} else {
			psfMeanX[i] = kNullFloat;
			psfMeanY[i] = kNullFloat;
			psfFlux[i] = kNullFloat;
		}
	}
	truthMetaData->setCellPsfMeanX(psfMeanX);
	truthMetaData->setCellPsfMeanY(psfMeanY);
	truthMetaData->setCellPsfFlux(psfFlux);
	for (unsigned i=0; i<kPsfTruthSize; i++) {
		if (psfMeanX[i] == kNullFloat) {
			truthMetaData->setNullPsfMeanX(i);
			truthMetaData->setNullPsfMeanY(i);
			truthMetaData->setNullPsfFlux(i);
		}
	}

	//Fill the Cosmic truth information
	unsigned cosmicTruthSize = fullFrame ? kCosmicTruthSize_fullFrame : kCosmicTruthSize;
	vector<TruthData::CosmicPixel> cosmicPixels = truthData->getCosmicPixels();
	int32_t cosmicXpixel[cosmicTruthSize], cosmicYpixel[cosmicTruthSize], cosmicNelectrons[cosmicTruthSize];
	for (unsigned i=0; i<cosmicTruthSize; i++) {
		cosmicXpixel[i] = kNullInt;
		cosmicYpixel[i] = kNullInt;
		cosmicNelectrons[i] = kNullInt;
		if (i<cosmicPixels.size()) {
			int ix = cosmicPixels[i].m_xPixel;
			int iy = cosmicPixels[i].m_yPixel;
			//Shift x and y values from full frame to sub-array coordinates
			float xpos = ix - imageXOffset;
			float ypos = iy - imageYOffset;
			if (!imagette) {
				//Account for CCD margins, effectively placing them adjacent to the sub-array
				if (ix<0) xpos = ix;
				if (ix>=Image::kXDim) xpos = ix - Image::kXDim + imageXSize;
				if (iy>=Image::kYDim) ypos = iy - Image::kYDim + imageYSize;
			}
			if (!imagette || (xpos>=0 && xpos<imageXSize && ypos>=0 && ypos<imageYSize)) {
				cosmicXpixel[i] = xpos;
				cosmicYpixel[i] = ypos;
				cosmicNelectrons[i] = cosmicPixels[i].m_nElectrons;
			}
		}
	}
	truthMetaData->setCellCosmicXpixel(cosmicXpixel);
	truthMetaData->setCellCosmicYpixel(cosmicYpixel);
	truthMetaData->setCellCosmicNelectrons(cosmicNelectrons);
	for (unsigned i=0; i<cosmicTruthSize; i++) {
		if (cosmicXpixel[i] == kNullInt) {
			truthMetaData->setNullCosmicXpixel(i);
			truthMetaData->setNullCosmicYpixel(i);
			truthMetaData->setNullCosmicNelectrons(i);
		}
	}

	//Fill the hot pixel information
	unsigned hotTruthSize = fullFrame ? kHotTruthSize_fullFrame : kHotTruthSize;
	vector<TruthData::HotPixel> hotPixels = truthData->getHotPixels();
	int32_t hotXpixel[hotTruthSize], hotYpixel[hotTruthSize], hotNelectrons[hotTruthSize], hotType[hotTruthSize];
	for (unsigned i=0; i<hotTruthSize; i++) {
		hotXpixel[i] = kNullInt;
		hotYpixel[i] = kNullInt;
		hotNelectrons[i] = kNullInt;
		hotType[i] = kNullInt;
		if (i<hotPixels.size()) {
			int ix = hotPixels[i].m_xPixel;
			int iy = hotPixels[i].m_yPixel;
			//Shift x and y values from full frame to sub-array coordinates
			float xpos = ix - imageXOffset;
			float ypos = iy - imageYOffset;
			if (!imagette) {
				//Account for CCD margins, effectively placing them adjacent to the sub-array
				if (ix<0) xpos = ix;
				if (ix>=Image::kXDim) xpos = ix - Image::kXDim + imageXSize;
				if (iy>=Image::kYDim) ypos = iy - Image::kYDim + imageYSize;
			}
			if (!imagette || (xpos>=0 && xpos<imageXSize && ypos>=0 && ypos<imageYSize)) {
				hotXpixel[i] = xpos;
				hotYpixel[i] = ypos;
				hotNelectrons[i] = hotPixels[i].m_nElectrons;
				hotType[i] = static_cast<int32_t>(hotPixels[i].m_type);
			}
		}
	}
	truthMetaData->setCellHotXpixel(hotXpixel);
	truthMetaData->setCellHotYpixel(hotYpixel);
	truthMetaData->setCellHotNelectrons(hotNelectrons);
	truthMetaData->setCellHotType(hotType);
	for (unsigned i=0; i<hotTruthSize; i++) {
		if (hotXpixel[i] == kNullInt) {
			truthMetaData->setNullHotXpixel(i);
			truthMetaData->setNullHotYpixel(i);
			truthMetaData->setNullHotNelectrons(i);
			truthMetaData->setNullHotType(i);
		}
	}

	//Fill the dead pixel information
	unsigned deadTruthSize = fullFrame ? kDeadTruthSize_fullFrame : kDeadTruthSize;
	vector<TruthData::DeadPixel> deadPixels = truthData->getDeadPixels();
	int32_t deadXpixel[deadTruthSize], deadYpixel[deadTruthSize];
	float deadQE[deadTruthSize];
	for (unsigned i=0; i<deadTruthSize; i++) {
		deadXpixel[i] = kNullInt;
		deadYpixel[i] = kNullInt;
		deadQE[i] = kNullFloat;
		if (i<deadPixels.size()) {
			float xpos = deadPixels[i].m_xPixel - imageXOffset;
			float ypos = deadPixels[i].m_yPixel - imageYOffset;
			if (xpos>=0 && xpos<imageXSize && ypos>=0 && ypos<imageYSize) {
				deadXpixel[i] = xpos;
				deadYpixel[i] = ypos;
				deadQE[i] = deadPixels[i].m_quantumEfficiency;
			}
		}
	}
	truthMetaData->setCellDeadXpixel(deadXpixel);
	truthMetaData->setCellDeadYpixel(deadYpixel);
	truthMetaData->setCellDeadQe(deadQE);
	for (unsigned i=0; i<deadTruthSize; i++) {
		if (deadXpixel[i] == kNullInt) {
			truthMetaData->setNullDeadXpixel(i);
			truthMetaData->setNullDeadYpixel(i);
			truthMetaData->setNullDeadQe(i);
		}
	}

	//Fill the frame transfer smear trail information
	unsigned smearRowLength = fullFrame ? Image::kXDim : imageXSize;
	vector<double> smearTrailRow_vec = truthData->getSmearTrailRow();
	float smearTrailRow[smearRowLength];
	for (unsigned i=0; i<smearRowLength; i++) {
		if (smearTrailRow_vec.size() == smearRowLength) {//Imagettes are the only case where this is not true
			smearTrailRow[i] = smearTrailRow_vec[i];
		} else {
			smearTrailRow[i] = kNullFloat;
		}
	}
	truthMetaData->setCellSmearRow(smearTrailRow);

	truthMetaData->WriteRow();

}

void ImageWriter::writeTruthSmearTrails(Data * data, int imageCount, int cubeLayer) const {

	//Stack the images for each exposure containing the up and down smear trails
	vector<Image*> imagesToSmearUp = data->getImagesToSmearUp();
	vector<Image*> imagesToSmearDown = data->getImagesToSmearDown();
	Image * firstImage = *(imagesToSmearUp.begin());
	int imageXOffset = firstImage->getXOffset();
	int imageYOffset = firstImage->getYOffset();
	int imageXSize = firstImage->getXDim();
	int imageYSize = firstImage->getYDim();
	Image * stackedImage_smear = new Image(imageXSize,imageYSize,imageXOffset,imageYOffset);
	for (vector<Image*>::const_iterator image = imagesToSmearUp.begin(); image!=imagesToSmearUp.end(); ++image) {
		*stackedImage_smear += **image;
	}
	for (vector<Image*>::const_iterator image = imagesToSmearDown.begin(); image!=imagesToSmearDown.end(); ++image) {
		*stackedImage_smear += **image;
	}

	//Get the image cube
	SimRawDoubleprecisionsubarray * imageCube = nullptr;
	data->getFitsImageCube_smear(&imageCube);
	if (imageCube == nullptr) throw runtime_error("Error in ImageWriter::writeTruthSmearTrails: image cube has not been initialized");
	imageCube->setKeyBunit("Electrons");

	//Get the sub-array dimensions
	Image * subarray = *data->getImages().begin();
	imageXOffset = subarray->getXOffset();
	imageYOffset = subarray->getYOffset();
	imageXSize = subarray->getXDim();
	imageYSize = subarray->getYDim();

	//Write the data
	for (int j = 0; j<Image::kYDim; j++) {
		for (int i = 0; i<Image::kXDim; i++) {

			if (j>=imageYOffset && j<imageYOffset+imageYSize &&
				i>=imageXOffset && i<imageXOffset+imageXSize) {

				//cout << j-imageYOffset << " " << i-imageXOffset << " " << stackedImage_smear->getPixelValue(i,j) << endl;
				((*imageCube)[cubeLayer])[j-imageYOffset][i-imageXOffset] = stackedImage_smear->getPixelValue(i,j);

			}
		}
	}

	delete stackedImage_smear;

}

void ImageWriter::writeTruthCosmicImage(Data * data, int imageCount, int cubeLayer) const {

	//Get the image cube
	SimRawDoubleprecisionsubarray * imageCube = nullptr;
	data->getFitsImageCube_cosmics(&imageCube);
	if (imageCube == nullptr) throw runtime_error("Error in ImageWriter::writeTruthCosmicImage: image cube has not been initialized");
	imageCube->setKeyBunit("Electrons");

	//Get the cosmic ray truth image
	Image * image = data->getCosmicRayImage();

	//Write the data
	for (int j = 0; j<Image::kYDim; j++) {
		for (int i = 0; i<Image::kXDim; i++) {

			if (j>=image->getYOffset() && j<image->getYOffset()+image->getYDim() &&
				i>=image->getXOffset() && i<image->getXOffset()+image->getXDim()) {

				//cout << j-imageYOffset << " " << i-imageXOffset << " " << image->getPixelValue(i,j) << endl;
				((*imageCube)[cubeLayer])[j-image->getYOffset()][i-image->getXOffset()] = image->getPixelValue(i,j);

			}
		}
	}

}

void ImageWriter::writeFullFrameImage(Data * data, int timeStep, bool initialFullFrame) const {

	//Get time configuration
	TimeConfiguration timeConf = data->getTimeConfiguration();
	double imageDuration = timeConf.getExposureTimeAsDouble();
	double timeSinceStart = timeConf.getTimeSinceStart(timeStep);
	if (initialFullFrame) timeSinceStart = -imageDuration -5.; //Start of initial full frame image, with 5s gap between end of image and sub-array sequence start time
	long seconds = static_cast<long>(lround((trunc(timeSinceStart))));
	long milliseconds = static_cast<long>(lround((timeSinceStart - trunc(timeSinceStart))*1000.));
	boost::posix_time::ptime startTime = timeConf.getStartTime() + boost::posix_time::seconds(seconds) + boost::posix_time::milliseconds(milliseconds);

	//Get start time, mid time and end time of the image
	UTC startTime_utc = UTC(to_iso_extended_string(startTime));
	UTC midTime_utc = startTime_utc + DeltaTime(imageDuration/2.);
	UTC endTime_utc = startTime_utc + DeltaTime(imageDuration);

	//Open output file
	if (!boost::filesystem::exists(m_dirname+"/data")) boost::filesystem::create_directory(m_dirname+"/data");
	SciRawFullarray * fullFrameImage = createFitsFile<SciRawFullarray>(m_dirname+"/data/", midTime_utc,
			                           VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), PassId());
	setImageKeywords<SciRawFullarray>(fullFrameImage,data,startTime,Data::fullframe);

	//(Re-)initialize the margin data structures
	data->resetFitsImageCube(Data::stacked);
	initializeCcdMargins(fullFrameImage,data,Image::kXDim,Image::kYDim,startTime,1,Data::fullframe);

    //Set end time keywords
	fullFrameImage->setKeyVStopU(endTime_utc);
	fullFrameImage->setKeyVStopM(endTime_utc.getMjd());
	fullFrameImage->setKeyTStopU(midTime_utc);
	fullFrameImage->setKeyTStopM(midTime_utc.getMjd());
	fullFrameImage->setKeyTStopO(midTime_utc.getObt());
	setCcdMarginImageEndTime<SciRawBlankleft>(data,endTime_utc,midTime_utc);
	setCcdMarginImageEndTime<SciRawBlankright>(data,endTime_utc,midTime_utc);
	setCcdMarginImageEndTime<SciRawDarkleft>(data,endTime_utc,midTime_utc);
	setCcdMarginImageEndTime<SciRawDarkright>(data,endTime_utc,midTime_utc);
	setCcdMarginImageEndTime<SciRawDarktop>(data,endTime_utc,midTime_utc);
	if (m_redundantHardware) {
		setCcdMarginImageEndTime<SciRawOverscanright>(data,endTime_utc,midTime_utc);
	} else {
		setCcdMarginImageEndTime<SciRawOverscanleft>(data,endTime_utc,midTime_utc);
	}
	setCcdMarginImageEndTime<SciRawOverscantop>(data,endTime_utc,midTime_utc);

	//Get the last unstacked image in memory (check that there is only one)
	if (data->getImages().size() != 1) throw runtime_error("Error in ImageWriter::writeFullFrameImage: There are multiple full frame images in memory.");
	Image * image = data->getImages().back();

	//Apply ADC saturation
	image->saturate(16);

	array2D * overscanLeftImage = new array2D(boost::extents[Image::kNOverscanCols][image->getYDim()]);
	array2D * blankLeftImage = new array2D(boost::extents[Image::kNBlankCols][image->getYDim()]);
	array2D * darkLeftImage = new array2D(boost::extents[Image::kNDarkCols][image->getYDim()]);
	array2D * blankRightImage = new array2D(boost::extents[Image::kNBlankCols][image->getYDim()]);
	array2D * darkRightImage = new array2D(boost::extents[Image::kNDarkCols][image->getYDim()]);
	array2D * darkTopImage = new array2D(boost::extents[image->getXDim()][Image::kNDarkRows]);
	array2D * overscanTopImage = new array2D(boost::extents[image->getXDim()][Image::kNOverscanRows]);

	//Fill fits image array
	for (int j = 0; j<Image::kYTotal; j++) {
		for (int i = 0; i<Image::kXTotal; i++) {

			int ix = i-Image::kLeftMargin;
			int iy = j;
			double pixelValue = image->getPixelValue(ix,iy);

			if (iy>=0 && iy<Image::kYDim &&
				ix>=0 && ix<Image::kXDim) {

				(*fullFrameImage)[iy][ix] = pixelValue>0. ? static_cast<uint16_t>(lround(pixelValue)+0.5) : 0;

			} else {

				if (iy>=0 && iy<Image::kYDim) {

					if (i<Image::kNOverscanCols) {
						(*overscanLeftImage)[i][j] = pixelValue;
					} else if (i < Image::kNOverscanCols+Image::kNBlankCols) {
						(*blankLeftImage)[i-Image::kNOverscanCols][j] = pixelValue;
					} else if (i < Image::kLeftMargin) {
						(*darkLeftImage)[i-Image::kNOverscanCols-Image::kNBlankCols][j] = pixelValue;
					} else if (ix>=Image::kXDim && ix<Image::kXDim+Image::kNDarkCols) {
						(*darkRightImage)[ix-Image::kXDim][j] = pixelValue;
					} else if (ix>=Image::kXDim+Image::kNDarkCols) {
						(*blankRightImage)[ix-Image::kXDim-Image::kNDarkCols][j] = pixelValue;
					}

				} else if (ix>=0 && ix<Image::kXDim) {

					if (iy>=Image::kYDim && iy<Image::kYDim+Image::kNDarkRows) {
						(*darkTopImage)[ix][j-Image::kYDim] = pixelValue;
					} else if (iy>=Image::kYDim+Image::kNDarkRows) {
						(*overscanTopImage)[ix][j-Image::kYDim-Image::kNDarkRows] = pixelValue;
					}

				}

			}

		}
	}

	Averages blankLeftAverages = ccdMarginAverages(blankLeftImage);
	Averages blankRightAverages = ccdMarginAverages(blankRightImage);
	Averages darkLeftAverages = ccdMarginAverages(darkLeftImage);
	Averages darkRightAverages = ccdMarginAverages(darkRightImage);
	Averages darkTopAverages = ccdMarginAverages(darkTopImage);
	Averages overscanLeftAverages = ccdMarginAverages(overscanLeftImage);
	Averages overscanTopAverages = ccdMarginAverages(overscanTopImage);

	writeCcdMargin<SciRawBlankleft>(data,blankLeftImage,blankLeftAverages,0,Data::fullframe,ImageWriter::image);
	writeCcdMargin<SciRawBlankright>(data,blankRightImage,blankRightAverages,0,Data::fullframe,ImageWriter::image);
	writeCcdMargin<SciRawDarkleft>(data,darkLeftImage,darkLeftAverages,0,Data::fullframe,ImageWriter::image);
	writeCcdMargin<SciRawDarkright>(data,darkRightImage,darkRightAverages,0,Data::fullframe,ImageWriter::image);
	writeCcdMargin<SciRawDarktop>(data,darkTopImage,darkTopAverages,0,Data::fullframe,ImageWriter::image);
	if (m_redundantHardware) {
		//For the redundant hardware case, write the overscan data to SCI_RAW_OverscanRight (the overscan data is in the left-most
		//column of the Image array, hence overscanLeftImage,overscanLeftAverages, but we write it to SCI_RAW_OverscanRight)
		writeCcdMargin<SciRawOverscanright>(data,overscanLeftImage,overscanLeftAverages,0,Data::fullframe,ImageWriter::image);
	} else {
		writeCcdMargin<SciRawOverscanleft>(data,overscanLeftImage,overscanLeftAverages,0,Data::fullframe,ImageWriter::image);
	}
	writeCcdMargin<SciRawOverscantop>(data,overscanTopImage,overscanTopAverages,0,Data::fullframe,ImageWriter::image);

	delete overscanLeftImage;
	delete blankLeftImage;
	delete darkLeftImage;
	delete blankRightImage;
	delete darkRightImage;
	delete darkTopImage;
	delete overscanTopImage;

	//Set metadata header keywords
	SciRawImagemetadata * metaData = Append_SciRawImagemetadata(fullFrameImage);
	setVisitKeywords(metaData);
	setTargetKeywords(metaData);
	metaData->setKeyVStrtU(fullFrameImage->getKeyVStrtU());
	metaData->setKeyVStrtM(fullFrameImage->getKeyVStrtM());
	metaData->setKeyVStopU(fullFrameImage->getKeyVStopU());
	metaData->setKeyVStopM(fullFrameImage->getKeyVStopM());
	metaData->setKeyProcChn(fullFrameImage->getKeyProcChn()); //kept outside setTargetKeywords because not used for truth metadata
	metaData->setKeyRdMode("full frame");

	//Compression Entity keywords
	metaData->setKeyAcqMode(5); //1: DUMP (sub-array) 2: DIGIT 3: FULL
	metaData->setKeyOversamp(true); //ToDo: Check this is correct. Documentation says: if true then averaging of several exposures is done
	metaData->setKeyFSource(0); //0: CCD 1: PATTERN 2:SIMULATION
	metaData->setKeyRepetit(timeConf.getRepetitionPeriod());

	//Fill metadata table row

	metaData->setCellObtTime(fullFrameImage->getKeyTStrtO());
	metaData->setCellUtcTime(fullFrameImage->getKeyTStrtU());
	metaData->setCellMjdTime(fullFrameImage->getKeyTStrtM());

	//Get satellite data
	double telescopeTemperature = SatelliteData::kDefaultTelescopeTemperature;
	double dpuTemperature = SatelliteData::kDefaultDpuTemperature;
	if (data->hasSatelliteData()) {
		SatelliteData * satData = data->getSatelliteData(timeStep);
		dpuTemperature = satData->getDpuTemperature();
		telescopeTemperature = satData->getTelescopeTemperature();
		metaData->setCellLosToSunAngle(satData->getSunAngle());
		metaData->setCellLosToMoonAngle(satData->getMoonAngle());
		metaData->setCellLosToEarthAngle(satData->getEarthLimbAngle());
		metaData->setCellLatitude(satData->getLatitude());
		metaData->setCellLongitude(satData->getLongitude());
	} else {
		metaData->setNullLosToSunAngle();
		metaData->setNullLosToMoonAngle();
		metaData->setNullLosToEarthAngle();
		metaData->setNullLatitude();
		metaData->setNullLongitude();
	}

	//Get HK data
	Data::HKData hkData = data->getClosestAveragedHKData((midTime_utc-timeConf.getVisitStartTimeUTC()).getSeconds());

	int ceCounter = m_ceCounterOffset+1;
	if (!initialFullFrame) {
		//CE counter for sub-array images starts from 1
		//unless counter=1 is already taken by a full frame image for timeStep=-1 in which case the counter starts from 2
		int indexOffset = m_ceCounterOffset + (data->doFullFrame() ? 2 : 1);
		ceCounter = data->stackedImageCount()+data->getNumberOfDiscardedImages()+indexOffset;
	}
	metaData->setCellCeCounter(ceCounter);
	metaData->setCellCeIntegrity(0);
	metaData->setCellPixDataOffset(static_cast<uint16_t>(lround(data->getBiasOffset())+0.5));
	metaData->setCellCcdTimingScript(readoutScript(data->getVisit().m_readMode,data->getTimeConfiguration().getExposureTimeAsDouble(),Data::fullframe));
	metaData->setCellHkSource("hk tm");
	metaData->setCellHkTempFeeCcd(hkData.m_ccdTemp-273.15);
	metaData->setCellHkTempFeeBias(hkData.m_biasTemp-273.15);
	metaData->setCellHkTempFeeAdc(hkData.m_adcTemp-273.15);
	metaData->setCellAdcTemp1(dpuTemperature-273.15);
	metaData->setCellAdcN5v(SatelliteData::kDefaultDpuVoltage);
	metaData->setCellHkVoltFeeVod(hkData.m_vod);
	metaData->setCellHkVoltFeeVrd(hkData.m_vrd);
	metaData->setCellHkVoltFeeVog(hkData.m_vog);
	metaData->setCellHkVoltFeeVss(hkData.m_vss);
	metaData->setCellThermaft4(telescopeTemperature-273.15+0.4+0.5);
	metaData->setCellThermaft3(telescopeTemperature-273.15+0.3+0.5);
	metaData->setCellThermaft2(telescopeTemperature-273.15+0.2+0.5);
	metaData->setCellThermaft1(telescopeTemperature-273.15+0.1+0.5);
	metaData->setCellThermfront1(telescopeTemperature-273.15+0.5);
	metaData->setCellThermfront2(telescopeTemperature-273.15-0.1+0.5);
	metaData->setCellThermfront3(telescopeTemperature-273.15-0.2+0.5);
	metaData->setCellThermfront4(telescopeTemperature-273.15-0.3+0.5);
	metaData->setCellLeftDarkColMask(65535);
	metaData->setCellRightDarkColMask(65535);
	metaData->WriteRow();

	//Set unstacked metadata header keywords
	SciRawUnstackedimagemetadata * metaData_unstacked = Append_SciRawUnstackedimagemetadata(fullFrameImage);
	setVisitKeywords(metaData_unstacked);
	setTargetKeywords(metaData_unstacked);
	metaData_unstacked->setKeyVStrtU(fullFrameImage->getKeyVStrtU());
	metaData_unstacked->setKeyVStrtM(fullFrameImage->getKeyVStrtM());
	metaData_unstacked->setKeyVStopU(fullFrameImage->getKeyVStopU());
	metaData_unstacked->setKeyVStopM(fullFrameImage->getKeyVStopM());
	metaData_unstacked->setKeyProcChn(fullFrameImage->getKeyProcChn()); //kept outside setVisitKeywords because not used for truth metadata

	//Fill unstacked metadata table row
	metaData_unstacked->setCellObtTime(fullFrameImage->getKeyTStrtO());
	metaData_unstacked->setCellUtcTime(fullFrameImage->getKeyTStrtU());
	metaData_unstacked->setCellMjdTime(fullFrameImage->getKeyTStrtM());
	metaData_unstacked->setCellCeCounter(ceCounter);
	metaData_unstacked->setCellGain0(m_nominalGain);
	metaData_unstacked->setCellBias(overscanLeftAverages.mean);
	metaData_unstacked->setCellBias0(data->getBiasOffset());
	metaData_unstacked->setCellCeTempFeeCcd(hkData.m_ccdTemp-273.15);
	metaData_unstacked->setCellCeVoltFeeVod(hkData.m_vod);
	metaData_unstacked->setCellCeVoltFeeVrd(hkData.m_vrd);
	metaData_unstacked->setCellCeVoltFeeVog(hkData.m_vog);
	metaData_unstacked->setCellCeVoltFeeVss(hkData.m_vss);
	metaData_unstacked->WriteRow();

	//Write truth information
	if (m_writeTruthData) {

		SimTruFullarray * truthMetaData = createFitsFile<SimTruFullarray>(m_dirname+"/truth/", midTime_utc,
				                                                          VisitId(m_visit.m_progType,m_visit.m_progId,m_visit.m_reqId,m_visit.m_visitCtr), PassId());
		setVisitKeywords(truthMetaData);
		setTargetKeywords(truthMetaData);
		truthMetaData->setKeyVStrtU(fullFrameImage->getKeyVStrtU());
		truthMetaData->setKeyVStrtM(fullFrameImage->getKeyVStrtM());
		truthMetaData->setKeyVStopU(fullFrameImage->getKeyVStopU());
		truthMetaData->setKeyVStopM(fullFrameImage->getKeyVStopM());

		truthMetaData->setCellObtTime(fullFrameImage->getKeyTStrtO());
		truthMetaData->setCellUtcTime(fullFrameImage->getKeyTStrtU());
		truthMetaData->setCellMjdTime(fullFrameImage->getKeyTStrtM());

		vector<SatelliteData*> satData;
		if (data->hasSatelliteData()) {
			satData.push_back(data->getSatelliteData(timeStep));
			truthMetaData->setCellValidAocs(satData[0]->validAocs());
			truthMetaData->setCellValidScience(satData[0]->validScience());
		}
		truthMetaData->setCellFullWellSaturated(image->getTruthData()->isFullWellSaturated());
		truthMetaData->setCellAdcSaturated(image->getTruthData()->isAdcSaturated());
		truthMetaData->setCellGlobalThroughput(image->getTruthData()->getGlobalThroughput());
		truthMetaData->setCellGain(image->getTruthData()->getGain());
		truthMetaData->setCellStrayLight(image->getTruthData()->getStrayLight());
		truthMetaData->setCellZodiacalLight(image->getTruthData()->getZodiacalLight());

		fillTruthData<SimTruFullarray>(truthMetaData,image->getTruthData(),satData,0,0,Image::kXDim,Image::kYDim,true,false);

		delete truthMetaData;

	}

	//Write the centroid information
	if (m_writeCentroid) {

		//Get the pointer to the centroid file
		if (data->getFitsCentroid() == nullptr) initializeCentroid(data);
		SciRawCentroid * centroid = data->getFitsCentroid();

		//Determine PSF barycentre using truth information
		//Calculation without using truth information not currently implemented for full frame image
		Data::SubarrayDimensions subarray = data->getSubarrayDimensions();
		pair<double,double> barycentre = make_pair(0.,0.);
		if (m_truthBarycentre) {
			vector<PSF> truthPSFs = image->getTruthData()->getPSFs();
			if (truthPSFs.size()>0) {
				barycentre = make_pair(truthPSFs[0].getMeanXPosition() - subarray.m_targetLocationX,
									   truthPSFs[0].getMeanYPosition() - subarray.m_targetLocationY);
			}
		}

		//Fill centroid table
		centroid->setCellUtcStart(startTime_utc);
		centroid->setCellMjdStart(startTime_utc.getMjd());
		centroid->setCellObtStart(startTime_utc.getObt());
		centroid->setCellUtcStop(endTime_utc);
		centroid->setCellMjdStop(endTime_utc.getMjd());
		centroid->setCellObtStop(endTime_utc.getObt());
		centroid->setCellFullFrame(true);
		centroid->setCellCeCounter(ceCounter);
		if (m_targetLocationFromJitter) {
			centroid->setCellOffsetX(0.);
			centroid->setCellOffsetY(0.);
			centroid->setCellLocationX(subarray.m_targetLocationX + barycentre.first);
			centroid->setCellLocationY(subarray.m_targetLocationY + barycentre.second);
		} else {
			centroid->setCellOffsetX(barycentre.first);
			centroid->setCellOffsetY(barycentre.second);
			centroid->setCellLocationX(subarray.m_targetLocationX);
			centroid->setCellLocationY(subarray.m_targetLocationY);
		}
		centroid->setCellDataCadence(timeConf.getRepetitionPeriod());
		centroid->setCellValidity(1);

		centroid->WriteRow();

	}

	delete metaData;
	delete metaData_unstacked;
	delete fullFrameImage;

}

double ImageWriter::getImageDuration(TimeConfiguration timeConf, Data::IMAGETYPE type) const {

	double repetitionPeriod = timeConf.getRepetitionPeriod();
	double exposuresPerStack = timeConf.getExposuresPerStack();
	double exposureTime = timeConf.getExposureTimeAsDouble();
	double imagetteStackingNumber = timeConf.getImagetteStackingNumber();

	//For stacked images/imagettes, the duration is from the start of the first exposure to the end of the last exposure
	//It does not include the gap after the end of the last exposure in the case where the repetition period is longer than the exposure duration
	double duration;
	if (type==Data::stacked) {
		duration = (exposuresPerStack-1)*repetitionPeriod + exposureTime;
	} else if (type==Data::imagette && imagetteStackingNumber>1) {
		duration = (imagetteStackingNumber-1)*repetitionPeriod + exposureTime;
	} else {
		duration = exposureTime;
	}

	return duration;

}
