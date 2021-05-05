/*
 * DarkCurrentGenerator.cxx
 *
 *  Created on: 11 Mar 2014
 *      Author: futyand
 */

#include <fstream>

#include "boost/algorithm/string.hpp"
#include "boost/filesystem.hpp"

#include "REF_APP_DarkFrame.hxx"
#include "REF_APP_GainCorrection.hxx"

#include "DarkCurrentGenerator.hxx"

DarkCurrentGenerator::DarkCurrentGenerator() : Module("DarkCurrentGenerator",timeLoop) {

	for (int ix=0; ix<Image::kXDim+2*Image::kNDarkCols; ix++) {
		for (int iy=0; iy<Image::kYDim+Image::kNDarkRows; iy++) {
			m_isHotPixel[ix][iy] = false;
			m_isWarmPixel[ix][iy] = false;
		}
	}

}

DarkCurrentGenerator::~DarkCurrentGenerator() {
	delete m_exponentialRandomGenerator;
	delete m_uniformRandomGeneratorX;
	delete m_uniformRandomGeneratorY;
}

void DarkCurrentGenerator::initialize(const ModuleParams& params) {

	m_empiricalDarkFrame = params.GetAsBool("empiricalDarkFrame");
	m_darkFrameFilename = params.GetAsString("darkFrameFilename");
	m_darkFrameScaleFactor = params.GetAsDouble("darkFrameScaleFactor");

	m_meanDarkCurrent233K = params.GetAsDouble("meanDarkCurrent233K");
	if (params.GetAsBool("includeReadoutTime")) {
		m_rowDownShiftTime = params.GetAsDouble("rowDownShiftTime")/1.E6;
		m_serialShiftTime = 1./(params.GetAsDouble("serialReadRate")*1000.);
	} else {
		m_rowDownShiftTime = 0.;
		m_serialShiftTime = 0.;
	}

	unsigned hotPixelPositionSeed = params.GetAsInt("hotPixelPositionSeed");
	unsigned telegraphicTransitionSeed = params.GetAsInt("telegraphicTransitionSeed");
	m_telegraphicTimeConstant = params.GetAsDouble("telegraphicTimeConstant");

	m_doHotPixels = params.GetAsBool("doHotPixels");
	double fracHotPixels = params.GetAsDouble("fracHotPixels");
	m_hotPixelRelativeDarkCurrent = params.GetAsDouble("hotPixelRelativeDarkCurrent");

	m_doWarmPixels = params.GetAsBool("doWarmPixels");
	double fracWarmPixels = params.GetAsDouble("fracWarmPixels");
	m_warmPixelRelativeDarkCurrent = params.GetAsDouble("warmPixelRelativeDarkCurrent");

	m_doTelegraphicPixels = params.GetAsBool("doTelegraphicPixels");
	double fracTelegraphicPixels = params.GetAsDouble("fracTelegraphicPixels");
	m_telegraphicPixelRelativeDarkCurrent = params.GetAsDouble("telegraphicPixelRelativeDarkCurrent");

	m_doManualHotPixels = params.GetAsBool("doManualHotPixels");

	//Convert the fraction to the number of hot and telegraphic pixels corresponding to full frame plus dark reference columns/rows
	int nPixels = (Image::kXDim+2*Image::kNDarkCols) * (Image::kYDim+Image::kNDarkRows);
	m_nHotPixels = static_cast<int>(lround(fracHotPixels * nPixels));
	m_nWarmPixels = static_cast<int>(lround(fracWarmPixels * nPixels));
	m_nTelegraphicPixels = static_cast<int>(lround(fracTelegraphicPixels * nPixels));

	//Initialize the random number generators
	boost::mt19937 randomNumberEngine2(hotPixelPositionSeed);
	boost::mt19937 randomNumberEngine3(hotPixelPositionSeed+1);
    boost::uniform_int<int> uniform_randomX(0,Image::kXDim+2*Image::kNDarkCols-1);
    boost::uniform_int<int> uniform_randomY(0,Image::kYDim+Image::kNDarkRows-1);
    m_uniformRandomGeneratorX = new RANDOM_UNIFORM(randomNumberEngine2,uniform_randomX);
    m_uniformRandomGeneratorY = new RANDOM_UNIFORM(randomNumberEngine3,uniform_randomY);

	boost::mt19937 randomNumberEngine4(telegraphicTransitionSeed);
    boost::exponential_distribution<float> exponentialDistribution(1./m_telegraphicTimeConstant);
    m_exponentialRandomGenerator = new RANDOM_EXPONENTIAL(randomNumberEngine4,exponentialDistribution);
    boost::bernoulli_distribution<> bernoulliDistribution(0.5);
    m_bernoulliRandomGenerator = new RANDOM_BERNOULLI(randomNumberEngine4,bernoulliDistribution);

	//Read in the empirical dark frame
    if (m_empiricalDarkFrame) readDarkFrame();

	//Get the manually generated hot/warm/telegraphic pixels
	if (m_doManualHotPixels) {

		string manualHotPixelX_str = params.GetAsString("manualHotPixelX");
		string manualHotPixelY_str = params.GetAsString("manualHotPixelY");
		string manualHotPixelRate_str = params.GetAsString("manualHotPixelRate");
		string manualHotPixelIsTelegraphic_str = params.GetAsString("manualHotPixelIsTelegraphic");

		vector<string> substrings_x,substrings_y,substrings_rate,substrings_isTelegraphic;
		boost::split(substrings_x,manualHotPixelX_str,boost::is_any_of(","));
		boost::split(substrings_y,manualHotPixelY_str,boost::is_any_of(","));
		boost::split(substrings_rate,manualHotPixelRate_str,boost::is_any_of(","));
		boost::split(substrings_isTelegraphic,manualHotPixelIsTelegraphic_str,boost::is_any_of(","));

		if (substrings_x.size() != substrings_y.size() ||
			substrings_x.size() != substrings_rate.size() ||
			substrings_x.size() != substrings_isTelegraphic.size()) {
			throw runtime_error("Error in DarkCurrentGenerator::initialize: inconsistent array sizes for position, rate and telegraphic flag for manually defined hot pixels");
		}

		for (unsigned i=0; i<substrings_x.size(); i++) {
			int ix = boost::lexical_cast<int>(substrings_x[i]);
			int iy = boost::lexical_cast<int>(substrings_y[i]);
			double rate = boost::lexical_cast<double>(substrings_rate[i]);
			double isTelegraphic = substrings_isTelegraphic[i]=="T";
			m_manualHotPixels.push_back(ManualHotPixel(ix,iy,rate,isTelegraphic));
		}

	}

}

void DarkCurrentGenerator::readDarkFrame() {

	//Open the dark frame fits file
	cout << "Dark frame file: " << (string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_darkFrameFilename) << endl;
	RefAppDarkframe * darkFrame = new RefAppDarkframe(string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_darkFrameFilename,"READONLY");
	RefAppDarkframeleft * darkLeft = Open_RefAppDarkframeleft(darkFrame);
	RefAppDarkframeright * darkRight = Open_RefAppDarkframeright(darkFrame);
	RefAppDarkframetop * darkTop = Open_RefAppDarkframetop(darkFrame);

	ofstream referenceFilesList("reference_files.txt", ios::app);
	referenceFilesList << m_darkFrameFilename << endl;
	referenceFilesList.close();

	//Fill the dark frame array with pixel values from the empirical frame and calculate the mean value
	double meanDarkCurrent = 0.;
	unsigned nPixels = 0;
	//unsigned nWarm=0;
	//unsigned nHot=0;
	for (int i=0; i<Image::kXDim+2*Image::kNDarkCols; i++) {
		for (int j=0; j<Image::kYDim+Image::kNDarkRows; j++) {
			double pixel_value = 0.;
			if (i<Image::kNDarkCols && j<Image::kYDim) {
				pixel_value = (*darkLeft)[0][j][i];
			} else if (i>=Image::kXDim+Image::kNDarkCols && j<Image::kYDim) {
				pixel_value = (*darkRight)[0][j][i-Image::kXDim-Image::kNDarkCols];
			} else if (i>=Image::kNDarkCols && i<Image::kXDim+Image::kNDarkCols && j>=Image::kYDim) {
				pixel_value = (*darkTop)[0][j-Image::kYDim][i-Image::kNDarkCols];
			} else  if (i>=Image::kNDarkCols && i<Image::kXDim+Image::kNDarkCols && j<Image::kYDim) {
				pixel_value = (*darkFrame)[0][j][i-Image::kNDarkCols];
			}
			if (pixel_value < 0.) pixel_value=0.;
			pixel_value *= m_darkFrameScaleFactor;
			meanDarkCurrent += pixel_value;
			nPixels += 1;
			m_darkFrame[Image::kNOverscanCols+Image::kNBlankCols+i][j] = pixel_value;
//			if (pixel_value>m_hotPixelRelativeDarkCurrent*m_meanDarkCurrent233K) {
//				nHot++;
//				cout << i-Image::kNDarkCols-2 << " " << j-2 << " " << pixel_value << endl;
//			} else if (pixel_value>m_warmPixelRelativeDarkCurrent*m_meanDarkCurrent233K) {
//				if (i-Image::kNDarkCols>=412 && i-Image::kNDarkCols<612 &&j>=412 && j<612) {
//					nWarm++;
//					cout << i-Image::Image::kNDarkCols-2 << " " << j-2 << " " << pixel_value << endl;
//				}
//			}
		}
	}
	meanDarkCurrent /= nPixels;
	//cout << nWarm << " " << nHot << " " << m_hotPixelRelativeDarkCurrent*m_meanDarkCurrent233K << " " << m_warmPixelRelativeDarkCurrent*m_meanDarkCurrent233K << endl;

	//Assign the dark frame values for blank and overscan columns to the mean of the values for all physical pixels
	for (int ix=0; ix<Image::kXTotal; ix++) {
		for (int iy=0; iy<Image::kYTotal; iy++) {
			if (ix<Image::kNOverscanCols+Image::kNBlankCols || ix>=Image::kXTotal-Image::kNBlankCols || iy>=Image::kYDim+Image::kNDarkRows) {
				m_darkFrame[ix][iy] = meanDarkCurrent;
			}
		}
	}

	delete darkFrame;
	delete darkLeft;
	delete darkRight;
	delete darkTop;

}

void DarkCurrentGenerator::doBegin(Data *data, bool fullFrame) {

	if (m_serialShiftTime>0.) data->setSerialReadRate((1./m_serialShiftTime)/1000.); //Used to set the RO_FREQU keyword in SCI_RAW_SubArray and SCI_RAW_FullArray in ImageWriter

	//Use faint readout mode (read out all pixels) for exposure times >11.296s.
	//For shorter exposures, use faint fast readout mode (read out only the 200 rows corresponding to the sub-array)
	double exposureTime = data->getTimeConfiguration().getExposureTimeAsDouble();
	bool faintFastMode = (exposureTime < 11.296);

	int yDim = data->getSubarrayDimensions().m_yDim;
	int yOffset = data->getSubarrayDimensions().m_yOffset;

	//Determine the total exposure time for each pixel, taking into account the time between the end of the frame transfer and the pixel being read out
	//For blank and overscan columns, exposure time corresponds to the read time for the current row only because these are not physical pixels,
	//except for the blanks in the first row (y==0), for which the full exposure time applies because these are physical pixels.
	//Time for frame transfer is assumed to be negligible
	double readExposureTime = 0.;
	for (int iy=0; iy<Image::kYTotal; iy++) {
		double rowExposureTime = 0.;
		for (int ix=Image::kXTotal-1; ix>=0; ix--) {

			//Accumulate exposure time during readout
			//In bright readout mode, only the 200 rows corresponding to the sub-array are read out
			if (!faintFastMode || (iy>=yOffset && iy<yOffset+yDim)) {
				readExposureTime += m_serialShiftTime;
				//Accumulate exposure time during the current row to determine read delay for blank and overscan pixels in the row
				if (ix>=Image::kNOverscanCols) rowExposureTime += m_serialShiftTime; //omit non-physical overscan pixels
			}

			if ((ix>=Image::kNOverscanCols+Image::kNBlankCols && ix<Image::kXTotal-Image::kNBlankCols && iy<Image::kYDim+Image::kNDarkRows) ||
					(iy==0 && ix>=Image::kNOverscanCols)) {
				//physical pixels include exposure time plus accumulated read time
				m_exposureTimeWithReadout[ix][iy] = readExposureTime + static_cast<double>(exposureTime);
			} else if (iy>=Image::kYTotal-Image::kNOverscanRows) {
				//Top overscan rows include accumulated read time for the whole CCD but not the exposure time
				m_exposureTimeWithReadout[ix][iy] = readExposureTime;
			} else {
				//Overscan and blank columns include read time for the current row only
				m_exposureTimeWithReadout[ix][iy] = rowExposureTime;
			}

		}

		readExposureTime += m_rowDownShiftTime;

	}

	data->setReadoutTime(readExposureTime);
	cout << "Readout in " << ( (faintFastMode && exposureTime>2.326) ? "faint fast" : data->getVisit().m_readMode) << " mode. Total time to read out the CCD: " << " " << readExposureTime << " seconds" << endl;

	data->initializeBadPixelMap(m_darkFrameFilename);

	//Write the list of hot pixels from the empirical dark frame to the bad pixels reference file.
	//Pixel value is 1 for hot (dark current>5e-/s), otherwise 0.
	for (int i=0; i<Image::kXTotal; i++) {
		for (int j=0; j<Image::kYTotal; j++) {
			if (m_darkFrame[i][j]>5.) {
				data->setBadPixel(i,j,1);
				cout << "empirical hot pixel at (" << i-Image::kLeftMargin << "," << j << "): " << m_darkFrame[i][j] << " electrons/s" << endl;
			}
		}
	}

    //Randomly generate the positions of the hot pixels
    if (m_doHotPixels) {
    	for (int ihp=0; ihp<m_nHotPixels; ihp++) {
    		int ix_hp = (*m_uniformRandomGeneratorX)();
    		int iy_hp = (*m_uniformRandomGeneratorY)();
    		m_isHotPixel[ix_hp][iy_hp] = true;
    		//cout << "hot pixel: " << ix_hp << " " << iy_hp << " " << m_isHotPixel[ix_hp][iy_hp] << endl;

			//Add the hot pixel to the bad pixel reference file
    		data->setBadPixel(Image::kNOverscanCols+Image::kNBlankCols+ix_hp,iy_hp,1);
    	}
    }

    //Randomly generate the positions of the warm pixels
    if (m_doWarmPixels) {
    	for (int ihp=0; ihp<m_nWarmPixels; ihp++) {
    		int ix_hp = (*m_uniformRandomGeneratorX)();
    		int iy_hp = (*m_uniformRandomGeneratorY)();
    		m_isWarmPixel[ix_hp][iy_hp] = true;
    		//cout << "warm pixel: " << ix_hp << " " << iy_hp << " " << m_isWarmPixel[ix_hp][iy_hp] << endl;
    	}
    }

    //Add the manually defined hot pixels to the bad pixel reference file
    if (m_doManualHotPixels) {
    	RefAppGaincorrection * gainCorrection_file = new RefAppGaincorrection(string(getenv("CHEOPS_SW"))+"/resources/"+data->getGainFilename(),"READONLY");
    	for (unsigned ihp=0; ihp<m_manualHotPixels.size(); ihp++) {
    		cout << "manually defined hot pixel: " << m_manualHotPixels[ihp].m_xPixel << " " << m_manualHotPixels[ihp].m_yPixel << " " << m_manualHotPixels[ihp].m_rate << " " << (m_manualHotPixels[ihp].m_isTelegraphic ? "telegraphic":"") << endl;
    		int type = 1;
    		if (m_manualHotPixels[ihp].m_rate*exposureTime > pow(2,16)/gainCorrection_file->getKeyGainNom()) type = 2;
    		if (m_manualHotPixels[ihp].m_isTelegraphic) type = 3;
    		data->setBadPixel(Image::kLeftMargin+data->getSubarrayDimensions().m_xOffset+m_manualHotPixels[ihp].m_xPixel,
    						  data->getSubarrayDimensions().m_yOffset+m_manualHotPixels[ihp].m_yPixel,type);
    	}
    	delete gainCorrection_file;
    }

	//Initialize the telegraphic pixels
    if (m_doTelegraphicPixels || m_doManualHotPixels) setTelegraphicPixels(data);

}

void DarkCurrentGenerator::setTelegraphicPixels(Data* data) {

	if (m_doTelegraphicPixels) {
		//Randomly generate the positions of the telegraphic pixels
		//unsigned n_tp_window=0;
		for (int itp=0; itp<m_nTelegraphicPixels; itp++) {
			int ix_tp = (*m_uniformRandomGeneratorX)();
			int iy_tp = (*m_uniformRandomGeneratorY)();
			//if (ix_tp-Image::kNDarkCols>=412 && ix_tp-Image::kNDarkCols<612 &&iy_tp>=412 && iy_tp<612) n_tp_window++;
			double meanActiveRate = m_telegraphicPixelRelativeDarkCurrent*m_meanDarkCurrent233K;
			bool active = (*m_bernoulliRandomGenerator)(); //Randomly assign initial state
			vector<float> transitions = telegraphicTransitions(data->getTimeConfiguration().getDuration());
			Data::TelegraphicPixel telegraphicPixel(ix_tp-Image::kNDarkCols,iy_tp,meanActiveRate,active,transitions);
			data->addTelegraphicPixel(telegraphicPixel);
		}
		//cout << "N telegraphic in 200x200 = " << n_tp_window << endl;
	}

	if (m_doManualHotPixels) {
		//Add the manually defined telegraphic pixels
		for (unsigned itp=0; itp<m_manualHotPixels.size(); itp++) {
			if (m_manualHotPixels[itp].m_isTelegraphic) {
				int ix_tp = m_manualHotPixels[itp].m_xPixel;
				int iy_tp = m_manualHotPixels[itp].m_yPixel;
				double meanActiveRate = m_manualHotPixels[itp].m_rate;
				bool active = (*m_bernoulliRandomGenerator)(); //Randomly assign initial state
				vector<float> transitions = telegraphicTransitions(data->getTimeConfiguration().getDuration());
				Data::TelegraphicPixel telegraphicPixel(ix_tp,iy_tp,meanActiveRate,active,transitions);
				data->addTelegraphicPixel(telegraphicPixel);
			}
		}
	}

	//Add the telegraphic pixels to the bad pixel reference file
	vector<Data::TelegraphicPixel> * telegraphicPixels = data->getTelegraphicPixels();
	for (unsigned i=0; i<(*telegraphicPixels).size(); i++) {
		cout << "Adding telegraphic pixel at position (" << (*telegraphicPixels)[i].m_xPixel << "," << (*telegraphicPixels)[i].m_yPixel
														 << "), mean active rate " << (*telegraphicPixels)[i].m_meanActiveRate
														 << " electrons/s, initially " << ((*telegraphicPixels)[i].m_active ? "active" : "inactive") << endl;
		data->setBadPixel(Image::kLeftMargin+data->getSubarrayDimensions().m_xOffset+(*telegraphicPixels)[i].m_xPixel,
						  data->getSubarrayDimensions().m_yOffset+(*telegraphicPixels)[i].m_yPixel,3);
	}

}

vector<float> DarkCurrentGenerator::telegraphicTransitions(double duration) const {

	//Assign the list of times for which a transition occurs (active->inactive or inactive->active),
	//by drawing randomly from an exponential distribution with user defined time constant
	vector<float> transitions;
	float itime = 0;
	while (itime < duration) {
		float timeToNextTransition = (*m_exponentialRandomGenerator)();
		itime += timeToNextTransition;
		transitions.push_back(itime);
		//cout << itime << endl;
	}

	return transitions;

}

void DarkCurrentGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	//Get instances of the primary image, together with images extended
	//to cover the full frame in vertical direction,
	//corresponding to the first and last seconds of the exposure,
	//needed to generate the upward and downward smear trails in FrameTransferSmearer
	Image* image = data->getImages().back();
	Image * imageToSmearUp = nullptr;
	Image * imageToSmearDown = nullptr;
	if (data->doFrameTransferSmearing()) {
		if (data->getImagesToSmearUp().size() == data->getImages().size()) { //Smear trail images already created in PSFGenerator
			imageToSmearUp = data->getImagesToSmearUp().back();
			imageToSmearDown = data->getImagesToSmearDown().back();
		} else { //Smear trail images for do not yet exist for the current time step because PSFGenerator was not run
			imageToSmearUp = new Image(image->getXDim(),Image::kYDim,image->getXOffset(),0);
			imageToSmearDown = new Image(image->getXDim(),Image::kYDim,image->getXOffset(),0);
			data->addImageToSmearUp(imageToSmearUp);
			data->addImageToSmearDown(imageToSmearDown);
		}
	}

	//Get the current temperature of the CCD
	double temperature = data->getSatelliteData(timeStep)->getCcdTemperature();

	//Calculate the factor by which the temperature affects the dark current,
	//using the relation given in the e2v CCD-4720 specification document.
	double temperatureFactor = darkTemperatureDependence(temperature)/darkTemperatureDependence(SatelliteData::kDefaultCcdTemperature);

	//Get the list of telegraphic pixels
	vector<Data::TelegraphicPixel> * telegraphicPixels;
	if (m_doTelegraphicPixels || m_doManualHotPixels) telegraphicPixels = data->getTelegraphicPixels();
	double repetitionPeriod = data->getTimeConfiguration().getRepetitionPeriod();

	bool ctiEoL = data->doChargeTransferEoL();

	int xmin,xmax,ymin,ymax;
	if (data->doFrameTransferSmearing() && imageToSmearUp != nullptr){
		xmin = imageToSmearUp->getXOffset();
		xmax = imageToSmearUp->getXOffset()+imageToSmearUp->getXDim();
		ymin = imageToSmearUp->getYOffset();
		ymax = imageToSmearUp->getYOffset()+imageToSmearUp->getYDim();
	} else {
		xmin = image->getXOffset();
		xmax = image->getXOffset()+image->getXDim();
		ymin = ctiEoL ? 0. : image->getYOffset(); //Include pixels below the subarray if end of life CTI is to be simulated
		ymax = image->getYOffset()+image->getYDim();
	}

	//Generate dark electrons for full frame or sub-array image
	for (int ix=xmin; ix<xmax; ix++) {
		for (int iy=ymin; iy<ymax; iy++) {
			generateDarkElectrons(image,imageToSmearUp,imageToSmearDown,ix,iy,
								  temperatureFactor,telegraphicPixels,timeStep,repetitionPeriod,ctiEoL);
		}
		//Generate noise also for dark and overscan rows
		for (int iy=Image::kYDim; iy<Image::kYTotal; iy++) {
			generateDarkElectrons(image,imageToSmearUp,imageToSmearDown,ix,iy,
								  temperatureFactor,telegraphicPixels,timeStep,repetitionPeriod,ctiEoL);
		}
	}
	//Generate noise also for dark, blank and overscan columns
	for (int iy=ymin; iy<ymax; iy++) {
		for (int ix=-Image::kLeftMargin; ix<0; ix++) {
			generateDarkElectrons(image,imageToSmearUp,imageToSmearDown,ix,iy,
								  temperatureFactor,telegraphicPixels,timeStep,repetitionPeriod,ctiEoL);
		}
		for (int ix=Image::kXDim; ix<Image::kXTotal-Image::kLeftMargin; ix++) {
			generateDarkElectrons(image,imageToSmearUp,imageToSmearDown,ix,iy,
								  temperatureFactor,telegraphicPixels,timeStep,repetitionPeriod,ctiEoL);
		}
	}

}

void DarkCurrentGenerator::generateDarkElectrons(Image* image, Image * imageToSmearUp, Image * imageToSmearDown,
		int ix, int iy, double temperatureFactor,
		vector<Data::TelegraphicPixel> * telegraphicPixels, int timeStep, double repetitionPeriod, bool doChargeTransferEoL) const {

	//Define the mean dark current per second for the pixel at 233K either using the user input value (same for all pixels),
	//or using the empirical dark frame
	double meanDarkCurrentRate233K = m_empiricalDarkFrame ? m_darkFrame[ix+Image::kLeftMargin][iy] : m_meanDarkCurrent233K;

	//Calculate dark current per second depending on the temperature using the relation given in the e2v CCD-4720 specification document.
	double meanDarkCurrentRate = meanDarkCurrentRate233K * temperatureFactor;

	//Multiply the dark current per second by the exposure duration to get the expected number of dark electrons in the pixel
	double meanDarkCurrent = m_exposureTimeWithReadout[ix+Image::kLeftMargin][iy] * meanDarkCurrentRate;

	bool isManualHot = false;
	double manualHotRate = 0.;
	if (m_doManualHotPixels) {
		//Determine whether or not the current pixel is a manually defined hot pixel
		for (unsigned i=0; i<m_manualHotPixels.size(); i++) {
			if (image->getXOffset()+m_manualHotPixels[i].m_xPixel == ix &&
				image->getYOffset()+m_manualHotPixels[i].m_yPixel == iy &&
				!m_manualHotPixels[i].m_isTelegraphic) {
				isManualHot = true;
				manualHotRate = m_manualHotPixels[i].m_rate;
			}
		}
	}

	bool isTelegraphic = false;
	double telegraphicFractionActive = 0.;
	double telegraphicActiveRate = meanDarkCurrent;
	if (m_doTelegraphicPixels || m_doManualHotPixels) {
		//Determine whether or not the current pixel is an active telegraphic pixel
		for (unsigned i=0; i<(*telegraphicPixels).size(); i++) {
			if (image->getXOffset()+(*telegraphicPixels)[i].m_xPixel == ix &&
				image->getYOffset()+(*telegraphicPixels)[i].m_yPixel == iy) {

				isTelegraphic = true;
				telegraphicActiveRate = (*telegraphicPixels)[i].m_meanActiveRate * m_exposureTimeWithReadout[ix+Image::kLeftMargin][iy];

				//Get the list of transitions for this pixel within the exposure, as a fraction to the exposure duration
				vector<double> transitionsWithinExposure;
				for (unsigned itrans=0; itrans<(*telegraphicPixels)[i].m_transitions.size(); itrans++) {
					if ((*telegraphicPixels)[i].m_transitions[itrans] >= timeStep*repetitionPeriod &&
						(*telegraphicPixels)[i].m_transitions[itrans] < (timeStep+1)*repetitionPeriod) {
						double exposureFraction = ((*telegraphicPixels)[i].m_transitions[itrans] - timeStep*repetitionPeriod)/repetitionPeriod;
						transitionsWithinExposure.push_back(exposureFraction);
					} else if ((*telegraphicPixels)[i].m_transitions[itrans] >= (timeStep+1)*repetitionPeriod) {
						break;
					}
				}

				//For each transition increment the fraction of the exposure for which the pixel is active
				//If there are no transitions within the exposure, the fraction is 0 or 1 depending on the current state
				for (unsigned itrans=0; itrans<transitionsWithinExposure.size()+1; itrans++) {
					if ((*telegraphicPixels)[i].m_active) {
						double activeStart = (itrans==0 ? 0. : transitionsWithinExposure[itrans-1]);
						double activeEnd = (itrans==transitionsWithinExposure.size() ? 1. : transitionsWithinExposure[itrans]);
						telegraphicFractionActive += (activeEnd - activeStart);
						//if (transitionsWithinExposure.size()>0) cout << itrans << " " << activeStart << " " << activeEnd << " " << telegraphicFractionActive << endl;
					}
					if (itrans<transitionsWithinExposure.size()) {
						(*telegraphicPixels)[i].flip(); //flip the state
						//cout << "transition for telegraphic pixel " << i << endl;
					}
				}
				//cout << "telegraphicFractionActive for telegraphic pixel " << i << ": " << telegraphicFractionActive << endl;

			}
		}
	}

	double nDarkElectrons;
	//generate hot and telegraphic pixels, only for physical pixels
	if ((ix>=-Image::kNDarkCols && ix<Image::kXDim+Image::kNDarkCols && iy<Image::kYDim+Image::kNDarkRows) &&
			((m_doHotPixels && m_isHotPixel[ix+Image::kNDarkCols][iy]) ||
			 (m_doWarmPixels && m_isWarmPixel[ix+Image::kNDarkCols][iy]) || isTelegraphic || isManualHot ||
			 meanDarkCurrentRate233K > 2.)) {

		TruthData::HOTPIXEL_TYPE type;
		if (isTelegraphic) {
			type = telegraphicFractionActive>0. ? TruthData::telegraphicActive : TruthData::telegraphicInactive;
			nDarkElectrons = meanDarkCurrent + (telegraphicActiveRate-meanDarkCurrent)*telegraphicFractionActive;
		} else if (isManualHot) {
			type = TruthData::hot;
			nDarkElectrons = manualHotRate * m_exposureTimeWithReadout[ix+Image::kLeftMargin][iy];
		} else if (m_doHotPixels && m_isHotPixel[ix+Image::kNDarkCols][iy]) {
			type = TruthData::hot;
			nDarkElectrons = meanDarkCurrent * m_hotPixelRelativeDarkCurrent;
		} else if (m_doWarmPixels && m_isWarmPixel[ix+Image::kNDarkCols][iy]) {
			type = TruthData::warm;
			nDarkElectrons = meanDarkCurrent * m_warmPixelRelativeDarkCurrent;
		} else if (meanDarkCurrentRate233K > 2.) {
			type = meanDarkCurrentRate233K > 5. ? TruthData::hot : TruthData::warm;
			nDarkElectrons = meanDarkCurrent;
		}
		//cout << "hot pixel: " << ix << " " << iy << " " << nDarkElectrons << " type = " << type << endl;

		//Write to truth data if the pixel is within the 200x200 subarray and associated CCD margins
		if ((ix>=image->getXOffset() && ix<image->getXOffset()+image->getXDim() &&
				((iy>=image->getYOffset() && iy<image->getYOffset()+image->getYDim()) ||
						(iy>=Image::kYDim && iy<Image::kYDim+Image::kNDarkRows))) ||
				((iy>=image->getYOffset() && iy<image->getYOffset()+image->getYDim() &&
						((ix>=-Image::kNDarkCols && ix<0) ||
								(ix>=Image::kXDim && ix<Image::kXDim+Image::kNDarkCols))))) {

			image->getTruthData()->addHotPixel(ix,iy,nDarkElectrons,type);

		}

	} else {
		nDarkElectrons = meanDarkCurrent;
	}

	//Increment the pixel value of the primary image if the pixel is within the 200x200 subarray and associated CCD margins
	int ymin = doChargeTransferEoL ? 0. : image->getYOffset(); //Include pixels below the subarray if end of life CTI is to be simulated
	if ((ix>=image->getXOffset() && ix<image->getXOffset()+image->getXDim() &&
			((iy>=ymin && iy<image->getYOffset()+image->getYDim()) ||
					(iy>=Image::kYDim && iy<Image::kYDim+Image::kNDarkRows))) ||
			((iy>=ymin && iy<image->getYOffset()+image->getYDim() &&
					((ix>=-Image::kNDarkCols && ix<0) ||
							(ix>=Image::kXDim && ix<Image::kXDim+Image::kNDarkCols))))) {
		image->incrementPixelValue(ix,iy,nDarkElectrons);
	}

	//For the vertically images extended to be used for smearing, set the number of dark electrons corresponding to a 1 second exposure
	if (imageToSmearUp != nullptr && imageToSmearDown != nullptr) {
		imageToSmearUp->incrementPixelValue(ix,iy,nDarkElectrons/m_exposureTimeWithReadout[ix+Image::kLeftMargin][iy]);
		imageToSmearDown->incrementPixelValue(ix,iy,nDarkElectrons/m_exposureTimeWithReadout[ix+Image::kLeftMargin][iy]);
	}

	//if (ix==512 && iy==512) cout << " Dark current: " << meanDarkCurrentRate233K << " " << meanDarkCurrentRate << " " << meanDarkCurrent << " " << nDarkElectrons << endl;

}

double DarkCurrentGenerator::darkTemperatureDependence(double temperature) const {
	return 1.14E6 * pow(temperature,3) * exp(-9080/temperature);
}
