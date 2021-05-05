/*
 * BiasGenerator.cxx
 *
 *  Created on: 11 Mar 2014
 *      Author: futyand
 */

#include <fstream>
#include <regex>

#include "boost/math/tools/stats.hpp"

#include "BiasGenerator.hxx"

#include "REF_APP_BiasFrame.hxx"
#include "REF_APP_CCDLinearisation230.hxx"
#include "REF_APP_CCDLinearisation100.hxx"
#include "REF_APP_GainCorrection.hxx"

void BiasGenerator::initialize(const ModuleParams& params) {

	m_empiricalBiasFrame = params.GetAsBool("empiricalBiasFrame");
	m_biasFrameFilename = params.GetAsString("biasFrameFilename");
	m_biasMean = params.GetAsDouble("biasMean");
	m_biasWidth = params.GetAsDouble("biasWidth");
	m_applyAnalogElectronicsStability = params.GetAsBool("applyAnalogElectronicsStability");
	m_doCcdNonLinearity = params.GetAsBool("doCcdNonLinearity");
	m_ccdLinearityFilename = params.GetAsString("ccdNonLinearityFilename");

	//Set up the Gaussian random number generator for the bias frame (bias offset and readout noise)
	unsigned seed = params.GetAsInt("biasNoiseSeed");
	boost::mt19937 randomNumberEngine(seed);
    boost::normal_distribution<double> normalDistribution(m_biasMean,m_biasWidth);
    m_gaussianNoiseGenerator = new RANDOM_GAUSS(randomNumberEngine,normalDistribution);

    //Readout noise is always 7.2ADU for the full frame, so use a different Gaussian generator
    boost::normal_distribution<double> normalDistribution_fullFrame(m_biasMean,7.2);
    m_gaussianNoiseGenerator_fullFrame = new RANDOM_GAUSS(randomNumberEngine,normalDistribution_fullFrame);

//    m_applyAnalogChainRandomNoise = params.GetAsBool("applyAnalogChainRandomNoise");
//	  unsigned seed2 = params.GetAsInt("analogChainRandomNoiseSeed");
//	  boost::mt19937 randomNumberEngine2(seed2);
//    boost::poisson_distribution<int,double> poissonDistribution(params.GetAsDouble("meanAnalogChainRandomNoise"));
//    m_poissonNoiseGenerator = new RANDOM_POISSON(randomNumberEngine2,poissonDistribution);

	m_ccdNonlinearityExtrapSlope = 0.;
	m_ccdNonlinearityExtrapIntercept = 0.;

}

void BiasGenerator::doBegin(Data* data, bool fullFrame) {

	static bool firstCall = true;

	m_aduNonLinearity = false;
	if (m_doCcdNonLinearity) {
		if (fullFrame || data->getSerialReadRate()>165.) {
			m_ccdLinearityExtension = "[REF_APP_CCDLinearisation230]";
			readCcdNonLinearity<RefAppCcdlinearisation230>(m_ccdLinearityFilename+m_ccdLinearityExtension);
		} else {
			m_ccdLinearityExtension = "[REF_APP_CCDLinearisation100]";
			readCcdNonLinearity<RefAppCcdlinearisation100>(m_ccdLinearityFilename+m_ccdLinearityExtension);
		}
	}

	//Read in the empirical bias frame
    if (m_empiricalBiasFrame) {
    	if (m_biasFrameFilename[m_biasFrameFilename.length()-10] != 'V') {
    		throw runtime_error("name of REF_APP_BiasFrame filename must end with Vxxxx.fits");
    	}
    	//Get the CCD temperature
    	double ccdTemperature = SatelliteData::kDefaultCcdTemperature;
    	if (data->hasSatelliteData()) {
    		ccdTemperature = data->getSatelliteData(0,false,false)->getCcdTemperature();
    	}
    	readBiasFrame(fullFrame ? 230. : data->getSerialReadRate(), data->redundantHardware(), ccdTemperature);
    }

    cout << "mean bias offset = " << m_biasMean << ", mean RON = " << m_biasWidth << endl;

    data->setBiasOffset(m_biasMean); // To be used to set PIX_DATA_OFFSET in the image metadata

    if (firstCall) {

    	if (m_doCcdNonLinearity) data->setGainNonLinearity(); // To be used to set the NonLinCorr keyword in SCI_RAW_Subarray and SCI_RAW_FullArray

    	readGainCorrection(data->getGainFilename()); //if statement to avoid executing twice unnecessarily
        data->initializeBadPixelMap(data->getGainFilename());

		ofstream referenceFilesList("reference_files.txt", ios::app);
		referenceFilesList << data->getGainFilename() << endl;
		referenceFilesList << m_ccdLinearityFilename+m_ccdLinearityExtension << endl;
		referenceFilesList << m_biasFrameFilename << endl;
		referenceFilesList.close();

		data->getIdealLightCurveParams()->m_biasOffset = m_biasMean;
		data->getIdealLightCurveParams()->m_readNoise = m_biasWidth;

		firstCall = false;

	}

}

void BiasGenerator::readBiasFrame(double serialReadRate, bool redundantHardware, double ccdTemperature) {

	cout << "Bias frame file: " << (string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_biasFrameFilename) << endl;

	//Open the bias frame fits file and its extensions
	RefAppBiasframe * biasFrame = new RefAppBiasframe(string(getenv("CHEOPS_TESTDATA"))+"/resources/"+m_biasFrameFilename,"READONLY");
	RefAppBiasdarkleftframe * darkLeft = Open_RefAppBiasdarkleftframe(biasFrame);
	RefAppBiasdarkrightframe * darkRight = Open_RefAppBiasdarkrightframe(biasFrame);
	RefAppBiasdarktopframe * darkTop = Open_RefAppBiasdarktopframe(biasFrame);
	RefAppBiasblankleftframe * blankLeft = Open_RefAppBiasblankleftframe(biasFrame);
	RefAppBiasblankrightframe * blankRight = Open_RefAppBiasblankrightframe(biasFrame);
	RefAppBiasoverscanleftframe * overscanLeft = nullptr;
	RefAppBiasoverscanrightframe * overscanRight = nullptr;
	if (redundantHardware) {
		overscanRight = Open_RefAppBiasoverscanrightframe(biasFrame);
	} else {
		overscanLeft = Open_RefAppBiasoverscanleftframe(biasFrame);
	}
	RefAppBiasoverscantopframe * overscanTop = Open_RefAppBiasoverscantopframe(biasFrame);
	RefAppBiasframemetadata * biasMetadata = Open_RefAppBiasframemetadata(biasFrame);
	RefAppBiasoffset * biasOffsets = Open_RefAppBiasoffset(biasFrame);

	//Determine the layer of the fits image cube to read given the readout channel (main or redundant) and frequency (100 or 230 kHz)
	//The layers are defined in the REF_APP_BiasFrameMetadata extension
	int iBias=-1;
	int iRON=-1;
	int i=0;
	while (biasMetadata->ReadRow()) {
		//cout << biasMetadata->getCellRoFrequ()/1000. << " " << serialReadRate << " " << (biasMetadata->getCellRoHw() == "redundant")
		//	 << " " << redundantHardware << " " << biasMetadata->getCellDataType() << endl;
		if (lround(biasMetadata->getCellRoFrequ()/1000.) == lround(serialReadRate) &&
				(biasMetadata->getCellRoHw() == "redundant") == redundantHardware) {
			if (biasMetadata->getCellDataType()=="BIAS") iBias = i;
			if (biasMetadata->getCellDataType()=="RON") iRON = i;
		}
		i++;
	}
	if (iBias==-1 or iRON==-1) {
		throw runtime_error("No row in REF_APP_BiasFrameMetadata corresponding to RO frequency "+to_string(serialReadRate)+" and "
				+(redundantHardware?"redundant":"main")+" RO hardware");
	}

	//Fill the bias offset and RON arrays with values read from the file
	double bias_mean = 0.;
	double ron_mean = 0.;
	unsigned nPixels = 0;
	for (int i=0; i<Image::kXTotal; i++) {
		for (int j=0; j<Image::kYTotal; j++) {
			double bias = 0.;
			double ron = 0.;
			if (j<Image::kYDim) {
				if (i<Image::kNOverscanCols) {
					if (redundantHardware) {
						bias = (*overscanRight)[iBias][j][i];
						ron = (*overscanRight)[iRON][j][i];
					} else {
						bias = (*overscanLeft)[iBias][j][i];
						ron = (*overscanLeft)[iRON][j][i];
					}
				} else if (i<Image::kNOverscanCols+Image::kNBlankCols) {
					bias = (*blankLeft)[iBias][j][i-Image::kNOverscanCols];
					ron = (*blankLeft)[iRON][j][i-Image::kNOverscanCols];
				} else if (i<Image::kLeftMargin) {
					bias = (*darkLeft)[iBias][j][i-Image::kNOverscanCols-Image::kNBlankCols];
					ron = (*darkLeft)[iRON][j][i-Image::kNOverscanCols-Image::kNBlankCols];
				} else if (i>=Image::kLeftMargin+Image::kXDim+Image::kNDarkCols) {
					bias = (*blankRight)[iBias][j][i-Image::kLeftMargin-Image::kXDim-Image::kNDarkCols];
					ron = (*blankRight)[iRON][j][i-Image::kLeftMargin-Image::kXDim-Image::kNDarkCols];
				} else if (i>=Image::kLeftMargin+Image::kXDim) {
					bias = (*darkRight)[iBias][j][i-Image::kLeftMargin-Image::kXDim];
					ron = (*darkRight)[iRON][j][i-Image::kLeftMargin-Image::kXDim];
				} else {
					bias = (*biasFrame)[iBias][j][i-Image::kLeftMargin];
					ron = (*biasFrame)[iRON][j][i-Image::kLeftMargin];
				}
			} else if (i>=Image::kLeftMargin && i<Image::kXDim+Image::kLeftMargin){
				if (j>=Image::kYDim+Image::kNDarkRows) {
					bias = (*overscanTop)[iBias][j-Image::kYDim-Image::kNDarkRows][i-Image::kLeftMargin];
					ron = (*overscanTop)[iRON][j-Image::kYDim-Image::kNDarkRows][i-Image::kLeftMargin];
				} else {
					bias = (*darkTop)[iBias][j-Image::kYDim][i-Image::kLeftMargin];
					ron = (*darkTop)[iRON][j-Image::kYDim][i-Image::kLeftMargin];
				}
			}
			m_bias[i][j] = bias > 0. ? bias : 0.;
			m_RON[i][j] = ron > 0 ? ron : 0.;
			bias_mean += bias;
			ron_mean += ron;
			nPixels++;
		}
	}

	m_biasMean = bias_mean/nPixels;
	m_biasWidth = ron_mean/nPixels;

	//Read in the bias offsets from the REF_APP_BiasOffset extension according to readout channel (main or redundant) and frequency (100 or 230 kHz)
	double biasOffset = 0.;
	while (biasOffsets->ReadRow()) {
		//cout << biasOffsets->getCellRoFrequ()/1000. << " " << serialReadRate << " " << (biasOffsets->getCellRoHw() == "redundant")
		//	 << " " << redundantHardware << " " << fabs(biasOffsets->getCellCcdTemp()+273.15 - ccdTemperature) << endl;
		if (lround(biasOffsets->getCellRoFrequ()/1000.) == lround(serialReadRate) &&
				(biasOffsets->getCellRoHw() == "redundant") == redundantHardware &&
				fabs(biasOffsets->getCellCcdTemp()+273.15 - ccdTemperature) < 2.5) {
			biasOffset = biasOffsets->getCellBiasOffset();
		}
	}
	if (biasOffset==0.) {
		throw runtime_error("No row in REF_APP_BiasOffset corresponding to RO frequency "+to_string(serialReadRate)+", "
				+(redundantHardware?"redundant":"main")+" RO hardware and CCD temperature "+to_string(ccdTemperature));
	}

	//Scale the bias frame according to the bias offsets from in flight measurements
	for (int i=0; i<Image::kXTotal; i++) {
		for (int j=0; j<Image::kYTotal; j++) {
			m_bias[i][j] *= biasOffset/m_biasMean;
		}
	}
	//cout << "m_biasMean = " << m_biasMean << ", biasOffset = " << biasOffset << ", read rate = " << serialReadRate << ", redundant = " << redundantHardware << endl;
	m_biasMean = biasOffset;

	delete darkLeft;
	delete darkRight;
	delete darkTop;
	delete blankLeft;
	delete blankRight;
	delete overscanLeft;
	delete overscanTop;
	delete biasMetadata;
	delete biasOffsets;
	delete biasFrame;

}

void BiasGenerator::readGainCorrection(string gainCorrectionFilename) {

	cout << "Gain correction file: " << (string(getenv("CHEOPS_SW"))+"/resources/"+gainCorrectionFilename) << endl;

	RefAppGaincorrection * gainCorrection_file = new RefAppGaincorrection(string(getenv("CHEOPS_SW"))+"/resources/"+gainCorrectionFilename,"READONLY");

	m_nominalGain = gainCorrection_file->getKeyGainNom();

	m_gainFormula.clear();
	while (gainCorrection_file->ReadRow()) {
		double constant = gainCorrection_file->getCellFactor();
		double vodExponent = gainCorrection_file->getCellExpVod();
		double vrdExponent = gainCorrection_file->getCellExpVrd();
		double vogExponent = gainCorrection_file->getCellExpVog();
		double vssExponent = gainCorrection_file->getCellExpVss();
		double tempExponent = gainCorrection_file->getCellExpTemp();
		m_gainFormula.push_back(GainFormulaLine(constant,vodExponent,vrdExponent,vogExponent,vssExponent,tempExponent));
	}

	delete gainCorrection_file;

}

template <class T> void BiasGenerator::readCcdNonLinearity(string ccdNonLinearityFilename) {

	cout << "CCD non-linearity file: " << (string(getenv("CHEOPS_TESTDATA"))+"/resources/"+ccdNonLinearityFilename) << endl;

	T * ccdNonLinearity_file = new T(string(getenv("CHEOPS_TESTDATA"))+"/resources/"+ccdNonLinearityFilename,"READONLY");

	//Read in the polynomial coefficients for each spline interval, together with the boundaries for each interval
	double offset = 0.;
	unsigned i = 0;
	m_nIntervals = 0;
	while (ccdNonLinearity_file->ReadRow()) {
		m_ccdNonlinearityBoundaries[i] = ccdNonLinearity_file->getCellBoundary();
		if (!ccdNonLinearity_file->isNullCoef0()) {
			m_nIntervals++;
			m_ccdNonlinearityCoeffs[0][i] = ccdNonLinearity_file->getCellCoef0();
			m_ccdNonlinearityCoeffs[1][i] = ccdNonLinearity_file->getCellCoef1();
			m_ccdNonlinearityCoeffs[2][i] = ccdNonLinearity_file->getCellCoef2();
			m_ccdNonlinearityCoeffs[3][i] = ccdNonLinearity_file->getCellCoef3();
		}

		//Shift all elements of the spline function by the constant term of the first interval, so that the curve goes through (0,0)
		if (i==0) offset = m_ccdNonlinearityCoeffs[0][0];
		m_ccdNonlinearityCoeffs[0][i] -= offset;

		i++;
	}

//	for (unsigned j=0; j<m_nIntervals; j++) {
//		cout << m_ccdNonlinearityBoundaries[j] << " " << m_ccdNonlinearityBoundaries[j+1] << " ";
//		for (unsigned k=0; k<4; k++) {
//			cout << m_ccdNonlinearityCoeffs[k][j] << " ";
//		}
//		cout << endl;
//	}

	delete ccdNonLinearity_file;

	//Support for older REF_APP_CCDLinearisation using NLC correction based on ADU rather than nElectrons
	if (m_aduNonLinearity) {

		//Cubic function for last interval has a maximum at this value, resulting in two possible solutions on either side of the maximum
		//Avoid this for setting the upper edge of the range for the last interval to the location of the maximum
		m_ccdNonlinearityBoundaries[28] = 62500;

		//The spline function is not valid for uncorrected values which correspond to corrected values > 2^16=65536
		//i.e. it is not valid for uncorrected values > 58000 (upper edge of boundary of last interval).
		//CHEOPSim should be able to output uncorrected ADU values up to 2^16. In order to achieve this, an extrapolation
		//of the spline function must be made to cover uncorrected ADU values in the range 58000 to 65536.
		//Up to uncorrected ADU=71000, it is assumed that the polynomial from the last interval is valid.
		//Beyond 71000, the following code defines a linear extrapolation as a tangent to the polynomial
		double adu1 = 70000;
		double adu2 = 71000;
		applyCcdNonLinearity(adu1);
		applyCcdNonLinearity(adu2);
		m_ccdNonlinearityExtrapSlope = (adu2-adu1)/1000.;
		m_ccdNonlinearityExtrapIntercept = adu1 - m_ccdNonlinearityExtrapSlope*70000.;

	}

	//The following can be uncommented to generate a table of values before and after the inverse correction for debugging purposes
//	for (double nElectrons = 0.; nElectrons<160000.; nElectrons+=100.) {
//		double nElectrons_before = nElectrons;
//		double nElectrons_after = nElectrons;
//		applyCcdNonLinearity(nElectrons_after);
//		cout << nElectrons_before << " " << nElectrons_after << endl;
//	}

}

void BiasGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	Image* image = data->getImages().back();

	int xmin = image->getXOffset();
	int xmax = image->getXOffset()+image->getXDim();
	int ymin = image->getYOffset();
	int ymax = image->getYOffset()+image->getYDim();

	//Get the CCD temperature
	double ccdTemperature = SatelliteData::kDefaultCcdTemperature;
	if (data->hasSatelliteData()) {
		ccdTemperature = data->getSatelliteData(timeStep)->getCcdTemperature();
	}

	//Nominal electronic gain
	double gain = m_nominalGain;
	double sum = 0.;

	//Analog electronics stability: Gain dependency on bias voltages
	if (m_applyAnalogElectronicsStability) {
		double timeSinceVisitStart = data->getTimeConfiguration().getTimeSinceVisitStart(timeStep);
		double vod = data->getVoltageWithDrift(timeSinceVisitStart,SatelliteData::VOD);
		double vrd = data->getVoltageWithDrift(timeSinceVisitStart,SatelliteData::VRD);
		double vog = data->getVoltageWithDrift(timeSinceVisitStart,SatelliteData::VOG);
		double vss = data->getVoltageWithDrift(timeSinceVisitStart,SatelliteData::VSS);
		double vod_nominal = data->getNominalVoltage(SatelliteData::VOD);
		double vrd_nominal = data->getNominalVoltage(SatelliteData::VRD);
		double vog_nominal = data->getNominalVoltage(SatelliteData::VOG);
		double vss_nominal = data->getNominalVoltage(SatelliteData::VSS);
		double temp_nominal = data->getNominalVoltage(SatelliteData::TEMP);
		//cout << setprecision(10) << vod << " " << vrd << " " << vog << " " << vss << " " << ccdTemperature-273.15 << " " << vod_nominal << " " << vrd_nominal << " " << vog_nominal << " " << vss_nominal << " " << temp_nominal << endl;
		//cout << (vod-vss)-(vod_nominal-vss_nominal) << " " << (vrd-vss)-(vrd_nominal-vss_nominal) << " " << (vog-vss)-(vog_nominal-vss_nominal) << " " << vss-vss_nominal << " " << ccdTemperature-273.15-temp_nominal << endl;
		for (unsigned i=0; i<m_gainFormula.size(); i++) {
			sum += m_gainFormula[i].m_constant *
					pow((vod-vss)-(vod_nominal-vss_nominal),m_gainFormula[i].m_vodExponent) *
					pow((vrd-vss)-(vrd_nominal-vss_nominal),m_gainFormula[i].m_vrdExponent) *
					pow((vog-vss)-(vog_nominal-vss_nominal),m_gainFormula[i].m_vogExponent) *
					pow(vss-vss_nominal,m_gainFormula[i].m_vssExponent) *
					pow(ccdTemperature-273.15-temp_nominal,m_gainFormula[i].m_tempExponent);
			//cout << m_gainFormula[i].m_constant << " " << m_gainFormula[i].m_vodExponent << " " << m_gainFormula[i].m_vrdExponent << " "
			//	   << m_gainFormula[i].m_vogExponent << " " << m_gainFormula[i].m_vssExponent << " " << m_gainFormula[i].m_tempExponent << endl;
		}
	}
	gain *= (1 + sum);
	//cout << m_nominalGain << " " << sum << " " << gain << endl;

	boost::math::tools::stats<double> ronStats;
	double pixelADU_max = 0.;
	for (int ix=xmin; ix<xmax; ix++) {
		for (int iy=ymin; iy<ymax; iy++) {
			generateBiasNoise(image,ix,iy,gain,pixelADU_max,fullFrame);
    		ronStats.add(image->getPixelValue(ix,iy));
		}
		//Generate noise also for dark and overscan rows
		for (int iy=Image::kYDim; iy<Image::kYDim+Image::kNDarkRows+Image::kNOverscanRows; iy++) {
			generateBiasNoise(image,ix,iy,gain,pixelADU_max,fullFrame);
		}
	}
	//Generate noise also for dark, blank and overscan columns
	for (int iy=ymin; iy<ymax; iy++) {
		for (int ix=-Image::kLeftMargin; ix<0; ix++) {
			generateBiasNoise(image,ix,iy,gain,pixelADU_max,fullFrame);
		}
		for (int ix=Image::kXDim; ix<Image::kXDim+Image::kNDarkCols+Image::kNBlankCols; ix++) {
			generateBiasNoise(image,ix,iy,gain,pixelADU_max,fullFrame);
		}
	}
    //cout << "after BiasGenerator: mean=" << ronStats.mean()
    //     << ", standard deviation=" << sqrt(ronStats.variance1()) << endl;

	//Set the truth information
	image->getTruthData()->setGain(gain);

	//cout << "max ADU count in image = " << round(pixelADU_max) << endl;

}

void BiasGenerator::generateBiasNoise(Image* image, int ix, int iy, double gain, double & pixelADU_max, bool fullFrame) const {

	//Get the number of electrons in the pixel
	double nElectrons = image->getPixelValue(ix,iy);

	//Apply the inverse of the CCD non-linearity correction modeled using a quadratic spline (parameterized in terms of nElectrons)
	if (!m_aduNonLinearity && m_doCcdNonLinearity) {
		//Do not apply for overscan and blank columns which do not correspond to physical pixels (see WG B3 minutes 7th Nov 2018)
		if (ix>=-Image::kNDarkCols && ix<Image::kXDim+Image::kNDarkCols) applyCcdNonLinearity(nElectrons);
	}

	//Convert electrons to ADU by applying the gain
	double pixelADU = nElectrons*gain;
	//if (iy==512) cout << nElectrons << " " << gain << " " << pixelADU << " ";

	//Apply the inverse of the CCD non-linearity correction modeled using a cubic spline (parameteized in terms of ADU)
	if (m_aduNonLinearity && m_doCcdNonLinearity) {
		//Do not apply for overscan and blank columns which do not correspond to physical pixels (see WG B3 minutes 7th Nov 2018)
		if (ix>=-Image::kNDarkCols && ix<Image::kXDim+Image::kNDarkCols) applyCcdNonLinearity(pixelADU);
	}

	//Generate bias offset and read noise and add it to the pixel value
	if (!m_empiricalBiasFrame) {
		if (fullFrame) {
			pixelADU += (*m_gaussianNoiseGenerator_fullFrame)();
		} else {
			pixelADU += (*m_gaussianNoiseGenerator)();
		}
	} else {
		//if (iy==1027) cout << ix << " " << iy << " " << m_bias[ix+Image::kLeftMargin][iy] << " " << m_RON[ix+Image::kLeftMargin][iy] << endl;
		if (fullFrame) { //using separate RNG for full frame allows sub-array content to be the same whether or not full frame is present
			m_gaussianNoiseGenerator_fullFrame->distribution().param(boost::normal_distribution<double>::param_type(m_bias[ix+Image::kLeftMargin][iy],m_RON[ix+Image::kLeftMargin][iy]));
			pixelADU += (*m_gaussianNoiseGenerator_fullFrame)();
		} else {
			m_gaussianNoiseGenerator->distribution().param(boost::normal_distribution<double>::param_type(m_bias[ix+Image::kLeftMargin][iy],m_RON[ix+Image::kLeftMargin][iy]));
			pixelADU += (*m_gaussianNoiseGenerator)();
		}
		//cout << m_bias[ix+Image::kLeftMargin][iy] << " " << m_RON[ix+Image::kLeftMargin][iy] << " " << pixelADU << endl;
	}

	//Generate analog chain random noise and add it to the pixel value
	//if (m_applyAnalogChainRandomNoise) pixelADU += (*m_poissonNoiseGenerator)();

	//Determine the pixel with the highest ADU value in the image (used for diagnostic printout only)
	if (pixelADU > pixelADU_max) pixelADU_max = pixelADU;

	//Update the image
	image->setPixelValue(ix,iy,pixelADU);

}

void BiasGenerator::applyCcdNonLinearity(double& nElectrons) const {

	if (m_aduNonLinearity) {
		//For corrected ADU>71000, use a linear extrapolation since the polynomial is no longer valid.
		//(see comments in readCcdNonLinearity)
		if (nElectrons>71000) {
			nElectrons = m_ccdNonlinearityExtrapSlope * nElectrons + m_ccdNonlinearityExtrapIntercept;
			return;
		}
	}

	//Determine which interval of the spline function to use
	//The BOUNDARIES column of the REF_APP_CCDLinearisation file (m_ccdNonlinearityBoundaries)
	//are the interval boundaries in uncorrected electrons. We are starting from corrected electrons so need
	//the boundaries in corrected electrons. These correspond to the constant term of the polynomial
	//in each interval (m_ccdNonlinearityCoeffs[0][interval])
	unsigned interval = 0;
	while (nElectrons > m_ccdNonlinearityCoeffs[0][interval+1] && interval<m_nIntervals-1) {interval++;};

	//Get the polynomial coefficients corresponding to the current interval
	double d = m_ccdNonlinearityCoeffs[0][interval];
	double c = m_ccdNonlinearityCoeffs[1][interval];
	double b = m_ccdNonlinearityCoeffs[2][interval];
	double a = m_ccdNonlinearityCoeffs[3][interval];

	double uncorrectedElectrons = -1.;

	//We need to solve the equation ax^3 + bx^2 + cx + d = nElectrons
	//i.e. ax^3 + bx^2 + cx + (d - nElectrons) = 0
	d -= nElectrons;

	if (a == 0.) {

		//solve bx^2 + cx + d = 0
		if (interval==m_nIntervals-1 && (c*c - 4*b*d)<0) { //polynomial for last interval is only valid up to around 400k electrons
			uncorrectedElectrons = nElectrons;
		} else {
			double root = sqrt(c*c - 4*b*d);
			double solution1 = (-c + root)/(2*b);
			double solution2 = (-c - root)/(2*b); //It seems that this is never the valid solution, in which case the if block below could be simplified
			//if (interval>3) cout << interval << " " << a << " " << b << " " << c << " " << d << " " << b*solution1*solution1 + c*solution1 + d << " " << b*solution2*solution2 + c*solution2 + d << endl;
			solution1 += m_ccdNonlinearityBoundaries[interval];
			solution2 += m_ccdNonlinearityBoundaries[interval];
			//if (interval>3) cout << " " << nElectrons << " " << solution1 << " " << solution2 << endl;
			if (solution1>m_ccdNonlinearityBoundaries[interval]-1 &&
					solution1<m_ccdNonlinearityBoundaries[interval+1]+1) {
				uncorrectedElectrons = solution1;
			} else if (solution2>m_ccdNonlinearityBoundaries[interval]-1 &&
					solution2<m_ccdNonlinearityBoundaries[interval+1]+1) {
				uncorrectedElectrons = solution2;
			} else if (interval == m_nIntervals-1) {
				double lastBound = m_ccdNonlinearityBoundaries[m_nIntervals];
				if (solution1 > lastBound && solution2 > lastBound) {
					uncorrectedElectrons = min(solution1,solution2);
				} else if (solution1 > lastBound && solution2 < lastBound) {
					uncorrectedElectrons = solution1;
				} else if (solution2 > lastBound && solution1 < lastBound) {
					uncorrectedElectrons = solution2;
				}
			} else {
				throw runtime_error("Error in BiasGenerator::applyCcdNonLinearity: solutions "
						+to_string(solution1)+" and "+to_string(solution2)+" are both outside the expected interval ("
						+to_string(m_ccdNonlinearityBoundaries[interval]-1)+","
						+to_string(m_ccdNonlinearityBoundaries[interval+1]+1));
			}
		}

	} else {

		//The equation is solved as described here:
		//https://en.wikipedia.org/wiki/Cubic_function#General_solution_to_the_cubic_equation_with_real_coefficients
		//First we determine the discriminant:
		double delta = 18*a*b*c*d - 4*pow(b,3)*d + b*b*c*c -4*a*pow(c,3) - 27*a*a*d*d;
		//cout << nElectrons << " " << interval << " " << delta << endl;

		if (delta<0) {

			//If the discriminant is negative, then there is a unique solution and we use the general formula from
			//https://en.wikipedia.org/wiki/Cubic_function#General_solution_to_the_cubic_equation_with_real_coefficients
			double delta0 = b*b - 3*a*c;
			double delta1 = 2*pow(b,3) - 9*a*b*c + 27*a*a*d;
			double C = pow((delta1 + sqrt(-27*a*a*delta))/2,1./3.);
			uncorrectedElectrons = (b + C + delta0/C)/(-3*a) + m_ccdNonlinearityBoundaries[interval];

		} else {

			//If the discriminant is positive, then there are three solutions. In this case the general formula involves taking
			//the cube root of a complex number. To avoid this, we use the trigonometric solution by reduction to a depressed cubic
			//as described here:
			//https://en.wikipedia.org/wiki/Cubic_function#Trigonometric_and_hyperbolic_solutions
			double p = (3*a*c - b*b)/(3*a*a);
			double q = (2*pow(b,3) - 9*a*b*c + 27*a*a*d)/(27*pow(a,3));
			for (unsigned k=0; k<3; k++) {
				double t = 2*sqrt(-p/3) * cos(acos( (3*q/(2*p)) * sqrt(-3/p) )/3 - 2*M_PI*k/3);
				double root = t - b/(3*a) + m_ccdNonlinearityBoundaries[interval];
				//cout << root << " " << m_ccdNonlinearityBoundaries[interval] << " " << m_ccdNonlinearityBoundaries[interval+1] << endl;
				//Of the three solutions (3 values of k), take the one which falls within the current interval
				if (root>m_ccdNonlinearityBoundaries[interval]-1 &&
						root<m_ccdNonlinearityBoundaries[interval+1]+1) {
					uncorrectedElectrons = root;
					break;
				}
			}

		}

		if (uncorrectedElectrons < m_ccdNonlinearityBoundaries[interval]-1 ||
			uncorrectedElectrons > m_ccdNonlinearityBoundaries[interval+1]+1) {
			throw runtime_error("Error in BiasGenerator::applyCcdNonLinearity: result "
					+to_string(uncorrectedElectrons)+" is outside the expected interval ("
					+to_string(m_ccdNonlinearityBoundaries[interval]-1)+","
					+to_string(m_ccdNonlinearityBoundaries[interval+1]+1));
		}

	}

	//cout << nElectrons << " " << uncorrectedElectrons << endl << endl;
	if (std::isnan(uncorrectedElectrons)) throw runtime_error("Error in BiasGenerator::applyCcdNonLinearity: result is nan");

	nElectrons = uncorrectedElectrons;

}

void BiasGenerator::correctCcdNonLinearity(double& nElectrons) const {

	unsigned interval = 0;
	while (nElectrons > m_ccdNonlinearityBoundaries[interval+1] && interval<m_nIntervals-1) {interval++;};

	//Get the polynomial coefficients corresponding to the current interval
	double d = m_ccdNonlinearityCoeffs[0][interval];
	double c = m_ccdNonlinearityCoeffs[1][interval];
	double b = m_ccdNonlinearityCoeffs[2][interval];
	double a = m_ccdNonlinearityCoeffs[3][interval];

	double x = nElectrons - m_ccdNonlinearityBoundaries[interval];
	nElectrons = a*pow(x,3) + b*x*x + c*x + d;

}
