/*
 * PSFGenerator.cxx
 *
 *  Created on: 10 Feb 2014
 *      Author: futyand
 */

#include <iostream>
#include <fstream>
#include <iomanip>

#include "boost/filesystem.hpp"

#include "REF_APP_WhitePSF.hxx"
#include "REF_APP_WhiteCCDLocationPSF.hxx"
#include "REF_APP_WhiteCCDLocationPSFMetadata.hxx"
#include "REF_APP_ColouredPSF.hxx"
#include "REF_APP_OversampledWhitePSF.hxx"
#include "REF_APP_OversampledColouredPSF.hxx"

#include "detector/include/FlatFieldGenerator.hxx"
#include "detector/include/CosmicRayGenerator.hxx"
#include "simulator/include/CommonTools.hxx"

#include "PSFGenerator.hxx"

void PSFGenerator::initialize(const ModuleParams& params) {

	m_psfFilename = params.GetAsString("PSF_filename");
	m_oversampleJitter = params.GetAsBool("oversampleJitter");
	string targetStarPositioning = params.GetAsString("targetStarPositioning");
	string backgroundStarPositioning = params.GetAsString("backgroundStarPositioning");
	m_convertFluxToPhotons = params.GetAsBool("convertFluxToPhotons");
	m_monochromaticPSF = params.GetAsBool("monochromaticPSF");
	m_monochromaticPSFwavelength = params.GetAsInt("monochromaticPSFwavelength");
	string thermalMap = params.GetAsString("thermalMap");

	if (thermalMap=="fixed") {
		m_thermalMap = fixed;
	} else if (thermalMap=="cold") {
		m_thermalMap = cold;
	} else if (thermalMap=="hot1") {
		m_thermalMap = hot1;
	} else if (thermalMap=="hot2") {
		m_thermalMap = hot2;
	} else if (thermalMap=="breathing") {
		m_thermalMap = breathing;
	} else {
		throw runtime_error("Error in PSFGenerator::initialize: thermalMap must be minusTen, cold, hot1, hot2 or breathing");
	}

	if (targetStarPositioning=="Oversampling") {
		m_targetStarPositioning = Oversampling;
	} else if (targetStarPositioning=="Interpolation") {
		m_targetStarPositioning = Interpolation;
	} else if (targetStarPositioning=="SnapToGrid") {
		m_targetStarPositioning = SnapToGrid;
	} else {
		throw runtime_error("Error in PSFGenerator::initialize: targetStarPositioning must be Oversampling, Interpolation or SnapToGrid");
	}

	if (backgroundStarPositioning=="Oversampling") {
		m_backgroundStarPositioning = Oversampling;
	} else if (backgroundStarPositioning=="Interpolation") {
		m_backgroundStarPositioning = Interpolation;
	} else if (backgroundStarPositioning=="SnapToGrid") {
		m_backgroundStarPositioning = SnapToGrid;
	} else {
		throw runtime_error("Error in PSFGenerator::initialize: backgroundStarPositioning must be Oversampling, Interpolation or SnapToGrid");
	}

	//Initialize the PSF arrays
	for (int k=0; k<2; k++) {
		for (int j = 0; j<kOversampledPsfSize; j++) {
			for (int i = 0; i<kOversampledPsfSize; i++) {
				m_oversampledPsf[i][j][k] = 0.;
			}
		}
		for (int j = 0; j<kPsfSize; j++) {
			for (int i = 0; i<kPsfSize; i++) {
				m_psf[i][j][k] = 0.;
			}
		}
	}

}

void PSFGenerator::doBegin(Data *data, bool fullFrame) {

	if (m_psfFilename == "FLAT") {

		if (m_targetStarPositioning==Oversampling || m_backgroundStarPositioning==Oversampling) {
			for (int j = 0; j<kOversampledPsfSize; j++) {
				for (int i = 0; i<kOversampledPsfSize; i++) {
					m_oversampledPsf[i][j][0] = 1./static_cast<double>(kOversampledPsfSize*kOversampledPsfSize);
				}
			}
		}

		for (int j = 0; j<kPsfSize; j++) {
			for (int i = 0; i<kPsfSize; i++) {
				m_psf[i][j][0] = 1./static_cast<double>(kPsfSize*kPsfSize);
			}
		}

	} else {

		// Construct filename for non-oversampled PSF and determine whether or not it already exists
		string psfFilePath = string(getenv("CHEOPS_TESTDATA"))+"/resources/";
		string prefix = m_psfFilename.substr(0,m_psfFilename.find_last_of("."));
		string suffix = m_psfFilename.substr(m_psfFilename.find_last_of("."));
		string psfFilename_rebin = prefix+"_rebin"+suffix;
		if (m_psfFilename.find("REF_APP_WhitePSF") != std::string::npos || m_psfFilename.find("REF_APP_WhiteCCDLocationPSF") != std::string::npos) {
			psfFilename_rebin = m_psfFilename;
		}
		bool rebinned_exists = boost::filesystem::exists(psfFilePath+psfFilename_rebin);

		//For monochromatic case, determine the weight for each wavelength bin
		if (m_monochromaticPSF) {
			getWavelengthWeights(data);
		}

		//Read oversampled PSFs if they are needed
		if (m_targetStarPositioning==Oversampling || m_backgroundStarPositioning==Oversampling || !rebinned_exists) {
			readPSFs(true,data->getSubarrayDimensions());
		}

		//Read non-oversampled PSFs if they exist
		if (rebinned_exists) {

			readPSFs(false,data->getSubarrayDimensions());

		} else {

			//For case of user uploaded PSF, instead of reading the rebinned PSF from a file, do the rebinning on the fly
			for (int k=0; k<((m_thermalMap==breathing)?2:1); k++) {
				for (int j = 0; j<kPsfSize; j++) {
					for (int i = 0; i<kPsfSize; i++) {
						m_psf[i][j][k]=0.;
						for (int ii=0; ii<10; ii++) {
							for (int jj=0; jj<10; jj++) {
								m_psf[i][j][k] += m_oversampledPsf[i*10+ii][j*10+jj][k];
							}
						}
					}
				}
			}

		}

//		//Special 4 pixel square PSF for testing spatial interpolation
//		for (int j = 0; j<kPsfSize; j++) {
//			for (int i = 0; i<kPsfSize; i++) {
//				m_psf[i][j][0] = ((i==99||i==100) && (j==99||j==100)) ? 1. : 0.;
//			}
//		}

	}

}

void PSFGenerator::getWavelengthWeights(Data* data) {

	//Calculate weights taking into account the BB spectrum of the target star, QE and throughput,
	//as a function of wavelength in 0.5nm steps
	WavelengthDependence * wavelengthDependence = data->getWavelengthDependence();
	double effectiveTemperature = 5660.; //Default to effective temperature for G5 star if there are no stars in the FOV
	if (data->getFieldOfView()->getStars().size()>0) {
		effectiveTemperature = (data->getFieldOfView()->getStars()[0])->getEffectiveTemperature();
	}
	wavelengthDependence->generateWeights(effectiveTemperature);

	//For each wavelength for which a PSF exists in the input REF_APP_(Oversampled)ColouredPSF file,
	//define an associated wavelength range and define the weight for that PSF as the sum over weights
	//within that wavelength range
	double weightSum = 0.;
	for (unsigned iWavelength=0; iWavelength<kNpsfWavelength; iWavelength++) {
		double wavelength = kPsfFirstWavelength + kPsfWavelengthBinWidth*iWavelength;
		double wavelengthLow = wavelength - kPsfWavelengthBinWidth/2;
		double wavelengthHigh = wavelength + kPsfWavelengthBinWidth/2;
		//For the first wavelength bin, set the lower end of the range to the lowest wavelength considered for stellar flux
		if (iWavelength==0) wavelengthLow = Star::kWavelengthLowBound;
		//For the last wavelength bin, set the upper end of the range to the highest wavelength considered for stellar flux
		if (iWavelength==kNpsfWavelength-1) wavelengthHigh = Star::kWavelengthLowBound+10*Star::kNWavelength;
		if (m_monochromaticPSFwavelength==0) {
			m_weight[iWavelength] = wavelengthDependence->getWeight(wavelengthLow,wavelengthHigh);
		} else {
			m_weight[iWavelength] = (m_monochromaticPSFwavelength>=wavelengthLow && m_monochromaticPSFwavelength<wavelengthHigh) ? 1. : 0.;
		}
		weightSum += m_weight[iWavelength];
	}

	//Normalize so that the sum over the weights is 1
	cout << "  Wavelength weighting for monochromatic PSFs:" << endl;
	for (unsigned iWavelength=0; iWavelength<kNpsfWavelength; iWavelength++) {
		m_weight[iWavelength] /= weightSum;
		cout << "    " << (kPsfFirstWavelength + kPsfWavelengthBinWidth*iWavelength) << " " << m_weight[iWavelength] << endl;
	}

}

void PSFGenerator::readPSFs(bool oversampled, Data::SubarrayDimensions subarrayDimensions) {

	string psfFilePath = string(getenv("CHEOPS_TESTDATA"))+"/resources/";
	string prefix = m_psfFilename.substr(0,m_psfFilename.find_last_of("."));
	string suffix = m_psfFilename.substr(m_psfFilename.find_last_of("."));
	string psfFilename_rebin = prefix+"_rebin"+suffix;
	if (m_psfFilename.find("REF_APP_WhitePSF") != std::string::npos || m_psfFilename.find("REF_APP_WhiteCCDLocationPSF") != std::string::npos) {
		psfFilename_rebin = m_psfFilename;
	}
	bool ccdLocation = (m_psfFilename.find("CCDLocation") != std::string::npos);
	if (ccdLocation && m_thermalMap==breathing) {
		throw runtime_error("Cannot select PSF breathing for CCD location dependent PSF");
	}
	unsigned iClosestLocation = 0;

	//Open the PSF fits file
	FitsDalImage<double> * psf;
	if (oversampled) {
		if (m_monochromaticPSF) {
			psf = new RefAppOversampledcolouredpsf(psfFilePath+m_psfFilename,"READONLY");
		} else {
			psf = new RefAppOversampledwhitepsf(psfFilePath+m_psfFilename,"READONLY");
		}
	} else {
		if (m_monochromaticPSF) {
			psf = new RefAppColouredpsf(psfFilePath+psfFilename_rebin,"READONLY");
		} else {
			if (!ccdLocation) {
				psf = new RefAppWhitepsf(psfFilePath+psfFilename_rebin,"READONLY");
			} else {
				psf = new RefAppWhiteccdlocationpsf(psfFilePath+psfFilename_rebin,"READONLY");
				RefAppWhiteccdlocationpsf * ccdLocationPsf = new RefAppWhiteccdlocationpsf(psfFilePath+psfFilename_rebin,"READONLY");
				RefAppWhiteccdlocationpsfmetadata * ccdLocationMetadata = Open_RefAppWhiteccdlocationpsfmetadata(ccdLocationPsf);
				unsigned targLocx = subarrayDimensions.m_targetLocationX;
				unsigned targLocy = subarrayDimensions.m_targetLocationY;
				double mindist = 999.;
				unsigned iLoc = 0;
				while (ccdLocationMetadata->ReadRow()) {
					double locx = ccdLocationMetadata->getCellXOffFullArray() + kPsfSize/2.;
					double locy = ccdLocationMetadata->getCellYOffFullArray() + kPsfSize/2.;
					double dist = sqrt((locx-targLocx)*(locx-targLocx) + (locy-targLocy)*(locy-targLocy));
					if (dist < mindist) {
						iClosestLocation = iLoc;
						mindist = dist;
					}
					iLoc++;
				}
				cout << "Closest PSF from REF_APP_WhiteCCDLocationPSF is image " << iClosestLocation << ", at distance " << mindist << " from target location (" << targLocx << "," << targLocy << ")"<< endl;
				delete ccdLocationMetadata;
				delete ccdLocationPsf;
			}
		}
	}

	ofstream referenceFilesList("reference_files.txt", ios::app);
	if (m_psfFilename == "PSF_SelectedMeasurement_12px_-10C.fits") {
		referenceFilesList << "CH_TU2017-05-01T13-22-00_REF_APP_WhitePSF_V0002.fits" << endl;
	} else {
		referenceFilesList << m_psfFilename << endl;
	}
	referenceFilesList.close();

	//Read in the PSF(s)
	if (m_monochromaticPSF) {

		//Perform weighted integral over wavelength
		for (int iWavelength=0; iWavelength<kNpsfWavelength; iWavelength++) {
			for (int j = 0; j<(oversampled?kOversampledPsfSize:kPsfSize); j++) {
				for (int i = 0; i<(oversampled?kOversampledPsfSize:kPsfSize); i++) {
					if (oversampled) {
						if (m_thermalMap==breathing) {
							m_oversampledPsf[i][j][0] += (*psf)[hot1][iWavelength][j][i] * m_weight[iWavelength];
							m_oversampledPsf[i][j][1] += (*psf)[hot2][iWavelength][j][i] * m_weight[iWavelength];
						} else {
							m_oversampledPsf[i][j][0] += (*psf)[m_thermalMap][iWavelength][j][i] * m_weight[iWavelength];
						}
					} else {
						if (m_thermalMap==breathing) {
							m_psf[i][j][0] += (*psf)[hot1][iWavelength][j][i] * m_weight[iWavelength];
							m_psf[i][j][1] += (*psf)[hot2][iWavelength][j][i] * m_weight[iWavelength];
						} else {
							m_psf[i][j][0] += (*psf)[m_thermalMap][iWavelength][j][i] * m_weight[iWavelength];
						}
					}
				}
			}
		}

	} else {

		for (int j = 0; j<(oversampled?kOversampledPsfSize:kPsfSize); j++) {
			for (int i = 0; i<(oversampled?kOversampledPsfSize:kPsfSize); i++) {
				if (oversampled) {
					if (m_thermalMap==breathing) {
						m_oversampledPsf[i][j][0] = (*psf)[hot1][j][i];
						m_oversampledPsf[i][j][1] = (*psf)[hot2][j][i];
					} else {
						m_oversampledPsf[i][j][0] = (*psf)[m_thermalMap][j][i];
					}
				} else {
					if (!ccdLocation) {
						if (m_thermalMap==breathing) {
							m_psf[i][j][0] = (*psf)[hot1][j][i];
							m_psf[i][j][1] = (*psf)[hot2][j][i];
						} else {
							m_psf[i][j][0] = (*psf)[m_thermalMap][j][i];
						}
					} else {
						m_psf[i][j][0] = (*psf)[iClosestLocation][j][i];
					}
				}
			}
		}

	}

	//Normalize the PSFs to 1
	for (int k=0; k<((m_thermalMap==breathing)?2:1); k++) {
		double psf_sum = 0.;
		for (int j = 0; j<(oversampled?kOversampledPsfSize:kPsfSize); j++) {
			for (int i = 0; i<(oversampled?kOversampledPsfSize:kPsfSize); i++) {
				if (oversampled) {
					psf_sum += m_oversampledPsf[i][j][k];
				} else {
					psf_sum += m_psf[i][j][k];
				}
			}
		}
		for (int j = 0; j<(oversampled?kOversampledPsfSize:kPsfSize); j++) {
			for (int i = 0; i<(oversampled?kOversampledPsfSize:kPsfSize); i++) {
				if (oversampled) {
					m_oversampledPsf[i][j][k]/=psf_sum;
				} else {
					m_psf[i][j][k]/=psf_sum;
				}
			}
		}
	}

	delete psf;

}

void PSFGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	unsigned nJitterPerExposure = data->getTimeConfiguration().getJitterSamplingsPerExposure();

	//Create instances of the primary image, together with images to store the PSF positions
	//corresponding to the first and last jitter positions of the exposure,
	//needed to generate the upward and downward smear trails in FrameTransferSmearer
	Image * image = data->getImages().back();
	Image * imageToSmearUp = nullptr;
	Image * imageToSmearDown = nullptr;
	if (data->doFrameTransferSmearing()) {
		imageToSmearUp = new Image(image->getXDim(),Image::kYDim,image->getXOffset(),0);
		imageToSmearDown = new Image(image->getXDim(),Image::kYDim,image->getXOffset(),0);
	}

	//For breathing case, interpolate in temperature between the hot1 and hot2 thermal cases
	double ** temperatureInterpolatedPsf = nullptr;
	double ** temperatureInterpolatedOversampledPsf = nullptr;

	if (m_thermalMap==breathing) {

		//get the fraction between the hot1 (min temperature) and hot2 (max temperature)
		//thermal cases corresponding to the current temperature
		double temperature = data->getSatelliteData(timeStep)->getTelescopeTemperature();
		pair<double,double> breathingTemperatureRange = data->getBreathingTemperatureRange();
		if (breathingTemperatureRange.first==0.) {
			throw runtime_error("Error in PSFGenerator::process: breathingTemperatureRange undefined. Requires OrbitSimulator to be run.");
		}
		double temperatureFrac = (temperature - breathingTemperatureRange.first) /
				(breathingTemperatureRange.second - breathingTemperatureRange.first);
		//cout << breathingTemperatureRange.first-273.15 << " " << breathingTemperatureRange.second-273.15 << " " << temperatureFrac << endl;

		//Generate temperature interpolated PSF for the oversampled case
		if (m_targetStarPositioning==Oversampling || m_backgroundStarPositioning==Oversampling) {
			temperatureInterpolatedOversampledPsf = new double*[kOversampledPsfSize];
			for (int i = 0; i < kOversampledPsfSize; ++i) {
				temperatureInterpolatedOversampledPsf[i] = new double[kOversampledPsfSize];
			}
			for (int j = 0; j<kOversampledPsfSize; j++) {
				for (int i = 0; i<kOversampledPsfSize; i++) {
					temperatureInterpolatedOversampledPsf[i][j] = (1.-temperatureFrac)*m_oversampledPsf[i][j][0] +
																       temperatureFrac*m_oversampledPsf[i][j][1];
				}
			}
		}

		//Generate temperature interpolated PSF for the non-oversampled case
		if (!(m_targetStarPositioning==Oversampling) || !(m_backgroundStarPositioning==Oversampling)) {
			temperatureInterpolatedPsf = new double*[kPsfSize];
			for (int i = 0; i < kPsfSize; ++i) {
				temperatureInterpolatedPsf[i] = new double[kPsfSize];
			}
			for (int j = 0; j<kPsfSize; j++) {
				for (int i = 0; i<kPsfSize; i++) {
					temperatureInterpolatedPsf[i][j] = (1.-temperatureFrac)*m_psf[i][j][0] +
														    temperatureFrac*m_psf[i][j][1];
				}
			}
		}

	}

	//Get the image dimensions and offsets w.r.t. CCD edges
	int xOffset = image->getXOffset();
	int yOffset = image->getYOffset();
	int xDim = image->getXDim();
	int yDim = image->getYDim();
	//If end of life charge transfer is to be simulated, include also stars below the sub-array
	if (data->doChargeTransferEoL()) {
		yDim += yOffset;
		yOffset = 0.;
	}
	int xOffset_smear = xOffset;
	int yOffset_smear = yOffset;
	int xDim_smear = xDim;
	int yDim_smear = yDim;
	if (data->doFrameTransferSmearing()) {
		xOffset_smear = imageToSmearUp->getXOffset();
		yOffset_smear = imageToSmearUp->getYOffset();
		xDim_smear = imageToSmearUp->getXDim();
		yDim_smear = imageToSmearUp->getYDim();
	}

	//clock_t startTime = clock();
	//Loop over all stars in the field of view
	for (unsigned istar=0; istar<data->getFieldOfView()->getStars().size(); istar++) {

		if (fullFrame) cout << "processing star " << istar << endl;

		const Star * star = data->getFieldOfView()->getStars()[istar];
		StarData * starTimeData = star->getTimeSeries()[timeStep];

		//Get the star position on the focal plane
		if ((star->getFocalPlaneX().size()!=nJitterPerExposure && star->getFocalPlaneX().size()!=1) ||
			(star->getFocalPlaneY().size() != star->getFocalPlaneX().size())) {
			throw runtime_error("Error in PSFGenerator::process: size of focal plane position vector for star "
					            +to_string(istar)+" does not match the exposure duration");
		}
		double starPositionX = star->getFocalPlaneX()[0];
		double starPositionY = star->getFocalPlaneY()[0];
		int xpixel = static_cast<int>(floor(starPositionX));
		int ypixel = static_cast<int>(floor(starPositionY));
		//if (timeStep==0) cout << xpixel << " " << ypixel << endl;

		//Require that the star position is within a circular region enclosing the sub-array plus a margin
		//or within the vertical column above and below the sub-array (for frame transfer smearing)
		double deltaX = starPositionX - (xOffset+xDim/2);
		double deltaY = starPositionY - (yOffset+yDim/2);
		double radius = sqrt(deltaX*deltaX + deltaY*deltaY);
		double radius_max = sqrt(xDim*xDim + yDim*yDim)/2. + kPsfMargin;
		if (radius > radius_max && (xpixel<xOffset_smear || xpixel>=xOffset_smear+xDim_smear ||
									ypixel<yOffset_smear || ypixel>=yOffset_smear+yDim_smear)) {
			vector<double> empty;
			image->getTruthData()->addPSF(empty,empty,0.); //So that the PSF truth column vector is the same length for all images
			continue;
		}
		//if (istar==0) cout << "xpixel = " << xpixel << ", ypixel = " << ypixel << endl;

		//Get the the number of PSF positions over the exposure (due to FOV rotation and/or jitter)
		unsigned nJitter = m_oversampleJitter ? star->getFocalPlaneX().size() : 1;
		if (nJitter != nJitterPerExposure && nJitter != 1) {
			throw runtime_error("Error in PSFGenerator::process: FOV rotation vector size must equal exposure time or 1");
		}

		//Generate a separate image for each star
		Image * starImage;
		if (fullFrame) {
			//Full frame image
			starImage = new Image(Image::kXDim,Image::kYDim,0,0);
		} else {
			//Sub-frame image
			starImage = new Image(xDim,yDim,xOffset,yOffset);
		}
		Image * starImageToSmearUp = nullptr;
		Image * starImageToSmearDown = nullptr;
		if (data->doFrameTransferSmearing()) {
			starImageToSmearUp = new Image(starImage->getXDim(),Image::kYDim,starImage->getXOffset(),0);
			starImageToSmearDown = new Image(starImage->getXDim(),Image::kYDim,starImage->getXOffset(),0);
		}

		starPositionX = 0.;
		starPositionY = 0.;
		vector<double> starPositionsX,starPositionsY;

		//Loop over jitter within the exposure and increment the image of the star for each jitter position
		for (unsigned ijitter=0; ijitter<nJitter; ijitter++) {

			//Get the PSF centre position taking into account jitter and FOV rotation
			//(Note: it is already required that the ijitter=0 position is within the sub-array plus margin)
			double jitteredPositionX = star->getFocalPlaneX()[ijitter];
			double jitteredPositionY = star->getFocalPlaneY()[ijitter];
			//double jitteredPositionX = ijitter==0 ? 512-20 :512+20; //fixed jitter offset for testing frame transfer
			//double jitteredPositionY = ijitter==0 ? 512-20 :512+20; //fixed jitter offset for testing frame transfer
			starPositionX += jitteredPositionX;
			starPositionY += jitteredPositionY;
			starPositionsX.push_back(jitteredPositionX);
			starPositionsY.push_back(jitteredPositionY);

			//Use spatial oversampling of PSF to perform sub-pixel positioning
			if ((istar==0 && m_targetStarPositioning==Oversampling) || (istar!=0 && m_backgroundStarPositioning==Oversampling)) {

				//Get the PSF positioning on a grid at 10 times the pixel grid resolution
				//static_cast<int> will truncate the floating point position.
				xpixel = static_cast<int>(jitteredPositionX*10.);
				ypixel = static_cast<int>(jitteredPositionY*10.);
				//Loop over the pixels of the oversampled PSF
				for (int i=xpixel-(kOversampledPsfSize/2); i<xpixel+(kOversampledPsfSize/2); i++) {
					for (int j=ypixel-(kOversampledPsfSize/2); j<ypixel+(kOversampledPsfSize/2); j++) {

						//Get the image pixel corresponding to the current oversampled PSF pixel
						int ix_image = floor(static_cast<double>(i)/10.);
						int iy_image = floor(static_cast<double>(j)/10.);

						if ((ix_image>=xOffset && ix_image<(xOffset+xDim) &&
								iy_image>=yOffset && iy_image<(yOffset+yDim)) ||
								((ijitter==0 || ijitter==nJitter-1) &&
										ix_image>=xOffset_smear && ix_image<(xOffset_smear+xDim_smear) &&
										iy_image>=yOffset_smear && iy_image<(yOffset_smear+yDim_smear))) {

							//Get the PSF value for the pixel
							double pixel_value;
							if (m_thermalMap==breathing) {
								pixel_value = temperatureInterpolatedOversampledPsf[(kOversampledPsfSize/2)+(i-xpixel)][(kOversampledPsfSize/2)+(j-ypixel)];
							} else {
								pixel_value = m_oversampledPsf[(kOversampledPsfSize/2)+(i-xpixel)][(kOversampledPsfSize/2)+(j-ypixel)][0];
							}

							//Increment the pixel value for the current jitter
							if (ix_image>=xOffset && ix_image<(xOffset+xDim) && iy_image>=yOffset && iy_image<(yOffset+yDim)) {
								starImage->incrementPixelValue(ix_image,iy_image,pixel_value/double(nJitter));
							}


							if (data->doFrameTransferSmearing()) {
								//Separately retain the image corresponding to the first and last jitter postions
								//of the exposure for the frame transfer smearing
								if (ijitter==0) {
									starImageToSmearDown->incrementPixelValue(ix_image,iy_image,pixel_value);
								}
								if (ijitter==nJitter-1) {
									starImageToSmearUp->incrementPixelValue(ix_image,iy_image,pixel_value);
								}
							}

						}
					}
				}

			//Interpolation of PSF to sub-pixel position
			} else if ((istar==0 && m_targetStarPositioning==Interpolation) || (istar!=0 && m_backgroundStarPositioning==Interpolation)) {

				//Interpolate the PSF to an array shifted w.r.t. the pixel grid according to the jitter
				double ** interpolatedPSF = new double*[kPsfSize];
				for (int i = 0; i < kPsfSize; ++i) interpolatedPSF[i] = new double[kPsfSize];
				interpolatePSF(jitteredPositionX,jitteredPositionY,temperatureInterpolatedPsf,interpolatedPSF);

				//static_cast<int> will truncate the floating point position
				xpixel = static_cast<int>(jitteredPositionX);
				ypixel = static_cast<int>(jitteredPositionY);
				//cout << ijitter << " " << jitteredPositionX-xOffset << " " << jitteredPositionY-xOffset << " " << xpixel-xOffset << " " << ypixel-xOffset << endl;

				//Add the interpolated PSF onto the image
				for (int i=xpixel-(kPsfSize/2)+1; i<xpixel+(kPsfSize/2)-1; i++) {
					for (int j=ypixel-(kPsfSize/2)+1; j<ypixel+(kPsfSize/2)-1; j++) {
						if ((i>=xOffset && i<(xOffset+xDim) && j>=yOffset && j<(yOffset+yDim)) ||
								((ijitter==0 || ijitter==nJitter-1) &&
										i>=xOffset_smear && i<(xOffset_smear+xDim_smear) &&
										j>=yOffset_smear && j<(yOffset_smear+yDim_smear))) {

							//Map the current pixel to the relevant pixel of the interpolated PSF
							int xpixel_psf = i + (kPsfSize/2) - xpixel;
							int ypixel_psf = j + (kPsfSize/2) - ypixel;

							//Get the PSF value for the pixel
							double pixel_value = interpolatedPSF[xpixel_psf][ypixel_psf];

							//Increment the pixel value for the current jitter
							if (i>=xOffset && i<(xOffset+xDim) && j>=yOffset && j<(yOffset+yDim)) {
								//if (istar==0 && j==ypixel) cout << i << " " << pixel_value << endl;
								starImage->incrementPixelValue(i,j,pixel_value/double(nJitter));
							}

							if (data->doFrameTransferSmearing()) {
								//Separately retain the image corresponding to the first and last jitter postions
								//of the exposure for the frame transfer smearing
								if (ijitter==0) {
									starImageToSmearDown->incrementPixelValue(i,j,pixel_value);
								}
								if (ijitter==nJitter-1) {
									starImageToSmearUp->incrementPixelValue(i,j,pixel_value);
								}
							}

						}
					}
				}

				for(int i = 0; i < kPsfSize; ++i) delete [] interpolatedPSF[i];
				delete [] interpolatedPSF;

			//Position the PSF to the nearest pixel grid position
			} else {

				//static_cast<int> will truncate the floating point position
				xpixel = static_cast<int>(jitteredPositionX);
				ypixel = static_cast<int>(jitteredPositionY);
				//cout << ijitter << " " << jitteredPositionX-xOffset << " " << jitteredPositionY-xOffset << " " << xpixel-xOffset << " " << ypixel-xOffset << endl;

				//Add the PSF onto the image
				for (int i=xpixel-(kPsfSize/2); i<xpixel+(kPsfSize/2); i++) {
					for (int j=ypixel-(kPsfSize/2); j<ypixel+(kPsfSize/2); j++) {
						if ((i>=xOffset && i<(xOffset+xDim) && j>=yOffset && j<(yOffset+yDim)) ||
								((ijitter==0 || ijitter==nJitter-1) &&
										i>=xOffset_smear && i<(xOffset_smear+xDim_smear) &&
										j>=yOffset_smear && j<(yOffset_smear+yDim_smear))) {

							//Get the PSF value for the pixel
							//if (timeStep==0) cout << (kPsfSize/2)+(i-xpixel) << " " << (kPsfSize/2)+(j-ypixel) << " " << m_psf[(kPsfSize/2)+(i-xpixel)][(kPsfSize/2)+(j-ypixel)] << endl;
							double pixel_value;
							if (m_thermalMap==breathing) {
								pixel_value = temperatureInterpolatedPsf[(kPsfSize/2)+(i-xpixel)][(kPsfSize/2)+(j-ypixel)];
							} else {
								pixel_value = m_psf[(kPsfSize/2)+(i-xpixel)][(kPsfSize/2)+(j-ypixel)][0];
							}
							//if (i==xpixel && j==ypixel) cout << iTemp << " " << pixel_value << endl;

							//Increment the pixel value for the current jitter
							if (i>=xOffset && i<(xOffset+xDim) && j>=yOffset && j<(yOffset+yDim)) {
								starImage->incrementPixelValue(i,j,pixel_value/double(nJitter));
							}

							if (data->doFrameTransferSmearing()) {
								//Separately retain the images corresponding to the first and last jitter postions
								//of the exposure for the frame transfer smearing
								if (ijitter==0) {
									starImageToSmearDown->incrementPixelValue(i,j,pixel_value);
								}
								if (ijitter==nJitter-1) {
									starImageToSmearUp->incrementPixelValue(i,j,pixel_value);
								}
							}

						}
					}
				}
			}

		}

		starPositionX /= static_cast<double>(nJitter);
		starPositionY /= static_cast<double>(nJitter);

		//Calculate the mean flux within the ring of pixels whose centres lie at radii between 95 and 96 pixels
		//from the mean centre of the jittered PSF (use 95-96 instead of 99-100 to avoid edge effects due to jitter)
		vector<double> ring95pixels;
		int nInRing = 0;
		//static_cast<int> will truncate the floating point position
		xpixel = static_cast<int>(starPositionX);
		ypixel = static_cast<int>(starPositionY);
		for (int i=xpixel-(kPsfSize/2); i<xpixel+(kPsfSize/2); i++) {
			for (int j=ypixel-(kPsfSize/2); j<ypixel+(kPsfSize/2); j++) {
				if (i>=-Image::kLeftMargin && i<Image::kXTotal-Image::kLeftMargin && j>=0 && j<Image::kYTotal) {
					double deltaX = static_cast<double>(i)+0.5 - starPositionX; //add 0.5 to shift to pixel centre
					double deltaY = static_cast<double>(j)+0.5 - starPositionY; //add 0.5 to shift to pixel centre
					double radius = sqrt(deltaX*deltaX + deltaY*deltaY);
					if (radius>=95. && radius<96.) {
						ring95pixels.push_back(starImage->getPixelValue(i,j));
						nInRing += 1;
					}
				}
			}
		}

		//Use the median value for pixels within the ring to define the scale for the 1/r^3 function that will
		//define the tails of the PSF beyond the edges of the 200x200 pixels for which the PSF is measured
		double rCubedConstant = pow(95.5,3)*CommonTools::median(ring95pixels);

		//Obtain the flux of the star
		double flux = m_convertFluxToPhotons ? star->getMeanPhotonFlux() : star->getMeanFlux();
		double fluxFactor = starTimeData->getCombinedFluxFactor();
		double exposureTimeAsDouble = data->getTimeConfiguration().getExposureTimeAsDouble();
		flux = flux * fluxFactor * kTelescopeArea * exposureTimeAsDouble;

		//Subtract the halo flux from the flux used to fill the PSF so that
		//the overall PSF+halo normalization will be correct after the halo is added
		//NOTE: if HaloGenerator is run, it must be run before PSFGenerator otherwise star->getHaloFlux() will return 0
		//double flux0 = flux;
		flux -= star->getHaloFlux();
		//cout << istar << " " << flux0 << " " << star->getHaloFlux() << " " << flux << endl;

		//cout << starTimeData->getVariationFluxFactor() << " " << flux << endl;

		//Add the star to the image:
		for (int ix=xOffset_smear; ix<xOffset_smear+xDim_smear; ix++) {
			for (int iy=yOffset_smear; iy<yOffset_smear+yDim_smear; iy++) {
				//Set the pixel value by multiplying the PSF intensity for the pixel by the flux
				if (ix>=xOffset && ix<(xOffset+xDim) && iy>=yOffset && iy<(yOffset+yDim)) {

					//Calculate the distance from the mean jittered position of the star to the centre of the pixel
					double deltaX = static_cast<double>(ix)+0.5 - starPositionX; //add 0.5 to shift to pixel centre
					double deltaY = static_cast<double>(iy)+0.5 - starPositionY; //add 0.5 to shift to pixel centre
					double radius = sqrt(deltaX*deltaX + deltaY*deltaY);

					if (m_psfFilename!="FLAT" && radius>=96) {
						//if (fabs(deltaX)>=96 || fabs(deltaY)>=96) { //for square rather than circular PSF truncation
						//If the radius from the centre of the PSF is more than 96 pixels (not 100 in order to avoid
						//edge effects due to jitter), use a 1/r^3 dependence to model extended tails
						image->incrementPixelValue(ix,iy,rCubedConstant*flux/pow(radius,3));
					} else {
						image->incrementPixelValue(ix,iy,starImage->getPixelValue(ix,iy)*flux);
					}
				}
				if (data->doFrameTransferSmearing()) {
					//For the frame transfer smearing, set the flux to correspond to a 1s exposure
					//(to be rescaled to the frame transfer clock period in FrameTransferSmearer)
					imageToSmearUp->incrementPixelValue(ix,iy,starImageToSmearUp->getPixelValue(ix,iy)*flux/exposureTimeAsDouble);
					imageToSmearDown->incrementPixelValue(ix,iy,starImageToSmearDown->getPixelValue(ix,iy)*flux/exposureTimeAsDouble);
				}
			}
		}

		delete starImage;
		if (data->doFrameTransferSmearing()) {
			delete starImageToSmearUp;
			delete starImageToSmearDown;
		}

		//Set the truth information for the PSF
		//cout << starPositionX-512 << " " << starPositionY-512 << endl;
		image->getTruthData()->addPSF(starPositionsX,starPositionsY,flux);
	}
	//cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC<< " seconds" << endl;

	if (temperatureInterpolatedPsf != nullptr) {
		for(int i = 0; i < kPsfSize; ++i) {
			delete [] temperatureInterpolatedPsf[i];
		}
		delete [] temperatureInterpolatedPsf;
	}

	if (temperatureInterpolatedOversampledPsf != nullptr) {
		for(int i = 0; i < kOversampledPsfSize; ++i) {
			delete [] temperatureInterpolatedOversampledPsf[i];
		}
		delete [] temperatureInterpolatedOversampledPsf;
	}

	if (data->doFrameTransferSmearing()) {
		data->addImageToSmearUp(imageToSmearUp);
		data->addImageToSmearDown(imageToSmearDown);
	}

}

void PSFGenerator::interpolatePSF(double jitteredPositionX, double jitteredPositionY,
								  double ** temperatureInterpolatedPsf, double ** interpolatedPSF) const {

	//Get the fractional offset in x and y of the PSF centre position w.r.t. the pixel boundaries
	double xfrac = jitteredPositionX-floor(jitteredPositionX);
	double yfrac = jitteredPositionY-floor(jitteredPositionY);

	//Interpolate the PSF to an array shifted w.r.t. the pixel grid by (xfrac,yfrac)
	//and calculate the PSF integral to ensure it is normalized to 1.
	double interpolatedPSF_sum = 0.;
	for (int ix=0; ix<kPsfSize; ++ix) {
		for (int iy=0; iy<kPsfSize; ++iy) {

			//Protect against boundaries: effectively add an extra row and column to the PSF which are copies
			//of the bottom row and leftmost column such that the PSF grid is 1 row+column larger than the target grid.
			//In this way all pixels in the target grid have 4 nearest neighbours in the PSF grid
			int ixminus1 = ix==0 ? 0 : ix-1;
			int iyminus1 = iy==0 ? 0 : iy-1;

			//Obtain the intensity of the four PSF pixels nearest to the current pixel for the interpolated PSF
			double q11,q12,q21,q22;
			if (m_thermalMap==breathing) {
				q11 = temperatureInterpolatedPsf[ixminus1][iyminus1];
				q12 = temperatureInterpolatedPsf[ixminus1][iy];
				q21 = temperatureInterpolatedPsf[ix][iyminus1];
				q22 = temperatureInterpolatedPsf[ix][iy];
			} else {
				q11 = m_psf[ixminus1][iyminus1][0];
				q12 = m_psf[ixminus1][iy][0];
				q21 = m_psf[ix][iyminus1][0];
				q22 = m_psf[ix][iy][0];
			}

			//Interpolate between the four PSF pixels to obtain the image pixel intensity
			double pixel_value = bilinearInterpolation(q11,q12,q21,q22,xfrac,yfrac);
			interpolatedPSF[ix][iy] = pixel_value;
			interpolatedPSF_sum += pixel_value;
			//if ((i==48+xOffset||i==49+xOffset||i==50+xOffset||i==51+xOffset) && (j==48+xOffset||j==49+xOffset||j==50+xOffset||j==51+xOffset))
			//	 cout << i-xOffset << " " << j-xOffset << " " << ix << " " << iy << " " << q11 << " " << q21 << " " << q12 << " " << q22 << " " << pixel_value << endl;

		}
	}

	//Normalization w.r.t. uninterpolated PSF is 0.99999974 with +/-1 variation in the last decimal place (0.01ppm)
	//cout << setprecision(15) << interpolatedPSF_sum << endl;

	//Normalize the PSF integral to 1
	for (int ix=0; ix<kPsfSize; ++ix) {
		for (int iy=0; iy<kPsfSize; ++iy) {
			interpolatedPSF[ix][iy] /= interpolatedPSF_sum;
		}
	}

}

double PSFGenerator::bilinearInterpolation(double q11, double q12, double q21, double q22, double xfrac, double yfrac) const {

	double r1 = xfrac*q11 + (1.-xfrac)*q21;
	double r2 = xfrac*q12 + (1.-xfrac)*q22;
	return yfrac*r1 + (1.-yfrac)*r2;

}
