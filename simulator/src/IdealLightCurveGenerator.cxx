/*
 * IdealLightCurveGenerator.cxx
 *
 *  Created on: 5 Jul 2017
 *      Author: futyand
 */

#include "IdealLightCurveGenerator.hxx"

void IdealLightCurveGenerator::initialize(const ModuleParams& params) {

	boost::mt19937 randomNumberEngine(2468);
    boost::poisson_distribution<int,double> poissonDistribution(1.);
    m_poissonNoiseGenerator = new RANDOM_POISSON(randomNumberEngine,poissonDistribution);

}

void IdealLightCurveGenerator::doBegin(Data* data, bool fullFrame) {

	//Initialize the photometric extraction
	m_photometry = new Photometry(data->getSubarrayDimensions(),data->getPhotometryParams());

}

void IdealLightCurveGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	//Make a copy of the current image
	Image * image_copy;
	int yDim = data->getSubarrayDimensions().m_yDim;
	int yOffset = data->getSubarrayDimensions().m_yOffset;
	int xDim = data->getSubarrayDimensions().m_xDim;
	int xOffset = data->getSubarrayDimensions().m_xOffset;
	if (fullFrame) {
		//Full frame image
		image_copy = new Image(Image::kXDim,Image::kYDim,0,0);
	} else {
		//Sub-frame image with default dimensions
		image_copy = new Image(xDim,yDim,xOffset,yOffset);
	}
	*image_copy += *(data->getImages().back());

	int	xmin = image_copy->getXOffset();
	int	xmax = image_copy->getXOffset()+image_copy->getXDim();
	int	ymin = image_copy->getYOffset();
	int	ymax = image_copy->getYOffset()+image_copy->getYDim();

	//Generate noise independently for each pixel
	for (int ix=xmin; ix<xmax; ix++) {
		for (int iy=ymin; iy<ymax; iy++) {
			generatePhotonNoise(image_copy,ix,iy);
		}
	}

	//Determine PSF barycentre using truth information
	pair<double,double> barycentre;
	vector<PSF> PSFs = image_copy->getTruthData()->getPSFs();
	if (PSFs.size()>0) {
		Data::SubarrayDimensions subarray = data->getSubarrayDimensions();
		barycentre = make_pair(PSFs[0].getMeanXPosition()-subarray.m_targetLocationX,
							   PSFs[0].getMeanYPosition()-subarray.m_targetLocationY);
	} else {
		barycentre = make_pair(0.,0.); //default to intended target location if there are no PSFs
	}

	//Extract the flux from the image array
	double flux = m_photometry->extractFlux(image_copy,barycentre);
	data->appendIdealLightCurve(flux);

	delete image_copy;

}

void IdealLightCurveGenerator::generatePhotonNoise(Image* image, int ix, int iy) const {

	double nPhotons = image->getPixelValue(ix,iy);

	//Draw from Poisson distribution with mean nPhotons
	if (nPhotons>0.) {

		m_poissonNoiseGenerator->distribution().param(boost::poisson_distribution<int,double>::param_type(nPhotons));
		int nPhotonsWithNoise = (*m_poissonNoiseGenerator)();
		image->setPixelValue(ix,iy,double(nPhotonsWithNoise));

	} else {

		if (nPhotons<0.) throw runtime_error("ERROR in IdealLightCurveGenerator: Negative flux at entry to IdealLightCurveGenerator for pixel ("+to_string(ix)+","+to_string(iy)+"): "+to_string(nPhotons));
		image->setPixelValue(ix,iy,0.);

	}

}
