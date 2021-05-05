/*
 * PhotonNoiseGenerator.cxx
 *
 *  Created on: 27 Feb 2014
 *      Author: futyand
 */

#include "boost/random/poisson_distribution.hpp"

#include "PhotonNoiseGenerator.hxx"

typedef boost::variate_generator<boost::mt19937,boost::poisson_distribution<int,double> > RANDOM_POISSON;

void PhotonNoiseGenerator::initialize(const ModuleParams& params) {

	//Initialize a uniform random number generator to generate seeds for the Poisson RNG.
	//This allows to reset the Poisson RNG and assign a new seed after each call.
	//This in turn guarantees that the output of the Poisson RNG for a given pixel
	//will always be the same for a given value for that pixel since it cannot depend on
	//previous states of the RNG. This means that the shot noise pattern will be identical
	//for non-hot pixels with and without hot pixels switched on.
	unsigned seed = params.GetAsInt("photonNoiseSeed");
	boost::mt19937 randomNumberEngine(seed);
	boost::uniform_int<int> uniform_random(0,1E9);
	m_uniformRandomGenerator = new RANDOM_UNIFORM(randomNumberEngine,uniform_random);

	//Initialize the Poisson RNG, seeding it by drawing from the uniform RNG
	m_randomNumberEngine = new boost::mt19937((*m_uniformRandomGenerator)());

}

void PhotonNoiseGenerator::doBegin(Data* data, bool fullFrame) {

	data->getIdealLightCurveParams()->m_photonNoise = true;

}

void PhotonNoiseGenerator::process(Data* data, int timeStep, bool fullFrame) const {

	Image* image = data->getImages().back();

	int xmin = image->getXOffset();
	int xmax = image->getXOffset()+image->getXDim();
	int ymin = data->doChargeTransferEoL() ? 0. : image->getYOffset(); //Include pixels below the subarray if end of life CTI is to be simulated
	int ymax = image->getYOffset()+image->getYDim();

	//Generate noise independently for each pixel
    //clock_t time0 = clock();
	for (int ix=xmin; ix<xmax; ix++) {
		for (int iy=ymin; iy<ymax; iy++) {
			generatePhotonNoise(image,ix,iy);
		}
		//Also generate noise also for dark and overscan rows (for transfer smear trail)
		for (int iy=Image::kYDim; iy<Image::kYDim+Image::kNDarkRows+Image::kNOverscanRows; iy++) {
			generatePhotonNoise(image,ix,iy);
		}
	}
    //cout << double( clock()- time0 ) / (double)CLOCKS_PER_SEC << " seconds" << endl;

}

void PhotonNoiseGenerator::generatePhotonNoise(Image* image, int ix, int iy) const {

	double nPhotons = image->getPixelValue(ix,iy);

	//Draw from Poisson distribution with mean nPhotons
	if (nPhotons>0.) {

		//Assign a new seed for the Poisson RNG drawn randomly from the uniform RNG (See comment in the initialize method)
		m_randomNumberEngine->seed((*m_uniformRandomGenerator)());

		//Set up the Poisson RNG, assigning the lambda parameter according to the expected value
		boost::poisson_distribution<int,double> poissonDistribution(nPhotons);
		RANDOM_POISSON poissonNoiseGenerator((*m_randomNumberEngine),poissonDistribution);

		//Replace the pixel value with a value randomly drawn from the Poisson RNG
		int nPhotonsWithNoise = poissonNoiseGenerator();
		image->setPixelValue(ix,iy,double(nPhotonsWithNoise));

	} else {

		if (nPhotons<0.) throw runtime_error("ERROR in PhotonNoiseGenerator: Negative flux at entry to PhotonNoiseGenerator for pixel ("+to_string(ix)+","+to_string(iy)+"): "+to_string(nPhotons));
		image->setPixelValue(ix,iy,0.);

	}

}
