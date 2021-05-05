/*
 * PhotonNoiseGenerator.hxx
 *
 *  Created on: 27 Feb 2014
 *      Author: futyand
 */

#ifndef _PHOTON_NOISE_GENERATOR_HXX_
#define _PHOTON_NOISE_GENERATOR_HXX_

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/random/variate_generator.hpp"

#include "simulator/include/Module.hxx"

typedef boost::variate_generator<boost::mt19937,boost::uniform_int<int> > RANDOM_UNIFORM;


/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module generates photon (or shot) noise, by drawing
 *  	   randomly from a Poisson distribution
 */

class PhotonNoiseGenerator: public Module {
public:
	PhotonNoiseGenerator() : Module("PhotonNoiseGenerator",timeLoop) {};
	virtual ~PhotonNoiseGenerator() {delete m_uniformRandomGenerator;
									 delete m_randomNumberEngine;}

	void initialize(const ModuleParams & params);
	void doBegin(Data *data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Generates photon noise for the specified pixel of the specified image,
	 *  	   using Poisson random number generation
	 *
	 *  @param [in] image  Pointer to the Image
	 *  @param [in] ix  x index of pixel, relative to the left edge of the sub-array
	 *  @param [in] iy  y index of pixel, relative to the bottom edge of the sub-array
	 */
	void generatePhotonNoise(Image * image, int ix, int iy) const;

	RANDOM_UNIFORM * m_uniformRandomGenerator; ///< Uniform random number generator used to generate seeds for m_poissonNoiseGenerator
	boost::mt19937 * m_randomNumberEngine; ///< Random number engine for m_poissonNoiseGenerator

};

#endif /* _PHOTON_NOISE_GENERATOR_HXX_ */
