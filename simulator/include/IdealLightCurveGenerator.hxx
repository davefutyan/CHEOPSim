/*
 * IdealLightCurveGenerator.hxx
 *
 *  Created on: 5 Jul 2017
 *      Author: futyand
 */

#ifndef _IDEAL_LIGHT_CURVE_GENERATOR_HXX_
#define _IDEAL_LIGHT_CURVE_GENERATOR_HXX_

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/poisson_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include "simulator/include/Module.hxx"
#include "simulator/include/Photometry.hxx"

typedef boost::variate_generator<boost::mt19937,boost::poisson_distribution<int,double> > RANDOM_POISSON;


/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module generates photon (or shot) noise, by drawing
 *  	   randomly from a Poisson distribution, applies it to a
 *  	   temporary copy of the image and performs photometry
 *  	   before application of downstream detector effects
 */

class IdealLightCurveGenerator: public Module {
public:
	IdealLightCurveGenerator() : Module("IdealLightCurveGenerator",timeLoop) {};
	virtual ~IdealLightCurveGenerator() {delete m_poissonNoiseGenerator;
										 delete m_photometry;}

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
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

	RANDOM_POISSON * m_poissonNoiseGenerator; ///< Poisson random number generator

	Photometry * m_photometry; ///< Pointer to an instance of the Photometry class, used to perform the photometric extraction

};

#endif /* _IDEAL_LIGHT_CURVE_GENERATOR_HXX_ */
