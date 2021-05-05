/*
 * CosmicRayGenerator.hxx
 *
 *  Created on: 12 May 2014
 *      Author: futyand
 */

#ifndef _COSMIC_RAY_GENERATOR_HXX_
#define _COSMIC_RAY_GENERATOR_HXX_

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/poisson_distribution.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/random/uniform_real.hpp"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#include "TH1F.h"
//#include "TFile.h"

#include "simulator/include/Module.hxx"

typedef boost::variate_generator<boost::mt19937,boost::uniform_real<double> > RANDOM_UNIFORM_REAL;
typedef boost::variate_generator<boost::mt19937,boost::poisson_distribution<int,double> > RANDOM_POISSON;

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module generates cosmic rays, with random incident positions and
 *  	   directions random directions and with an energy distribution
 *  	   corresponding to HST observations
 */

class CosmicRayGenerator: public Module {
public:

	static constexpr double kPixelDepth = 15; ///< Thickness of CCD Si substrate in microns
	static constexpr double kPixelWidth = 13.; ///< Width of CCD pixel in microns

	CosmicRayGenerator() : Module("CosmicRayGenerator",timeLoop) {};
	virtual ~CosmicRayGenerator();

	void initialize(const ModuleParams & params);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/// @brief Parameters defining the path of the cosmic ray
	struct CosmicPath {
		/** *************************************************************************
		 *  @brief Constructor for CosmicPath
		 *
		 *  @param [in] x  x position of impact on the surface of the pixel grid
		 *  @param [in] y  y position of impact on the surface of the pixel grid
		 *  @param [in] polar  Incident polar angle of the cosmic ray
		 *  @param [in] azimuth  Incident azimuthal angle of the cosmic ray
		 */
		CosmicPath(double x, double y, double polar, double azimuth) :
						 m_xstart(x), m_ystart(y), m_polar(polar), m_azimuth(azimuth) {
			m_length = kPixelDepth/kPixelWidth * tan(m_polar);
			m_length3D = m_length/sin(m_polar);
			m_xend = m_xstart + m_length*cos(m_azimuth);
			m_yend = m_ystart + m_length*sin(m_azimuth);

		};
		double m_xstart; ///< x position of entry on the pixel grid
		double m_ystart; ///< y position of entry on the pixel grid
		double m_xend; ///< x position of exit on the pixel grid
		double m_yend; ///< y position of exit on the pixel grid
		double m_polar; ///< Incident polar angle of the cosmic ray
		double m_azimuth; ///< Incident azimuthal angle of the cosmic ray
		double m_length; ///< Length of cosmic trail projected onto the 2D surface of the CCD in units of pixel widths
		double m_length3D; ///< Length of cosmic trail in 3D in units of pixel widths
	};

	/** *************************************************************************
	 *  @brief Propagate the cosmic ray through the pixel grid in 3D. Returns
	 *  	   false if the end of the trail has been passed
	 *
	 *  @param [in] cosmicPath  Parameters defining the path of the cosmic ray
	 *  @param [in] step  Propagation step size in units of pixel widths
	 *  @param [inout] x  Current x position of the CR in the pixel grid
	 *  @param [inout] y  Current y position of the CR in the pixel grid
	 *  @param [inout] z  Current z position of the CR in the pixel grid
	 */
	bool propagateCosmic(const CosmicPath & cosmicPath, double step, double & x, double & y, double & z) const;

	unsigned m_seed; ///< Seed for random number generation
	double m_cosmicsPerMinute; ///< Mean number of cosmic rays per minute on 200x200 pixels
	double m_SAAFluxFactor; ///< Factor by which the cosmic ray flux rate increases within the SAA
	double m_cosmicEnergyScaleFactor; ///< Factor by which to scale the cosmic ray energy distribution

	static constexpr double kPropagationStep = 0.001; ///< Pixel fraction for cosmic ray propagation steps

	RANDOM_UNIFORM_REAL * m_uniformRandomGenerator; ///< Uniform random number generator
	RANDOM_POISSON * m_poissonRandomGenerator; ///< Poisson random number generator for number of cosmics
	RANDOM_POISSON * m_poissonRandomGenerator2; ///< Poisson random number generator for electrons per pixel in cosmic trails
	gsl_rng * m_gslRandomGenerator; ///< GNU Scientific Library random number generator for Landau distribution

//	TFile * hist_file;
//	TH1F * hist_logNElectrons;
};

/// @ brief Cosmic ray pixel
class CosmicRayPixel {
public:
	/** *************************************************************************
	 *  @brief Constructor for cosmic ray pixel
	 *
	 *  @param [in] ix  x position of the pixel
	 *  @param [in] iy  y position of the pixel
	 *  @param [in] val  Number of electrons in the pixel
	 */
	CosmicRayPixel(int ix, int iy, double val) : m_ix(ix), m_iy(iy), m_value(val) {};
	~CosmicRayPixel() {};

	double getValue() const {return m_value;} ///< Returns the number of electrons in the pixel
	int getXIndex() const {return m_ix;} ///< Returns the x position of the pixel
	int getYIndex() const {return m_iy;} ///< Returns the y position of the pixel
	void setValue(double val) {m_value = val;}  ///< Sets the number of electrons in the pixel
	void incrementValue(double val) {m_value += val;} ///< Increases the number of electrons in the pixel by the specified amount
private:
	int m_ix; ///< x position of the pixel
	int m_iy; ///< y position of the pixel
	double m_value; ///< Number of electrons in the pixel
};

#endif /* _COSMIC_RAY_GENERATOR_HXX_ */
