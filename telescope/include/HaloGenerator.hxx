/*
 * HaloGenerator.hxx
 *
 *  Created on: 14 Jul 2016
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief Generates PSF halo due to scattering and ghosts
 *
 *  This module generates a halo due to scattering for each PSF in the FOV,
 *  according to a simplified analytic model (Peterson), and generates
 *  ghost flux uniformly over the CCD according to the position and
 *  flux of each star in the FOV
 */

#include "simulator/include/Module.hxx"

#ifndef TELESCOPE_INCLUDE_HALOGENERATOR_HXX_
#define TELESCOPE_INCLUDE_HALOGENERATOR_HXX_

class HaloGenerator: public Module {
public:
	HaloGenerator() : Module("HaloGenerator",timeLoop), m_oversampleJitter(false), m_doGhosts(false) {};
	virtual ~HaloGenerator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Generates ghost flux for each star uniformly across the CCD. The
	 *  	   total ghost flux for a given star is a fraction of the total PSF
	 *  	   flux, which depends on the angular distance of the star from the
	 *  	   pointing direction
	 *
	 *  @param [in] image  Pointer to the Image
	 *  @param [in] incidentFlux  Incident flux from the star in photons/cm2
	 *  @param [in] starPositionX  X position of the star on the focal plane in
	 *  						   pixel units w.r.t left edge of exposed part
	 *  						   of CCD
	 *  @param [in] starPositionY  Y position of the star on the focal plane in
	 *  						   pixel units w.r.t bottom edge of exposed part
	 *  						   of CCD
	 *  @param [in] doChargeTransferEoL  Flag to indicate whether or not end of
	 *  						   		 life CTI is to be simulated
	 */
	void generateGhostFlux(Image * image, double incidentFlux, double starPositionX, double starPositionY, bool doChargeTransferEoL) const;
	bool m_oversampleJitter; ///< Flag to indicate whether or not to generate the scattering halo once per second rather than once per exposure
	bool m_doGhosts; ///< Flag to indicate whether or not to include flux due to ghosts

};

#endif /* TELESCOPE_INCLUDE_HALOGENERATOR_HXX_ */
