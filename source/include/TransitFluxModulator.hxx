/*
 * TransitFluxModulator.hxx
 *
 *  Created on: Dec 18, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to define the modulation in the mean flux
///        of the target star resulting from a planetary transit
///
/// Parameter input for ModuleParams:
///		1. firstTransitTimeFraction (double) : Time of midpoint of first transit as a fraction of the simulation duration
///		2. planetRadius (double) : Planet radius as a multiple of the radius of Jupiter/Neptune/Earth
///		3. orbitPeriod (double) : Planet orbit period around the star in hours
///		4. impactParameter (double) : Impact parameter for the transit (smallest distance from planet centre to star centre divided by star radius)
///		5. doLimbDarkening (bool) : Set to true to include limb darkening
///
////////////////////////////////////////////////////////////////////////

#ifndef _TRANSIT_FLUX_MODULATOR_HXX_
#define _TRANSIT_FLUX_MODULATOR_HXX_

#include "boost/lexical_cast.hpp"

#include "simulator/include/Module.hxx"
#include "TransitModel.hxx"

class TransitFluxModulator: public Module {
public:

	static constexpr double kJupiterRadius = 69911.; ///< Radius of Jupiter in km
	static constexpr double kNeptuneRadius = 24622.; ///< Radius of Neptune in km
	static constexpr double kEarthRadius = 6371.; ///< Radius of Earth in km
	static constexpr double kSunRadius = 696342.; ///< Radius of the Sun in km
	static constexpr double kMuSun = 132712440018.; ///< Standard gravitational parameter of the Sun (G*M_sun) in km^3s^-2

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] starIndex  Index of the star (0 for target star)
	 */
	TransitFluxModulator(unsigned starIndex) : Module("TransitFluxModulator_star"+boost::lexical_cast<string>(starIndex),timeLoop),m_starIndex(starIndex) {};
	virtual ~TransitFluxModulator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	unsigned m_starIndex; ///< Index of Star (0 for target star)
	double m_firstTransitTimeFraction; ///< Time of midpoint of first transit as a fraction of the simulation duration
	double m_planetRadius; ///< Radius of the planet in km
	double m_starRadius; ///< Radius of the star in km
	double m_orbitPeriod; ///< Planet orbit period around the star in seconds
	double m_impactParameter; ///< Impact parameter for the transit (smallest distance from planet centre to star centre divided by star radius)
	bool m_doLimbDarkening; ///< Flag to indicate whether or not to model limb darkening
	bool m_doModification; ///< Flag to indicate whether or not to add a modification to the transit curve read from an input file
	double m_semiMajorAxis; ///< Semi major axis of the orbit of the planet
	BJD m_transitTime_BJD; ///< BJD time of first transit

	TransitModel * m_transitModel; ///< Pointer to an instance of TransitModel
	double m_transitModification[148]; ///< Modification to the transit curve read from an input file
};

#endif /* _TRANSIT_FLUX_MODULATOR_HXX_ */
