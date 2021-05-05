/*
 * StellarNoiseFluxModulator.hxx
 *
 *  Created on: 3 Feb 2014
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to define the modulation in the mean flux
///        of the target star resulting from stellar noise
///
////////////////////////////////////////////////////////////////////////


#ifndef _STELLAR_NOISE_FLUX_MODULATOR_HXX_
#define _STELLAR_NOISE_FLUX_MODULATOR_HXX_

#include "boost/lexical_cast.hpp"

#include "simulator/include/Module.hxx"

class StellarNoiseFluxModulator: public Module {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] starIndex  Index of the star (0 for target star)
	 */
	StellarNoiseFluxModulator(unsigned starIndex) : Module("StellarNoiseFluxModulator_star"+boost::lexical_cast<string>(starIndex),timeLoop),m_starIndex(starIndex) {};
	virtual ~StellarNoiseFluxModulator() {};

	void initialize(const ModuleParams & params) {};
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:
	unsigned m_starIndex; ///< Index of Star (0 for target star)
	double m_granulation[11518]; ///< Stellar granulation flux factor time series in 15 second intervals
};

#endif /* _STELLAR_NOISE_FLUX_MODULATOR_HXX_ */
