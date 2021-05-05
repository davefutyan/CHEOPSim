/*
 * StellarVariationFluxModulator.hxx
 *
 *  Created on: 22 Feb 2016
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to define the modulation in the mean flux
///        of the target star resulting from stellar variations
///
////////////////////////////////////////////////////////////////////////

#ifndef SOURCE_INCLUDE_STELLARVARIATIONFLUXMODULATOR_HXX_
#define SOURCE_INCLUDE_STELLARVARIATIONFLUXMODULATOR_HXX_

#include "boost/lexical_cast.hpp"

#include "simulator/include/Module.hxx"

class StellarVariationFluxModulator: public Module {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] starIndex  Index of the star (0 for target star)
	 */
	StellarVariationFluxModulator(unsigned starIndex) : Module("StellarVariationFluxModulator_star"+boost::lexical_cast<string>(starIndex),timeLoop),
														m_starIndex(starIndex), m_rotationPeriod(0.) {};
	virtual ~StellarVariationFluxModulator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	unsigned m_starIndex; ///< Index of Star (0 for target star)
	unsigned m_seed; ///< Seed for random number generation
	double m_rotationPeriod; ///< Rotation period of the star
	vector<double> m_time; ///< Times (seconds) since the start of the simulation corresponding to the values in m_variation
	vector<double> m_variation; ///< Stellar variation flux factor time series for times corresponding to m_time
};

#endif /* SOURCE_INCLUDE_STELLARVARIATIONFLUXMODULATOR_HXX_ */
