/*
 * UserFluxModifier.hxx
 *
 *  Created on: 16 Dec 2016
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module allows the user to define a time series for the variation of the incident flux from the specified star
///
////////////////////////////////////////////////////////////////////////

#ifndef SOURCE_INCLUDE_USERFLUXMODIFIER_HXX_
#define SOURCE_INCLUDE_USERFLUXMODIFIER_HXX_

#include "boost/lexical_cast.hpp"

#include "simulator/include/Module.hxx"

class UserFluxModifier: public Module {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] starIndex  Index of the star (0 for target star)
	 */
	UserFluxModifier(unsigned starIndex) : Module("UserFluxModifier_star"+boost::lexical_cast<string>(starIndex),timeLoop),m_starIndex(starIndex) {};
	virtual ~UserFluxModifier() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	string m_fluxFilename; ///< Filename for user defined flux time series
	unsigned m_starIndex; ///< Index of Star (0 for target star)
	vector<double> m_fluxFactor; ///< Time series for the flux value relative to the nominal flux of the star with index m_starIndex
	vector<double> m_time; ///< Time in seconds from the start of the simulation for each entry in m_fluxFactor
};

#endif /* SOURCE_INCLUDE_USERFLUXMODIFIER_HXX_ */
