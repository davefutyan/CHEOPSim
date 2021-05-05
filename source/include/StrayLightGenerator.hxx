/*
 * StrayLightGenerator.hxx
 *
 *  Created on: 20 Feb 2016
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to determine the stray light flux incident on the telescope
///		   as a function of orbit position, for a user selected scenario defined for
///		   a typical pointing direction at a certain time of year
///
////////////////////////////////////////////////////////////////////////

#ifndef SOURCE_INCLUDE_STRAYLIGHTGENERATOR_HXX_
#define SOURCE_INCLUDE_STRAYLIGHTGENERATOR_HXX_

#include "satellite/include/OrbitSimulator.hxx"

#include "simulator/include/Module.hxx"

class StrayLightGenerator: public Module {
public:

	StrayLightGenerator() : Module("StrayLightGenerator",timeLoop) {};
	virtual ~StrayLightGenerator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data *data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Assigns values for the m_timeFromFile and m_strayLightFromFile arrays
	 *
	 *  @param [in] time  value of time corresponding to index
	 *  @param [in] flux  value of stray light flux corresponding to index
	 */
	void setStrayLightValues(double time, double flux);

	/** *************************************************************************
	 *  @brief Returns the stray light value for the specified timestep
	 *
	 *  @param [in] timeStep  time step of the simulation
	 *  @param [in] timeConf  the time configuration of the simulation
	 */
	double getStrayLightFlux(int timeStep, TimeConfiguration timeConf) const;

	string m_strayLightFilename; ///< Filename of the file containing the stray light time series
	vector<double> m_strayLightFromFile; ///< Array of stray light values
	vector<double> m_timeFromFile; ///< Array of times corresponding to m_strayLightFromFile
	double m_integratedThroughput; ///< Wavelength integrated throughput, weighted according to the black body spectrum of the target star
	bool m_strayLightFromVisitConstraints; ///< Flag to indicate if the stray light should be read from the MPS visit constraints rather than from a selected/uploaded time series
	MpsPreVisitconstraints * m_visitConstraints; ///< Pointer to the MPS visit constraints

};

#endif /* SOURCE_INCLUDE_STRAYLIGHTGENERATOR_HXX_ */
