/*
 * GlobalThroughputGenerator.hxx
 *
 *  Created on: 3 Feb 2016
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module applies the global throughput by applying the combined
 *  	   effect of three wavelength dependent quantities: the stellar black
 *  	   body spectrum, the optical throughput and quantum efficiency
 *
 *  For a given wavelength, the module converts the number of incoming photons
 *  whose path arrives at a certain pixel to the number of signal electrons in
 *  that pixel
 */

#ifndef TELESCOPE_INCLUDE_GLOBALTHROUGHPUTGENERATOR_HXX_
#define TELESCOPE_INCLUDE_GLOBALTHROUGHPUTGENERATOR_HXX_

#include "simulator/include/Module.hxx"

class GlobalThroughputGenerator: public Module {
public:
	GlobalThroughputGenerator() : Module("GlobalThroughputGenerator",timeLoop) {};
	virtual ~GlobalThroughputGenerator() {};

	void initialize(const ModuleParams & params) {};
	void doBegin(Data *data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:
	double m_effectiveTemperature; ///< Black body effective temperature of the target star used to weight the throughput*QE wavelength integral
	string m_spectralType; ///< Spectral type of the target star

};

#endif /* TELESCOPE_INCLUDE_GLOBALTHROUGHPUTGENERATOR_HXX_ */
