/*
 * JitterProducer.hxx
 *
 *  Created on: 13 Feb 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module reads in the jitter time series to define the movement
 *  	   of the position of PSFs on the folca plane due to jitter
 *
 *  The module also defines the AOCS and science validity flags
 */

#ifndef _JITTER_PRODUCER_HXX_
#define _JITTER_PRODUCER_HXX_

#include "data/include/SatelliteData.hxx"
#include "simulator/include/Module.hxx"

class JitterProducer: public Module {
public:

	static const unsigned kNJitter = 174872; ///< Max size for input jitter: corresponds to 48 hours at 1s resolution

	JitterProducer() : Module("JitterProducer",begin) {};
	virtual ~JitterProducer() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const {};

private:
	string m_jitterFilename; ///< Filename for jitter time series
	double m_jitterFileOffset; ///< Number of seconds to omit from the start of the uploaded jitter file
	double m_jitterScale; ///< Scale factor to be applied to jitter RMS
	unsigned m_jitterSize; ///< Number of entries in the jitter file
	SatelliteData::APE m_APE[kNJitter]; ///< Array to contain APE values for each jitter file entry
	bool m_validAocs[kNJitter]; ///< Array to contain AOCS validity for each jitter file entry
	bool m_validScience[kNJitter]; ///< Array to contain science validity for each jitter file entry
};

#endif /* _JITTER_PRODUCER_HXX_ */
