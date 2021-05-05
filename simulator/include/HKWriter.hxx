/*
 * HKWriter.hxx
 *
 *  Created on: 6 Sep 2017
 *      Author: futyand
 */

#ifndef _HKWRITER_HXX_
#define _HKWRITER_HXX_

#include "simulator/include/Module.hxx"

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module generates houskeeping data over the duration of the
 *  	   simulation, generated using Hk2RawProcessing
 */

class HKWriter: public Module {
public:
	HKWriter() : Module("HKWriter",end) {};
	virtual ~HKWriter() {};

	void initialize(const ModuleParams & params) {m_sdsOffset = params.GetAsInt("sdsCounterOffset");};
	void process(Data * data, int timeStep, bool fullFrame=false) const {};
	void doEnd(Data *data) const;

private:
	unsigned m_sdsOffset; ///< SDS counter initial offset, intended for Dark Frame M&C data, for which SDS counter increments through consecutive visits

};

#endif /* _HKWRITER_HXX_ */
