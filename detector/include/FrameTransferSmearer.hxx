/*
 * FrameTransferSmearer.hxx
 *
 *  Created on: 22 Jul 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module generates smear trails due to frame transfer for each
 *  	   of the PSFs generated in PSFGenerator.
 *
 * The trails are generated based on the PSF positions at the start and end
 * of the exposure (i.e. first and last jitter values of the exposure).
 */

#ifndef _FRAME_TRANSFER_SMEARER_HXX_
#define _FRAME_TRANSFER_SMEARER_HXX_

#include "simulator/include/Module.hxx"

class FrameTransferSmearer: public Module {
public:
	FrameTransferSmearer() : Module("FrameTransferSmearer",timeLoop) {};
	virtual ~FrameTransferSmearer() {};

	void initialize(const ModuleParams & params);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:
	double m_transferClockPeriod; ///< Clock period for frame transfer in microseconds

};

#endif /* FRAMETRANSFERSMEARER_HXX_ */
