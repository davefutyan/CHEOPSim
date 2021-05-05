/*
 * ChargeTransferSimulator.hxx
 *
 *  Created on: 12 Mar 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module simulates the charge transfer efficiency (CTE)
 *
 *  The effect of CTE is assumed to be cumulative and linear, such that the
 *  number of electrons in a given pixel after charge transfer is equal to
 *  the initial number of electrons multiplied by the product of the CTE and
 *  the number of row/column shifts
 */

#ifndef _CHARGE_TRANSFER_SIMULATOR_HXX_
#define _CHARGE_TRANSFER_SIMULATOR_HXX_

#include "simulator/include/Module.hxx"

class ChargeTransferSimulator: public Module {
public:
	ChargeTransferSimulator() : Module("ChargeTransferSimulator",timeLoop) {};
	virtual ~ChargeTransferSimulator() {};

	void initialize(const ModuleParams & params);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Simulates the effect of charge transfer efficiency at beginning of
	 *  	   life, by simulating the shifting of rows and pixels performed
	 *  	   during frame transfer and readout, assuming a fixed fraction of
	 *  	   charge is not transferred during each shift
	 *
	 *  @param [in] image  Pointer to the Image
	 */
	void processBeginOfLife(Image * image) const;

	/** *************************************************************************
	 *  @brief Simulates the effect of charge transfer efficiency at end of
	 *  	   life, by generating exponential tails in the upward direction for
	 *  	   each pixel, containing a fraction of the pixel charge which
	 *  	   scales with the pixel intensity
	 *
	 *  @param [in] image  Pointer to the Image
	 */
	void processEndOfLife(Image * image) const;

	bool m_endOfLife; ///< Boolean to indicate whether or not to model end of life trails
	double m_cte_vertical; ///< Vertical charge transfer efficiency
	double m_cte_horizontal; ///< Horizontal charge transfer efficiency
	double m_cti_trailFraction; ///< Fraction of intensity transferred to CTI tails for an intensity of 10000 electrons
	double m_cti_intensityScaling; ///< Scaling exponent for fraction of intensity transferred to CTI tails as a function of intensity
	double m_exponentialDistribution[Image::kYDim+Image::kNDarkRows+Image::kNOverscanRows]; ///<Exponential distribution for CTI trails
};

#endif /* _CHARGE_TRANSFER_SIMULATOR_HXX_ */
