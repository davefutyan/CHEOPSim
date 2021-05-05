/*
 * FullWellSimulator.hxx
 *
 *  Created on: 5 Aug 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module models the saturation of pixels for which the number
 *  	   of accumulated electrons during an exposure exceeds the full well
 *  	   capacity.
 *
 *  For pixels with number of electrons exceeding the full well capacity,
 *  the excess electrons bleed to adjacent pixels in the vertical direction,
 *  assuming an equal probability to bleed in either direction
 */

#ifndef _FULL_WELL_SIMULATOR_HXX_
#define _FULL_WELL_SIMULATOR_HXX_

#include "simulator/include/Module.hxx"

class FullWellSimulator: public Module {
public:
	FullWellSimulator() : Module("FullWellSimulator",timeLoop) {};
	virtual ~FullWellSimulator() {};

	void initialize(const ModuleParams & params);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Simulates bleeding of charges in the vertical direction when the
	 *  	   number of charges exceeds the full well capacity. Starting from
	 *  	   the pixel with maximum charge in the specified column, charges are
	 *  	   repeatedly shifted in the specified direction (up or down).
	 *  	   If the destination pixel is below capacity, charges are added to
	 *  	   that pixel until it is full, with the remainder moved on to the
	 *  	   next pixel, until all the excess charges have been deposited.
	 *
	 *  @param [in] image  Pointer to the image
	 *  @param [in] ix  x position of the column containing at least one
	 *  				saturated pixel
	 *  @param [in] iy  y position of pixel within column ix which has the
	 *  				maximum number of charges
	 *  @param [in] nOverflow  Number of charges by which pixel ix,iy exceeds
	 *  					   full well capacity, divided by 2 (half the charges
	 *  					   will be moved in each direction)
	 *  @param [in] up  Boolean to indicate whether to shift the excess charges
	 *  				up or down
	 */
	void overflow(Image * image, int ix, int iy, double nOverflow, bool up) const;

	/** *************************************************************************
	 *  @brief Calculates the integrated exposure time for which the pixel
	 *  	   in the image with the highest number of electrons reaches full
	 *  	   well saturation, and writes the value to a file. The method
	 *  	   is not used by default, but the call to it from the process method
	 *  	   can be uncommented in order to calculate time to saturation for
	 *  	   different spectral types, to provide input of this information to
	 *  	   the web interface.
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void calculateTimeToSaturation(const Data * data) const;

	int m_fullWellCapacity; ///< Full well capacity for a pixel
};

#endif /* _FULL_WELL_SIMULATOR_HXX_ */
