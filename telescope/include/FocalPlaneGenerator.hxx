/*
 * FocalPlaneGenerator.hxx
 *
 *  Created on: 5 Feb 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module projects the positions of stars within the field of
 *  	   view onto the focal plane of the detector using the
 *  	   fits_world_to_pix World Coordinate System routine.
 *
 *  The gnomonic (tangent plane) projection and plate scale of 1 pixel per
 *  arcsecond are used. The projection algorithm uses the roll angle as a
 *  function of time defined in OrbitSimulator.
 */

#ifndef _FOCAL_PLANE_GENERATOR_HXX_
#define _FOCAL_PLANE_GENERATOR_HXX_

#include "REF_APP_PixelScale.hxx"

#include "simulator/include/Module.hxx"

class FocalPlaneGenerator: public Module {
public:
	FocalPlaneGenerator() : Module("FocalPlaneGenerator",timeLoop) {};
	virtual ~FocalPlaneGenerator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data *data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:
	int m_subArrayXDim; ///< Number of pixels in X dimension of image sub-array
	int m_subArrayYDim; ///< Number of pixels in Y dimension of image sub-array
	int m_subArrayXOffset; ///< Offset of first pixel of image sub-array in X dimension w.r.t. CCD edge
	int m_subArrayYOffset; ///< Offset of first pixel of image sub-array in Y dimension w.r.t. CCD edge
	double m_targetLocationX; ///< Location of centre of FOV rotation in X dimension w.r.t. CCD edge
	double m_targetLocationY; ///< Location of centre of FOV rotation in Y dimension w.r.t. CCD edge
	vector<double> m_targetLocationListX; ///< vector of target location x coordinates, with one entry per exposyre, intended for MandC data
	vector<double> m_targetLocationListY; ///< vector of target location y coordinates, with one entry per exposure, intended for MandC data
	double m_plateScale; ///< Plate scale in arcseconds per pixel read from REF_APP_PixelScale reference file
};

#endif /* FOCAL_PLANE_GENERATOR_HXX_ */
