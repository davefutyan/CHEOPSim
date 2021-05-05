/*
 * Photometry.hxx
 *
 *  Created on: 21 Nov 2014
 *      Author: futyand
 */

#ifndef PHOTOMETRY_HXX_
#define PHOTOMETRY_HXX_

#include "data/include/Data.hxx"
#include "data/include/Image.hxx"

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief Class to perform photometric extraction from an Image
 *
 *	Parameters for the photometric extraction are passed to the constructor.
 *	The setFlatField method can optionally be used to define a flat field
 *	correction during photometric extraction. Photometric extraction is
 *	performed using the extractFlux method. Methods are also provided
 *	to evaluate the PSF barycentre and the median background.
 */

class Photometry {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] subarrayDimensions  Subarray dimensions
	 *  @param [in] radius_barycentre  Radius of circle centred on intended target location used to determine PSF barycentre
	 *  @param [in] radius_bkgInner  Inner radius of annulus centred on intended target location used to evaluate the background
	 *  @param [in] radius_bkgOuter  Outer radius of annulus centred on intended target location used to evaluate the background
	 *  @param [in] radius_psf  Radius of circle centred on PSF barycentre used to extract the signal flux
	 *  @param [in] subtractBackground  Flag to indicate whether or not to perform background subtraction
	 */
	Photometry(const Data::SubarrayDimensions & subarrayDimensions,
			   double radius_barycentre, double radius_bkgInner, double radius_bkgOuter,
			   double radius_psf=0, bool subtractBackground=true) :
		m_targetLocationX(subarrayDimensions.m_targetLocationX), m_targetLocationY(subarrayDimensions.m_targetLocationY),
		m_radius_barycentre(radius_barycentre), m_radius_bkgInner(radius_bkgInner),
		m_radius_bkgOuter(radius_bkgOuter), m_radius_psf(radius_psf),
		m_subtractBackground(subtractBackground), m_flatField(nullptr) {};

	/** *************************************************************************
	 *  @brief Constructor from PhotometryParams
	 *
	 *  @param [in] subarrayDimensions  Subarray dimensions
	 *  @param [in] photometry_params  PhotometryParams struct containing the
	 *  							   input parameters for photometric extraction
	 */
	Photometry(const Data::SubarrayDimensions & subarrayDimensions, const Data::PhotometryParams & photometry_params) :
		m_targetLocationX(subarrayDimensions.m_targetLocationX), m_targetLocationY(subarrayDimensions.m_targetLocationY),
		m_radius_barycentre(photometry_params.m_radius_barycentre), m_radius_bkgInner(photometry_params.m_radius_bkgInner),
		m_radius_bkgOuter(photometry_params.m_radius_bkgOuter), m_radius_psf(photometry_params.m_radius_psf),
		m_subtractBackground(photometry_params.m_subtractBackground), m_flatField(nullptr) {};

	virtual ~Photometry() {};

	/** *************************************************************************
	 *  @brief Sets the pointer to an image used to define the flat field correction
	 *
	 *  @param [in] flatField  Pointer to an image defining the flat field correction
	 */
	void setFlatField(Image * flatField) {m_flatField = flatField;}

	/** *************************************************************************
	 *  @brief Performs photometric extraction, optionally with background
	 *  	   subtraction and flat field correction
	 *
	 *  @param [in] image  Pointer to the Image
	 *  @param [in] barycentre  Position of the barycentre used to define the
	 *  					    centre of the aperture
	 *  @return integrated flux within photometric extraction aperture
	 */
	double extractFlux(const Image * image, std::pair<double,double> barycentre) const;

	/** *************************************************************************
	 *  @brief Evaluates the position of the barycentre of the PSF
	 *
	 *  @param [in] image  Pointer to the Image
	 *  @return position of the PSF barycentre
	 */
	std::pair<double,double> getBarycentre(const Image * image) const;

	/** *************************************************************************
	 *  @brief Evaluates the median value of pixel fluxes within the user defined
	 *  	   annulus to be used for background subtraction
	 *
	 *  @param [in] image  Pointer to the Image
	 *  @param [in] xOffset  x offset of annulus centre w.r.t. centre of CCD
	 *  @param [in] yOffset  y offset of annulus centre w.r.t. centre of CCD
	 *  @return median background flux
	 */
	double medianBackground(const Image * image, double xOffset, double yOffset) const;

private:

	double m_targetLocationX; ///<Intended location of target in x direction w.r.t. the bottom edge of the exposed part of the CCD
	double m_targetLocationY; ///<Intended location of target in y direction w.r.t. the bottom edge of the exposed part of the CCD
	double m_radius_barycentre; ///< Radius of circle centred on image centre used to determine PSF barycentre
	double m_radius_bkgInner; ///< Inner radius of annulus centred on image centre used to evaluate the background
	double m_radius_bkgOuter; ///< Inner radius of annulus centred on image centre used to evaluate the background
	double m_radius_psf; ///< Radius of circle centred on PSF barycentre used to extract the signal flux
	bool m_subtractBackground; ///< Flag to indicate whether or not to perform background subtraction

	Image * m_flatField; ///< Pointer to an image used to define the flat field correction
};

#endif /* PHOTOMETRY_HXX_ */
