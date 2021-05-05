/*
 * Image.hxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief 2D image array
///
/// The vector of Images available the simulation
/// (each vector element corresponds to a different time step) can be
/// accessed from the Data class via Data::getImages()
///
////////////////////////////////////////////////////////////////////////

#ifndef _IMAGE_HXX_
#define _IMAGE_HXX_

#include "TruthData.hxx"

class Image {
public:

	static const int kXDim = 1024; ///< Number of pixels in X dimension for full frame
	static const int kYDim = 1024; ///< Number of pixels in Y dimension for full frame
	static const int kNDarkCols = 16; ///< Number of dark reference columns
	static const int kNBlankCols = 8; ///< Number of blank reference columns
	static const int kNOverscanCols = 4; ///< Number of overscan reference columns
	static const int kNDarkRows = 3; ///< Number of dark reference rows
	static const int kNOverscanRows = 6; ///< Number overscan reference rows
	static const int kXTotal = kXDim+2*kNDarkCols+2*kNBlankCols+kNOverscanCols; ///< Total x array size including reference columns
	static const int kYTotal = kYDim+kNDarkRows+kNOverscanRows; ///< Total y array size including reference rows
	static const int kLeftMargin = kNDarkCols+kNBlankCols+kNOverscanCols; ///< Number of reference columns to left of exposed pixels

	/** *************************************************************************
	 *  @brief Constructor from image dimensions and offset
	 *
	 *  @param [in] xDim  Number of pixels for x dimension (default 200)
	 *  @param [in] yDim  Number of pixels for y dimension (default 200)
	 *  @param [in] xOffset  Offset in x direction of the left edge of the image
	 *  					 w.r.t. the left edge of the exposed part of the CCD
	 *  					 (default 412)
	 *  @param [in] yOffset  Offset in y direction of the bottom edge of the image
	 *  					 w.r.t. the bottom edge of the exposed part of the CCD
	 *  					 (default 412)
	 */
	Image(int xDim=200,int yDim=200,int xOffset=412,int yOffset=412);
	virtual ~Image() {delete m_truthData;}

	/** *************************************************************************
	 *  @brief Addition operator: performs addition of two images, including
	 *  	   associated truth data
	 */
	Image & operator +=(const Image & image);

	/** *************************************************************************
	 *  @brief Multiplication operator: Multiplies the value of all pixels in the
	 *  	   image by the specified scale factor
	 *
	 *  @param [in] scaleFactor  Scale factor by which to multiply all pixels in
	 *  						 the image
	 */
	Image & operator *=(double scaleFactor);

	/** *************************************************************************
	 *  @brief Returns the value for the pixel at the specified position
	 *
	 *  @param [in] ix  x position of pixel w.r.t. the left edge of the exposed
	 *  				part of the CCD
	 *  @param [in] iy  y position of pixel w.r.t. the bottom edge of the exposed
	 *  				part of the CCD
	 *  @return the value for the pixel at the specified position
	 */
	double getPixelValue(int ix, int iy) const {return m_image[ix+kLeftMargin][iy];}

	/** *************************************************************************
	 *  @brief Assigns the value for the pixel at the specified position
	 *
	 *  @param [in] ix  x position of pixel w.r.t. the left edge of the exposed
	 *  				part of the CCD
	 *  @param [in] iy  y position of pixel w.r.t. the bottom edge of the exposed
	 *  				part of the CCD
	 *  @param [in] value  Value for the pixel at the specified position
	 */
	void setPixelValue(int ix, int iy, double value) {m_image[ix+kLeftMargin][iy] = value;}

	/** *************************************************************************
	 *  @brief Increments the value for the pixel at the specified position by
	 *  	   the specified amount
	 *
	 *  @param [in] ix  x position of pixel w.r.t. the left edge of the exposed
	 *  				part of the CCD
	 *  @param [in] iy  y position of pixel w.r.t. the bottom edge of the exposed
	 *  				part of the CCD
	 *  @param [in] value  Amount by which to increment the value for the pixel
	 *  				   at the specified position
	 */
	void incrementPixelValue(int ix, int iy, double value) {m_image[ix+kLeftMargin][iy] += value;}

	/** *************************************************************************
	 *  @brief Performs integer rounding for all pixels in the image.
	 *  	   Should be called before image stacking.
	 */
	void roundToIntegers();

	/** *************************************************************************
	 *  @brief Performs ADC saturation according to the specified number of bits.
	 *  	   Should be performed before integer rounding and image stacking.
	 */
	void saturate(unsigned nBits);

	/** *************************************************************************
	 *  @brief Assigns the truth data associated to the image
	 */
	void setTruthData(TruthData * truthData);

	/** *************************************************************************
	 *  @brief Returns the truth data associated to the image
	 */
	TruthData * getTruthData() const {return m_truthData;}

	/** *************************************************************************
	 *  @brief Returns the number of pixels for x dimension
	 */
	int getXDim() const {return m_xDim;}

	/** *************************************************************************
	 *  @brief Returns the number of pixels for y dimension
	 */
	int getYDim() const {return m_yDim;}


	/** *************************************************************************
	 * @brief Returns the offset in the x direction of the left edge of the
	 * 		  image w.r.t. the left edge of the exposed part of the CCD
	 */
	int getXOffset() const {return m_xOffset;}

	/** *************************************************************************
	 * @brief Returns the offset in the y direction of the bottom edge of the
	 * 		  image w.r.t. the bottom edge of the exposed part of the CCD
	 */
	int getYOffset() const {return m_yOffset;}

private:
	double m_image[kXTotal][kYTotal]; ///< Image array
	int m_xDim; ///< Number of pixels in X dimension
	int m_yDim; ///< Number of pixels in Y dimension
	int m_xOffset; ///< Offset of first pixel in X dimension w.r.t. CCD edge
	int m_yOffset; ///< Offset of first pixel in Y dimension w.r.t. CCD edge
	TruthData * m_truthData; ///< Truth information corresponding to the image
};

inline Image operator +(Image image1, const Image & image2) {
	image1+=image2;
	return image1;
}

#endif /* _IMAGE_HXX_ */
