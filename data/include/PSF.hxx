/*
 * PSF.hxx
 *
 *  Created on: 4 Dec 2014
 *      Author: futyand
 */

#ifndef PSF_HXX_
#define PSF_HXX_

#include <vector>
using namespace std;


/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief Class to hold PSF truth information: the flux and list of
 *  	   positions during the exposure
 */

class PSF {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] x  vector of x positions of PSF (one value for each second of the exposure) relative to left edge of exposed part of CDD
	 *  @param [in] y  vector of y positions of PSF (one value for each second of the exposure) relative to left edge of exposed part of CDD
	 *  @param [in] val  Incident flux value for PSF
	 */
	PSF(vector<double> x, vector<double> y, double val);
	virtual ~PSF() {};

	/** *************************************************************************
	 *  @brief Addition operator to add two PSFs, used to combine PSFs with
	 *  	   differing positions due to jitter
	 *
	 *  @param [in] psf  PSF to be added to the current PSF
	 */
	PSF & operator +=(const PSF & psf);

	double getFlux() const {return m_flux;} ///< @brief Returns the flux of the PSF
	double getMeanXPosition() const; ///< @brief Returns the mean x position of the PSF, averaged over the m_xPositions
	double getMeanYPosition() const; ///< @brief Returns the mean y position of the PSF, averaged over the m_yPositions
	vector<double> getXPositions() const {return m_xPositions;} ///< @brief Returns the vector of x positions of the PSF
	vector<double> getYPositions() const {return m_yPositions;} ///< @brief Returns the vector of y positions of the PSF

private:
	vector<double> m_xPositions; ///< Vector of x positions for the PSF at the set of jitter positions corresponding to the image exposure
	vector<double> m_yPositions; ///< Vector of y positions for the PSF at the set of jitter positions corresponding to the image exposure
	double m_flux; ///< Incident flux for the PSF
};

inline PSF operator +(PSF psf1, const PSF & psf2) {
	psf1+=psf2;
	return psf1;
}

#endif /* PSF_HXX_ */
