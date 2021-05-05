/*
 * StarData.hxx
 *
 *  Created on: Dec 19, 2013
 *      Author: futyand
 */

#ifndef _STAR_DATA_HXX_
#define _STAR_DATA_HXX_

#include <vector>
using namespace std;

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Container class for time dependent star data.
///
/// A vector of pointers to StarData (one element for each time step)
/// is provided in the Star class.
///
/// Each instance of StarData (corresponding to a given time step)
/// provides for that time step the multiplicative factors to be applied to the flux of a Star
/// (relative to Star::meanFlux() or Star::meanFlux(double wavelength)
///
////////////////////////////////////////////////////////////////////////

class StarData {
public:
	StarData();
	virtual ~StarData() {};

	/// @brief Returns the multiplicative factor to be applied to the flux of a Star resulting from a planetary transit
	double getTransitFluxFactor() const {return m_transitFluxFactor;}

	/// @brief  Returns the multiplicative factor to be applied to the flux of a Star resulting from stellar variation
	double getVariationFluxFactor() const {return m_variationFluxFactor;}

	/// @brief  Returns the multiplicative factor to be applied to the flux of a Star resulting from stellar noise
	double getNoiseFluxFactor() const {return m_noiseFluxFactor;}

	/// @brief  Returns the user defined multiplicative factor to be applied to the flux of a Star
	double getUserFluxFactor() const {return m_userFluxFactor;}

	/// @brief  Returns the combined multiplicative factor to be applied to the flux of a Star due to all simulated processes
	double getCombinedFluxFactor() const {return m_transitFluxFactor * m_variationFluxFactor * m_noiseFluxFactor * m_userFluxFactor;}

	/// @brief Sets the multiplicative factor to be applied to the flux of a Star resulting from a planetary transit
	void setTransitFluxFactor(double transitFluxFactor) {m_transitFluxFactor = transitFluxFactor;}

	/// @brief  Sets the multiplicative factor to be applied to the flux of a Star resulting from stellar variation
	void setVariationFluxFactor(double variationFluxFactor) {m_variationFluxFactor = variationFluxFactor;}

	/// @brief  Sets the multiplicative factor to be applied to the flux of a Star resulting from stellar noise
	void setNoiseFluxFactor(double noiseFluxFactor) {m_noiseFluxFactor = noiseFluxFactor;}

	/// @brief  Sets the user defined multiplicative factor to be applied to the flux of a Star
	void setUserFluxFactor(double userFluxFactor) {m_userFluxFactor = userFluxFactor;}

private:
	double m_transitFluxFactor; ///< Multiplicative factor to be applied to the flux of a Star resulting from a planetary transit
	double m_variationFluxFactor; ///< Multiplicative factor to be applied to the flux of a Star resulting from stellar variation
	double m_noiseFluxFactor; ///< Multiplicative factor to be applied to the flux of a Star resulting from stellar noise
	double m_userFluxFactor; ///< User specified multiplicative factor to be applied to the flux of a Star
};

#endif /* _STAR_DATA_HXX_ */
