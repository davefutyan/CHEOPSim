/*
 * FluxConverter.hxx
 *
 *  Created on: 2 Jun 2020
 *      Author: futyand
 */

#ifndef SOURCE_INCLUDE_FLUXCONVERTER_HXX_
#define SOURCE_INCLUDE_FLUXCONVERTER_HXX_

#include "data/include/WavelengthDependence.hxx"

using namespace std;

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This class provides methods to calculate the difference between
 *  	   CHEOPS magnitude and Gaia or V-band magnitude for a star with
 *  	   a given effective temperature, using the formulae documented in
 *  	   https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion
 *  	   and using parameters calculated using
 *  	   https://svn.isdc.unige.ch/svn-cheops/06_cheopsim/software/CHEOPSim/trunk/resources/FluxConversion.html
 */

class FluxConverter {
public:

	static constexpr double kDefaultBlackBodyOffset = -0.004795; ///< Default value of empirically derived offset correction for (CHEOPS mag - Gaia mag) for TEff>7200K, when black body spectra are used to model the stellar spectra
	static constexpr double kDefaultFitIntercept = 0.94604; ///< Default value of intercept of linear fit to (observed flux / predicted flux) vs effective temperature
	static constexpr double kDefaultFitSlope = -2.0074e-06; ///< Default value of slope of linear fit to (observed flux / predicted flux) vs effective temperature

	/** *************************************************************************
	 *  @brief Constructor: pass in the input parameters.
	 *
	 *  @param [in] blackBodyOffset		Empirically derived offset correction
	 *  								for (CHEOPS mag - Gaia mag) for TEff>7200K,
	 *  								when black body spectra are used to model
	 *  								the stellar spectra
	 *  @param [in] fitIntercept		Intercept of linear fit to
	 *  								(observed flux / predicted flux)
	 *  								vs effective temperature
	 *  @param [in] fitIntercept		Slope of linear fit to
	 *  								(observed flux / predicted flux)
	 *  								vs effective temperature
	 */
	FluxConverter(double blackBodyOffset=kDefaultBlackBodyOffset, double fitIntercept=kDefaultFitIntercept, double fitSlope=kDefaultFitSlope);

	virtual ~FluxConverter() {};

	/** *************************************************************************
	 *  @brief Returns the difference between CHEOPS magnitude and Gaia or V-band
	 *  	   magnitude for a star with the specified effective temperature.
	 *  	   The Gaia band is used rather than V-band if the gaiaBand
	 *  	   configuration parameter is set to true in runCHEOPSim.xml or
	 *  	   runFluxConversion.xml.
	 *
	 *  @param [in] refBand 			  Reference passband: GaiaBand or Vband
	 *  @param [in] effectiveTemperature  Black body effective temperature
	 *  @param [in] wavelengthDependence  Pointer to instance of
	 *  								  WavelengthDependence class
	 *  @param [in] target				  Flag to indicate whether or not
	 *  								  the Star is the target
	 */
	double cheopsMinusRefBandMagnitude(WavelengthDependence::REFERENCE_BAND refBand, double effectiveTemperature,
									   const WavelengthDependence * wavelengthDependence, bool target=false);

	/** *************************************************************************
	 *  @brief Returns the multiplicative correction factor to be applied to the
	 *  	   predicted flux to account for the empirical discrepancy
	 *  	   between observed and predicted fluxes
	 *
	 *  @param [in] effectiveTemperature  Black body effective temperature
	 */
	double fluxCorrectionFactor(double effectiveTemperature) {return m_fitIntercept + effectiveTemperature*m_fitSlope;}

private:

	double m_blackBodyOffset; ///< Empirically derived offset correction for (CHEOPS mag - Gaia mag) for TEff>7200K, when black body spectra are used to model the stellar spectra
	double m_fitIntercept; ///< Intercept of linear fit to (observed flux / predicted flux) vs effective temperature
	double m_fitSlope; ///< Slope of linear fit to (observed flux / predicted flux) vs effective temperature

};

#endif /* SOURCE_INCLUDE_FLUXCONVERTER_HXX_ */
