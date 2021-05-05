/*
 * TransitFluxCalculator.hxx
 *
 *  Created on: 31 Jan 2014
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Implementation of the transit light curve model developed by Mandel and Agol
///
/// This class implements the formulae of Mandel and Agol given in
/// http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/ for the case of quadratic limb
/// darkening.  The code is a conversion to C++ of C code written by Laura Kreidberg, available
/// from Eric Agol's website: http://www.astro.washington.edu/users/agol/transit.html
///
////////////////////////////////////////////////////////////////////////

#ifndef _TRANSIT_MODEL_HXX_
#define _TRANSIT_MODEL_HXX_

class TransitModel {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] radiusRatio  Ratio of the radius of the planet to that of the star
	 *  @param [in] doLimbDarkening  Flag to indicate whether or not to model limb darkening
	 *  @param [in] limbDarkeningCoeffs  Limb darkening coefficients
	 */
	TransitModel(double radiusRatio, bool doLimbDarkening, std::pair<double,double> limbDarkeningCoeffs);
	virtual ~TransitModel() {};

	/** *************************************************************************
	 *  @brief Returns the multiplicative factor by which the flux from the star
	 *  	   is reduced by the transiting planet
	 *
	 *  @param [in] starRadiusFraction  Separation between the star and the
	 *  								planet centres as a fraction of the
	 *  								stellar radius
	 * @ return Multiplicative factor by which the flux from the star
	 *  	   	is reduced by the transiting planet
	 */
	double fluxFactor(double starRadiusFraction) const;

private:

	double rc(double x, double y) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/
	double rj(double x, double y, double z, double p) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/
	double ellec(double k) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/
	double ellk(double k) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/
	double rf(double x, double y, double z) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/
	double max(double a, double b) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/
	double min(double a, double b) const; ///< See http://iopscience.iop.org/1538-4357/580/2/L171/fulltext/

	double m_radiusRatio; ///< Ratio of the radius of the planet to that of the star
	double m_doLimbDarkening; ///< Flag to indicate whether or not to model limb darkening
	std::pair<double,double> m_limbDarkeningCoeffs; ///< Limb darkening coefficients
};

#endif /* _TRANSIT_MODEL_HXX_ */
