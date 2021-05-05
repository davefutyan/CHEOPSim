/*
 * ZodiacalLightGenerator.hxx
 *
 *  Created on: 6 Jan 2016
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to calculate the Zodiacal light flux incident on the telescope
///		   based on the pointing direction and the position of the sun on the current date
///
////////////////////////////////////////////////////////////////////////

#ifndef ZODIACALLIGHTGENERATOR_HXX_
#define ZODIACALLIGHTGENERATOR_HXX_

#include "simulator/include/Module.hxx"

class ZodiacalLightGenerator: public Module {
public:

	ZodiacalLightGenerator() : Module("ZodiacalLightGenerator",timeLoop), m_zodiacalLightFlux(0.),
							   m_pointingEclipticLatitude(0.),m_pointingEclipticLongitude(0.) {};
	virtual ~ZodiacalLightGenerator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Returns the ecliptic longitude of the Sun at the specified date
	 *
	 *  @param [in] julianDate  Julian date
	 *  @return Ecliptic longitude of the Sun
	 */
	double sunEclipticLongitude(double julianDate) const;

	/** *************************************************************************
	 *  @brief Returns the Vband magnitude of the zodiacal light given the sun
	 *  	   position relative to the pointing direction
	 *
	 *  @param [in] sunEclipticLong  Ecliptic longitude of the Sun
	 *  @return Vband magnitude of the zodiacal light
	 */
	double zodiacalVbandMagnitude(double sunEclipticLong) const;

	/// @brief Returns the sine of an angle provided in degrees
	double sind(double angle) const {return sin(angle * M_PI/180.);}

	/// @brief Returns the cosine of an angle provided in degrees
	double cosd(double angle) const {return cos(angle * M_PI/180.);}

	/** *************************************************************************
	 *  @brief Performs a bilinear interpolation between values at 4 nearest
	 *  	   grid points
	 *
	 *  @param [in] q11  Value for nearest grid pixel below and left of
	 *  				 the target position
	 *  @param [in] q12  Value for nearest grid pixel above and left of
	 *  				 the target position
	 *  @param [in] q21  Value for nearest grid pixel below and right of
	 *  				 the target position
	 *  @param [in] q22  Value for nearest grid pixel above and right of
	 *  				 the target position
	 *  @param [in] xfrac  Target position in the horizontal direction as a
	 *  				   fraction of the distance between the two nearest
	 *  				   vertical grid lines
	 *  @param [in] yfrac  Target position in the vertical direction as a
	 *  				   fraction of the distance between the two nearest
	 *  				   horizontal grid lines
	 */
	double bilinearInterpolation(double q11, double q12, double q21, double q22, double xfrac, double yfrac) const;

	/// Zodiacal light flux in photons/s/pixel corresponding to a V-band surface brightness of 22.1 mag arcsec-2
	double m_zodiacalLightFlux;

	/// Ecliptic latitude of pointing direction
	double m_pointingEclipticLatitude;

	/// Ecliptic longitude of pointing direction
	double m_pointingEclipticLongitude;

	/// List of ecliptic latitudes for which ZL magnitude is provided
	double m_eclipticLatitude[7];

	/// List of ecliptic longitudes for which ZL magnitude is provided
	double m_eclipticLongitude[13];

	/// ZL magnitude for each point in the ecliptic latitude-longitude grid defined by m_eclipticLatitude and m_eclipticLatitude
	double m_eclipticDependence[13][7];
};

#endif /* ZODIACALLIGHTGENERATOR_HXX_ */
