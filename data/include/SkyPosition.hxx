/*
 * SkyPosition.hxx
 *
 *  Created on: Dec 19, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Class to specify a position in sky coordinates (RA,dec).
///
/// Includes a method to return the separation between two SkyPositions
////////////////////////////////////////////////////////////////////////

#ifndef _SKY_POSITION_HXX_
#define _SKY_POSITION_HXX_

#include <cmath>

class SkyPosition {
public:

	static constexpr double kEclipticObliquity = 23.4392911; ///< Obliquity of the ecliptic plane

	/// @brief Converts an angle in degrees to a value in the range 0 to 360
	static void setRange0to360(double & angle);

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] ra  Right ascension in degrees
	 *  @param [in] dec  Declination in degrees
	 *  @param [in] raErr  Right ascension error in degrees (optional)
	 *  @param [in] decErr  Declination error in degrees (optional)
	 */
	SkyPosition(double ra, double dec, double raErr=0., double decErr=0.) : m_ra(ra),m_dec(dec),m_raErr(raErr),m_decErr(decErr) {};

	/// @brief Copy constructor
	SkyPosition(const SkyPosition & skyPosition) : m_ra(skyPosition.m_ra),m_dec(skyPosition.m_dec),m_raErr(skyPosition.m_raErr),m_decErr(skyPosition.m_decErr) {};
	virtual ~SkyPosition() {};

	/// @brief Sets the right ascension
	void setRightAscension(double ra) {m_ra = ra;}

	/// @brief Sets the declination
	void setDeclination(double dec) {m_dec = dec;}

	/// @brief Returns the right ascension
	double getRightAscension() const {return m_ra;}

	/// @brief Returns the declination
	double getDeclination() const {return m_dec;}

	/// @brief Returns the right ascension
	double getRightAscensionError() const {return m_raErr;}

	/// @brief Returns the declination
	double getDeclinationError() const {return m_decErr;}

	/// @brief Returns the ecliptic latitude
	double getEclipticLatitude();

	/// @brief Returns the ecliptic longitude
	double getEclipticLongitude();

	/// @brief Returns the separation between two SkyPositions
	double separation(const SkyPosition & pos) const;

private:
	double m_ra; ///< Right ascension
	double m_dec; ///< Declination
	double m_raErr; ///< Right ascension error
	double m_decErr; ///< Declination error
};

#endif /* _SKY_POSITION_HXX_ */
