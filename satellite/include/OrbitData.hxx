/*
 * OrbitData.hxx
 *
 *  Created on: 4 Apr 2016
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Class to define the position and velocity of the spacecraft in its orbit
///
////////////////////////////////////////////////////////////////////////

#ifndef DATA_INCLUDE_ORBITDATA_HXX_
#define DATA_INCLUDE_ORBITDATA_HXX_

#include <cmath>

#include "Utc.hxx"

class OrbitData {
public:
	/** *************************************************************************
	 *  @brief Constructor for orbit data
	 *
	 *  @param [in] utcTime  UTC time
	 *  @param [in] x  x component of orbit position (optional)
	 *  @param [in] y  y component of orbit position (optional)
	 *  @param [in] z  z component of orbit position (optional)
	 *  @param [in] latitude  latitude of orbit position (optional)
	 *  @param [in] longitude  longitude of orbit position (optional)
	 */
	OrbitData(UTC utcTime, double x=0., double y=0., double z=0., double latitude=0., double longitude=0.) :
		m_utc(utcTime), m_x(x), m_y(y), m_z(z), m_vx(0.), m_vy(0.), m_vz(0.),
		m_latitude(latitude), m_longitude(longitude), m_moonAngle(0.), m_sunAngle(0.) {};
	virtual ~OrbitData() {};

	/** *************************************************************************
	 *  @brief Sets the orbit position
	 *
	 *  @param [in] x  x component of orbit position
	 *  @param [in] y  y component of orbit position
	 *  @param [in] z  z component of orbit position
	 */
	void setPosition(double x, double y, double z) {m_x=x; m_y=y; m_z=z;}

	/** *************************************************************************
	 *  @brief Sets the latitude and longitude of the orbit position
	 *
	 *  @param [in] latitude  latitude of orbit position
	 *  @param [in] longitude  longitude of orbit position
	 */
	void setLatLong(double latitude, double longitude) {m_latitude=latitude; m_longitude=longitude;}

	/** *************************************************************************
	 *  @brief Sets the orbit velocity
	 *
	 *  @param [in] vx  x component of orbit velocity
	 *  @param [in] vy  y component of orbit velocity
	 *  @param [in] vz  z component of orbit velocity
	 */
	void setVelocity(double vx, double vy, double vz) {m_vx=vx; m_vy=vy; m_vz=vz;}

	/** *************************************************************************
	 *  @brief Sets the angle between the pointing direction and the moon
	 *
	 *  @param [in] angle  Angle between the pointing direction and the moon
	 */
	void setMoonAngle(double angle) {m_moonAngle=angle;}

	/** *************************************************************************
	 *  @brief Sets the angle between the pointing direction and the sun
	 *
	 *  @param [in] angle  Angle between the pointing direction and the sun
	 */
	void setSunAngle(double angle) {m_sunAngle=angle;}

	UTC utc() const {return m_utc;} ///< @brief Returns the UTC time
	double x() const {return m_x;} ///< @brief Returns the x component of orbit position
	double y() const {return m_y;} ///< @brief Returns the y component of orbit position
	double z() const {return m_z;} ///< @brief Returns the z component of orbit position
	double vx() const {return m_vx;} ///< @brief Returns the x component of orbit velocity
	double vy() const {return m_vy;} ///< @brief Returns the y component of orbit velocity
	double vz() const {return m_vz;} ///< @brief Returns the z component of orbit velocity
	double latitude() const {return m_latitude;} ///< @brief Returns the latitude of orbit position
	double longitude() const {return m_longitude;} ///< @brief Returns the longitude of orbit position
	double moonAngle() const {return m_moonAngle;} ///< @brief Returns the angle between the pointing direction and the moon
	double sunAngle() const {return m_sunAngle;} ///< @brief Returns the angle between the pointing direction and the sun

private:
	UTC m_utc; ///< UTC time
	double m_x; ///< x component of orbit position
	double m_y; ///< y component of orbit position
	double m_z; ///< z component of orbit position
	double m_vx; ///< x component of orbit velocity
	double m_vy; ///< y component of orbit velocity
	double m_vz; ///< z component of orbit velocity
	double m_latitude; ///< latitude of orbit position
	double m_longitude; ///< longitude of orbit position
	double m_moonAngle; ///< Angle between the pointing direction and the moon
	double m_sunAngle; ///< Angle between the pointing direction and the sun
};

#endif /* DATA_INCLUDE_ORBITDATA_HXX_ */
