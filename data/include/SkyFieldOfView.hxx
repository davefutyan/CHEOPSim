/*
 * SkyFieldOfView.hxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Container class for field of view data.
///
/// Contains accessors for the following information:
///		- FOV pointing direction, radius
///		- minimum magnitude for stars to be included in the simulation
///		- list of Stars
///
/// A pointer to an instance of this class is provided in the Data class
///
////////////////////////////////////////////////////////////////////////

#ifndef _SKY_FIELD_OF_VIEW_HXX_
#define _SKY_FIELD_OF_VIEW_HXX_

#include "SkyPosition.hxx"
#include "Star.hxx"

class SkyFieldOfView {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] pointingDirection  Pointing direction
	 */
	SkyFieldOfView(const SkyPosition & pointingDirection) :
		m_FOVradius(759.), m_pointingDirection(pointingDirection), m_minMagnitude(99.) {}
	virtual ~SkyFieldOfView();

	/// @brief Returns the pointing direction
	void setPointingDirection(SkyPosition pointingDirection) {m_pointingDirection = pointingDirection;}

	/// @brief Sets the radius of the field of view
	void setFOVradius(double fovRadius) {m_FOVradius = fovRadius;}

	/// @brief Sets the radius of the field of view
	void setMinMagnitude(double minMagnitude) {m_minMagnitude = minMagnitude;}

	/// @brief Returns the pointing direction
	SkyPosition getPointingDirection() const {return m_pointingDirection;}

	/// @brief Returns the radius of the field of view
	double getFOVradius() const {return m_FOVradius;}

	/// @brief Returns the minimum magnitude for stars in the field of view
	double getMinMagnitude() const {return m_minMagnitude;}

	/// @brief Returns true if the position passed as argument is located within the field of viewx
	bool isInside(const SkyPosition & position);

	/// @brief Returns the list of stars within the field of view
	vector<Star*> getStars() const {return m_stars;}

	/// @brief Adds a star to the list of stars within the field of view
	void addStar(Star * star) {m_stars.push_back(star);}

	/// @brief Resets the list of stars
	void resetStars();

private:
	double m_FOVradius; ///< Radius of the field of view in number of pixels
	SkyPosition m_pointingDirection; ///< Pointing direction
	double m_minMagnitude; ///< Minimum magnitude for stars in the field of view

	vector<Star*> m_stars; ///< List of stars within the field of view
};

#endif /* _SKY_FIELD_OF_VIEW_HXX_ */
