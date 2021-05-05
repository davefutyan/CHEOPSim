/*
 * SkyFieldOfView.cxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

#include "SkyFieldOfView.hxx"

#include <stdexcept>

bool SkyFieldOfView::isInside(const SkyPosition& position) {
	return m_pointingDirection.separation(position)<m_FOVradius;
}

SkyFieldOfView::~SkyFieldOfView() {
    for (vector<Star*>::const_iterator it = m_stars.begin(); it!=m_stars.end(); ++it) delete (*it);
}

void SkyFieldOfView::resetStars() {
    for (vector<Star*>::const_iterator it = m_stars.begin(); it!=m_stars.end(); ++it) delete (*it);
    m_stars.clear();
}
