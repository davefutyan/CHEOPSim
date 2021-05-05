/*
 * PSF.cxx
 *
 *  Created on: 4 Dec 2014
 *      Author: futyand
 */

#include <numeric>      // std::accumulate

#include "PSF.hxx"

PSF::PSF(vector<double> x, vector<double> y, double val) : m_xPositions(x), m_yPositions(y), m_flux(val) {}

PSF& PSF::operator +=(const PSF& psf) {
	this->m_xPositions.insert(this->m_xPositions.end(),psf.m_xPositions.begin(),psf.m_xPositions.end());
	this->m_yPositions.insert(this->m_yPositions.end(),psf.m_yPositions.begin(),psf.m_yPositions.end());
	this->m_flux += psf.m_flux;
	return *this;
}

double PSF::getMeanXPosition() const {
	return (accumulate(m_xPositions.begin(), m_xPositions.end(), 0.)/m_xPositions.size());
}

double PSF::getMeanYPosition() const {
	return (accumulate(m_yPositions.begin(), m_yPositions.end(), 0.)/m_yPositions.size());
}
