/*
 * SkyPosition.cxx
 *
 *  Created on: Jan 7, 2016
 *      Author: futyand
 */

#include <iostream>

#include "SkyPosition.hxx"

void SkyPosition::setRange0to360(double& angle) {
	if (angle > 360.) {
		angle -= 360.*(floor(angle/360.));
	} else if (angle < 0.) {
		angle += 360.*(ceil(-angle/360.));
	}
}

double SkyPosition::getEclipticLatitude() {
	//Formula from https://books.google.co.uk/books?id=DwJfCtzaVvYC&printsec=frontcover#v=onepage&q&f=false p.42
	//("Practical Astronomy With Your Calculator" by Peter Duffett-Smith)
	//Checked using https://ned.ipac.caltech.edu/forms/calculator.html
	double epsilon_rads = kEclipticObliquity * M_PI/180.;
	double ra_rads = m_ra/3600. * M_PI/180.;
	double dec_rads = m_dec/3600. * M_PI/180.;
	return asin(sin(dec_rads)*cos(epsilon_rads) - cos(dec_rads)*sin(epsilon_rads)*sin(ra_rads)) * 180./M_PI;
}

double SkyPosition::getEclipticLongitude() {
	//Formula from https://books.google.co.uk/books?id=DwJfCtzaVvYC&printsec=frontcover#v=onepage&q&f=false p.42
	//("Practical Astronomy With Your Calculator" by Peter Duffett-Smith)
	//Checked using https://ned.ipac.caltech.edu/forms/calculator.html
	double epsilon_rads = kEclipticObliquity * M_PI/180.;
	double ra_rads = m_ra/3600. * M_PI/180.;
	double dec_rads = m_dec/3600. * M_PI/180.;
	double eclipticLongitude = atan2(sin(ra_rads)*cos(epsilon_rads)+tan(dec_rads)*sin(epsilon_rads),cos(ra_rads)) * 180/M_PI;
	setRange0to360(eclipticLongitude);
	return eclipticLongitude;
}

double SkyPosition::separation(const SkyPosition & pos) const {
	double dec_rads = m_dec/3600. * M_PI/180.;
	return std::sqrt(pow((m_ra-pos.m_ra)*cos(dec_rads),2) + pow(m_dec-pos.m_dec,2));
}
