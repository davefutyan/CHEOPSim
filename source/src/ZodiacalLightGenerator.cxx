/*
 * ZodiacalLightGenerator.cxx
 *
 *  Created on: 6 Jan 2016
 *      Author: futyand
 */

#include <iostream>
#include <fstream>

#include "boost/date_time/gregorian/gregorian_types.hpp"

#include "telescope/include/PSFGenerator.hxx"
#include "data/include/SkyPosition.hxx"
#include "data/include/WavelengthDependence.hxx"
#include "ZodiacalLightGenerator.hxx"

void ZodiacalLightGenerator::initialize(const ModuleParams& params) {

	//Read zodiacal light flux vs wavelength taken from
	//http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c09_exposuretime08.html
	ifstream zodiacal_vs_wavelength_in(string(getenv("CHEOPS_SW"))+"/resources/zodiacal_vs_wavelength.txt", ios::in);
	if (!zodiacal_vs_wavelength_in.good()) throw runtime_error("Error opening zodiacal_vs_wavelength.txt");
    string line;
	string dummy;
	int i=0;
	vector<double> wavelengths;
	vector<double> fluxes;
	while (getline(zodiacal_vs_wavelength_in, line)) {
		if(line[0] != '#') {
			double flux, wavelength;
    		stringstream(line) >> wavelength >> dummy >> flux >> dummy;
    		wavelengths.push_back(wavelength);
    		fluxes.push_back(flux);
			i++;
		}
	}
	zodiacal_vs_wavelength_in.close();

	//Calculate the integral of the distribution of zodiacal light flux vs wavelength
	//weighted by 1/photonEnergy in order to calculate the integrated photon flux
	//Only consider the same wavelength range as considered for the target star
	//Wavelength in angstroms, 100 angstrom steps
	m_zodiacalLightFlux = 0.;
	for (unsigned wavelength=Star::kWavelengthLowBound*10.; wavelength<(Star::kWavelengthLowBound*10. + 100*Star::kNWavelength); wavelength+=100) {
		double flux = 0.;
		//Find the entries in the vectors read from the file immediately above and below the current wavelength
		for (unsigned i=1; i<wavelengths.size(); i++) {
			if (wavelengths[i]>wavelength) {
				//interpolate between the fluxes read from the file to the current wavelength
				double frac = (wavelength - wavelengths[i-1])/(wavelengths[i] - wavelengths[i-1]);
				flux = fluxes[i-1]*(1.-frac) + fluxes[i]*frac;
				break;
			}
		}
		double photonEnergy = WavelengthDependence::kPlanck * WavelengthDependence::kSpeedOfLight / (wavelength*1.e-10);
		m_zodiacalLightFlux += (flux/photonEnergy)*100.;
	}
	m_zodiacalLightFlux *= PSFGenerator::kTelescopeArea;

	//Read dependence on the difference in ecliptic latitude and longitude between the pointing direction and the sun, taken from
	//http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c09_exposuretime08.html
	ifstream eclipticDependence_in(string(getenv("CHEOPS_SW"))+"/resources/zodiacal_vs_ecliptic_coords.txt", ios::in);
	if (!eclipticDependence_in.good()) throw runtime_error("Error opening zodiacal_vs_ecliptic_coords.txt");
	int ilong=-1;
	while (getline(eclipticDependence_in, line)) {
		if(line[0] != '#') {
			stringstream ss;
			ss << line;
    		if (ilong==-1) {
    			for (int ilat=0; ilat<7; ilat++) {
    				ss >> m_eclipticLatitude[ilat];
    			}
    		} else {
    			ss >> m_eclipticLongitude[ilong];
    			for (int ilat=0; ilat<7; ilat++) {
    				ss >> m_eclipticDependence[ilong][ilat];
    			}
    		}
			ilong++;
		}
	}
	eclipticDependence_in.close();

}

void ZodiacalLightGenerator::doBegin(Data* data, bool fullFrame) {

	//Convert the pointing direction (equitorial coordinates) to ecliptic polar coordinates
	SkyPosition pointingDirection = data->getFieldOfView()->getPointingDirection();
	m_pointingEclipticLatitude = pointingDirection.getEclipticLatitude();
	m_pointingEclipticLongitude = pointingDirection.getEclipticLongitude();

}

void ZodiacalLightGenerator::process(Data* data, int timeStep,
		bool fullFrame) const {

	//Get the julian date for the current time step
	TimeConfiguration timeConf = data->getTimeConfiguration();
	boost::posix_time::ptime currentTime = timeConf.getStartTime() +
										 boost::posix_time::seconds(static_cast<long>(lround(timeConf.getTimeSinceStart(timeStep))));
	double fractionOfDay = static_cast<double>(currentTime.time_of_day().total_seconds())/(3600.*24.);
	double julianDate = 2400000 + static_cast<double>(currentTime.date().modjulian_day())+fractionOfDay;

	//Calculate the ecliptic longitude of the sun for the current julian date
	double sunEclipticLong = sunEclipticLongitude(julianDate);
	//cout << m_pointingEclipticLongitude << " " << m_pointingEclipticLatitude << " " << sunEclipticLong << endl;

	//Calculate the Vband magnitude of the zodiacal light given the sun position relative to the pointing direction
	double magVband = zodiacalVbandMagnitude(sunEclipticLong);

	//Calculate how much to modify the flux due to the ecliptic coordinate dependency
	//See http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c09_exposuretime08.html:
	//"the values in Figure 9.1 and Table 9.3 correspond to the high value of V-band surface brightness of 22.1 mag arcsec-2"
	double deltaMagnitude = magVband - 22.1;
	double zodiacalLightFlux = m_zodiacalLightFlux * pow(10,-0.4*deltaMagnitude);
	//cout << magVband << " " << m_zodiacalLightFlux << " " << zodiacalLightFlux << endl;

	//Add the zodiacal light uniformly to each pixel of the image
	Image* image = data->getImages().back();
	int ymin = data->doChargeTransferEoL() ? 0. : image->getYOffset(); //Include pixels below the subarray if end of life CTI is to be simulated
	for (int ix=image->getXOffset(); ix<image->getXOffset()+image->getXDim(); ix++) {
		for (int iy=ymin; iy<image->getYOffset()+image->getYDim(); iy++) {
			image->incrementPixelValue(ix,iy,zodiacalLightFlux*timeConf.getExposureTimeAsDouble());
		}
	}

	//Set the truth information
	image->getTruthData()->setZodiacalLight(zodiacalLightFlux*timeConf.getExposureTimeAsDouble());

}

double ZodiacalLightGenerator::sunEclipticLongitude(double julianDate) const {

	//Calulation taken from http://www.meteorobs.org/maillist/msg09197.html, with modifications from
	//http://www.stargazing.net/kepler/sun.html and https://en.wikipedia.org/wiki/Position_of_the_Sun#Ecliptic_coordinates
	//Checked using http://www.satellite-calculations.com/Satellite/suncalc.htm and http://www.imo.net/data/solar

	//Julian Centuries of 36525 ephemeris days from the epoch J2000.0 (2000 January 1.5 TD):
	double T = (julianDate - 2451545.0) / 36525.;
	//cout << julianDate << " " << T << " " << T*36525 << endl;

	//Geometric mean longitude of the Sun, referred to the mean equinox of the date:
	double L0 = 280.46645 + 36000.76983*T + 0.0003032*pow(T,2);
	SkyPosition::setRange0to360(L0);

	//Mean anomaly of the Sun:
	double M = 357.52910 + 35999.05030*T - 0.0001559*pow(T,2) - 0.00000048*pow(T,3);
	SkyPosition::setRange0to360(M);

	//Sun's equation of center:
	double C = (1.914600 - 0.004817*T - 0.000014*pow(T,2)) * sind(M) + (0.01993 - 0.000101*T) * sind(2.*M) + 0.000290 * sind(3.*M);

	//Ecliptic longitude of the Sun:
	double S = L0 + C;
	SkyPosition::setRange0to360(S);

	//cout << L0 << " " << M << " " << S << endl;

	//The Sun's geometric equatorial coordinates:
//	double y = cosd(SkyPosition::kEclipticObliquity)*sind(S);
//	double x = cosd(S);
//	double RA = atan2(y,x)*180./M_PI;
//	double dec = asin(sind(SkyPosition::kEclipticObliquity)*sind(S))*180./M_PI;
//	cout << "RA and dec of Sun: " << RA << " " << dec << endl;

	return S;

}

double ZodiacalLightGenerator::zodiacalVbandMagnitude(double sunEclipticLong) const {

	//Calculate the difference in ecliptic latitude and longitude between the pointing direction and the sun
	double deltaLong = fabs(m_pointingEclipticLongitude - sunEclipticLong);
	if (deltaLong>180.) deltaLong = 360.-deltaLong;
	double deltaLat = fabs(m_pointingEclipticLatitude);

	//Determine which row and column to read from the table for the ecliptic coordinate dependency of the zodiacal light
	int ilong = 0;
	for (int i=1; i<13; i++) {
		if (deltaLong >= m_eclipticLongitude[i]) ilong = i;
	}
	int ilat = 0;
	for (int i=1; i<7; i++) {
		if (deltaLat >= m_eclipticLatitude[i]) ilat = i;
	}

	//Obtain the magnitudes of the four nearest long,lat grid points in the table
	//Protect against boundaries by effectively add an extra row and column to the table which are copies
	//of the last row and column.
	int ilongplus1 = ilong==12 ? 12 : ilong+1;
	int ilatplus1 = ilat==6 ? 6 : ilat+1;
	double q11 = m_eclipticDependence[ilong][ilat];
	double q12 = m_eclipticDependence[ilongplus1][ilat];
	double q21 = m_eclipticDependence[ilong][ilatplus1];
	double q22 = m_eclipticDependence[ilongplus1][ilatplus1];

	//Exit with error message if pointing direction is in solar avoidance zone.
	//This will be protected against in the web interface.
	if (q11==0. || q12==0. || q21==0. || q22==0.) {
		string message = "Error in ZodiacalLightGenerator: pointing direction is ";
		message += "within solar avoidance zone on the start date specified. ";
		message += "Please use a different start time or pointing direction, ";
		message += "or switch off the ZodiacalLightGenerator module.";
		throw runtime_error(message);
	}

	//Interpolate between the four nearest points. Step size in the table is 15 degrees.
	double longfrac = (deltaLong - m_eclipticLongitude[ilong])/15.;
	double latfrac = (deltaLat - m_eclipticLatitude[ilat])/15.;
	//cout << ilat << " " << ilong << " " << q11 << " " << q12 << " " << q21 << " " << q22 << " " << longfrac << " " << latfrac << endl;
	return bilinearInterpolation(q11,q12,q21,q22,latfrac,longfrac);

}

double ZodiacalLightGenerator::bilinearInterpolation(double q11, double q12, double q21, double q22, double xfrac, double yfrac) const {

	double r1 = (1.-xfrac)*q11 + xfrac*q21;
	double r2 = (1.-xfrac)*q12 + xfrac*q22;
	return (1.-yfrac)*r1 + yfrac*r2;

}
