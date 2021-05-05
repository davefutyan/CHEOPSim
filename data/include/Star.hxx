/*
 * Star.hxx
 *
 *  Created on: Dec 19, 2013
 *      Author: futyand
 */

#ifndef _STAR_HXX_
#define _STAR_HXX_

#include <vector>
#include <string>
using namespace std;

#include "SkyPosition.hxx"
#include "StarData.hxx"
#include "Bjd.hxx"

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Container class for Star data
///
/// Contains accessors for:
/// 	1. magnitude, position, mean flux (total and in wavelength bins)
///		2. vector<StarData*>:  Contains time dependent information for the star.
///							   Each vector element corresponds to a
///                            different time step (the vector is initially empty).
///		3. x and y coordinates in the focal plane for each second of the current
///		   time step (time dependence is due to jitter and roll angle)
///
/// The vector of Stars in the simulation can be accessed from the
/// Data class via Data::getFieldOfView() and SkyFieldOfView::getStars()
///
////////////////////////////////////////////////////////////////////////

class Star {
public:

	/// @brief Luminosity of the Sun in Watts
	static constexpr double kLSun = 3.828E26;
	/// @brief Stefan Boltzmann constant in kg s^-3 K^-4
	static constexpr double kSigma = 5.670367E-8;
	/// @brief Number of wavelength bins (bin width 10nm)
	static const int kNWavelength = 77;
	/// @brief Lower bound on wavelength range
	static constexpr double kWavelengthLowBound = 330.;

	/// @brief Spectral type of the star
	enum spectralType{Undefined,O7,O8,O9,
								B0,B1,B2,B3,B4,B5,B6,B7,B8,B9,
								A0,A1,A2,A3,A4,A5,A6,A7,A8,A9,
								F0,F1,F2,F3,F4,F5,F6,F7,F8,F9,
								G0,G1,G2,G3,G4,G5,G6,G7,G8,G9,
								K0,K1,K2,K3,K4,K5,K6,K7,K8,K9,
								M0,M1,M2,M3,M4,M5,M6,M7,M8,M9};
	/// @brief Spectral type of the star as a string
	static const string spectralTypeString[64];

	/// @brief Returns the effective temperature corresponding to a given spectral type
	static double getEffectiveTemperature(spectralType specType);

	/** *************************************************************************
	 *  @brief Constructor. Either specType or effectiveTemperature must be specified
	 *
	 *  @param [in] pos  Position (RA, dec) of the star in the sky
	 *  @param [in] VbandMag  V-band magnitude of the star
	 *  @param [in] cheopsMag  CHEOPS magnitude of the star
	 *  @param [in] gaiaMag  Gaia magnitude of the star
	 *  @param [in] specType  Spectral type of the star
	 *  @param [in] effectiveTemperature  Effective temperature of the star (optional)
	 *  @param [in] id  Catalogue ID of the star (optional)
	 */
	Star(SkyPosition pos, double VbandMag, double cheopsMag, double gaiaMag, spectralType specType, double effectiveTemperature=0., string id="");
	virtual ~Star();

	/// @brief Returns the V-band magnitude of the star
	double getVbandMagnitude() const {return m_Vmag;}

	/// @brief Returns the CHEOPS magnitude of the star
	double getCheopsMagnitude() const {return m_cheopsMag;}

	/// @brief Returns the Gaia band magnitude of the star
	double getGaiaMagnitude() const {return m_gaiaMag;}

	/// @brief Sets the error on the V-band magnitude of the star
	void setVbandMagnitudeError(double VbandMagErr) {m_VmagErr = VbandMagErr;}

	/// @brief Sets the error on the CHEOPS magnitude of the star
	void setCheopsMagnitudeError(double cheopsMagErr) {m_cheopsMagErr = cheopsMagErr;}

	/// @brief Sets the error on the V-band magnitude of the star
	void setGaiaMagnitudeError(double gaiaMagErr) {m_gaiaMagErr = gaiaMagErr;}

	/// @brief Returns the error on the V-band magnitude of the star
	double getVbandMagnitudeError() const {return m_VmagErr;}

	/// @brief Returns the error on the CHEOPS magnitude of the star
	double getCheopsMagnitudeError() const {return m_cheopsMagErr;}

	/// @brief Returns the error on the Gaia band magnitude of the star
	double getGaiaMagnitudeError() const {return m_gaiaMagErr;}

	/// @brief Returns the position (RA, dec) of the star in the sky
	const SkyPosition & getSkyPosition() const {return m_pos;}

	/// @brief Returns the spectral type of the star
	spectralType getSpectralType() const {return m_spectralType;}

	/// @brief Returns the spectral type of the star as a string
	string getSpectralTypeString() const {return spectralTypeString[m_spectralType];}

	/// @brief Returns the effective temperature of the star
	double getEffectiveTemperature() const {return m_effectiveTemperature;}

	/// @brief Sets the error on the effective temperature of the star
	void setEffectiveTemperatureError(double effectiveTemperatureErr) {m_effectiveTemperatureErr = effectiveTemperatureErr;}

	/// @brief Returns the error on the effective temperature of the star
	double getEffectiveTemperatureError() const {return m_effectiveTemperatureErr;}

	/// @brief Returns the catalogue ID of the star
	string getCatalogueId() const {return m_id;}

	/// @brief Returns the radius of the star in solar radii
	double getRadius() const {return m_radius;}

	/// @brief Returns the mass of the star in solar masses
	double getMass() const {return m_mass;}

	/// @brief Sets the central time of the first transit
	void setTransitTime(BJD transitTime) {m_transitTime = transitTime;}

	/// @brief Sets the time interval in days between two consecutive transits
	void setTransitPeriod(double period) {m_transitPeriod = period;}

	/// @brief Returns the central time of the first transit
	BJD getTransitTime() const {return m_transitTime;}

	/// @brief Returns the time interval in days between two consecutive transits
	double getTransitPeriod() const {return m_transitPeriod;}

	/// @brief Returns the limb darkening coefficients, as calculated using https://github.com/nespinoza/limb-darkening
	pair<double,double> getLimbDarkeningCoefficients() const {return m_limbDarkeningCoeffs;}

	/// @brief Returns the total mean flux integrated over all wavelengths
	double getMeanFlux() const;

	/// @brief Returns the total mean photon flux integrated over all wavelengths
	double getMeanPhotonFlux() const;

	/// @brief Sets the mean flux in bins of wavelength
	void setMeanFlux(double flux[kNWavelength]) {for (int i=0; i<kNWavelength; i++) m_meanFlux[i] = flux[i];}

	/// @brief Sets the mean photon flux in bins of wavelength
	void setMeanPhotonFlux(double photonFlux[kNWavelength]) {for (int i=0; i<kNWavelength; i++) m_meanPhotonFlux[i] = photonFlux[i];}

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Returns the mean flux in the wavelength bin (width 10nm)
	///        corresponding to the specified wavelength
	///        Units: ergs/s/cm2/10nm
	/// @param [in] wavelength:  Wavelength at which flux is to be evaluated
	double getMeanFlux(double wavelength) const;

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Returns the mean flux in the specified wavelength bin (width 10nm)
	///        Units: ergs/s/cm2/10nm
	/// @param [in] iWavelength:  Wavelength bin at which flux is to be evaluated
	double getMeanFlux(int iWavelength) const {return m_meanFlux[iWavelength];}

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Returns the mean photon flux in the wavelength bin (width 10nm)
	///        corresponding to the specified wavelength
	///        Units: photons/s/cm2/10nm
	/// @param [in] wavelength:  Wavelength at which flux is to be evaluated
	double getMeanPhotonFlux(double wavelength) const;

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Returns the mean photon flux in the specified wavelength bin (width 10nm)
	///        Units: photons/s/cm2/10nm
	/// @param [in] iWavelength:  Wavelength bin at which flux is to be evaluated
	double getMeanPhotonFlux(int iWavelength) const {return m_meanPhotonFlux[iWavelength];}

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Initializes the StarData for each time step with default values
	/// @param [in] nTimeSteps:  Number of time steps (defines the size of the vector)
	void initializeTimeSeries(int nTimeSteps);

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Returns time dependent information for the star.
	///        Each vector element corresponds to a different time step.
	vector<StarData*> getTimeSeries() const;

	/// @brief  Sets the list of x coordinates of the position in the focal plane for each second of the current time step of the simulation (centre is 0.0)
	void setFocalPlaneX(vector<double> focalPlaneX) {m_focalPlaneX.clear(); m_focalPlaneX = focalPlaneX;}

	/// @brief  Sets the list of y coordinates of the position in the focal plane for each second of the current time step of the simulation (centre is 0.0)
	void setFocalPlaneY(vector<double> focalPlaneY) {m_focalPlaneY.clear(); m_focalPlaneY = focalPlaneY;}

	/// @brief  Returns the list of x coordinates of the position in the focal plane for each second of the current time step of the simulation (centre is 0.0)
	vector<double> getFocalPlaneX() const {return m_focalPlaneX;}

	/// @brief  Returns the list of y coordinates of the position in the focal plane for each second of the current time step of the simulation (centre is 0.0)
	vector<double> getFocalPlaneY() const {return m_focalPlaneY;}

	/// @brief Sets halo flux integrated over 200x200 pixels, per unit incident flux
	void setHaloFlux(double haloFlux) {m_haloFlux = haloFlux;}

	/// @brief Returns halo flux intergated over 200x200 pixels, per unit incident flux
	double getHaloFlux() const {return m_haloFlux;}

private:
	SkyPosition m_pos; ///< Position of the star in the sky
	double m_Vmag; ///< V-band magnitude of the star
	double m_cheopsMag; ///< CHEOPS magnitude of the star
	double m_gaiaMag; ///< Gaia magnitude of the star
	double m_VmagErr; ///< Error on V-band magnitude of the star
	double m_cheopsMagErr; ///< Error on CHEOPS magnitude of the star
	double m_gaiaMagErr; ///< Error on Gaia magnitude of the star
	spectralType m_spectralType; ///< Spectral type of the star
	string m_spectralTypeString; ///< Spectral type of the star as a string
	string m_id; ///< Catalogue ID of the star
	double m_radius; ///< Radius of the star in solar masses
	double m_mass; ///< Mass of the star in solar radii
	double m_effectiveTemperature; ///< Effective temperature of the star
	double m_effectiveTemperatureErr; ///< Error on effective temperature of the star
	BJD m_transitTime; ///< Central time of first transit
	double m_transitPeriod; ///< Number of days between two consecutive transits
	pair<double,double> m_limbDarkeningCoeffs; ///< Quadratic coefficients for limb darkening
	double m_meanFlux[kNWavelength]; ///< Mean flux of the star binned in wavelength in ergs/s/cm2/10nm
	double m_meanPhotonFlux[kNWavelength]; ///< Mean photon flux of the star binned in wavelength in photons/s/cm2/10nm
	double m_wavelengthBinLowEdge[kNWavelength]; ///< Low edges of the wavelength bins
	vector<StarData*> m_timeData; ///< Time series data for the star (multiplicative factors to be applied to the flux as a function of time)
	vector<double> m_focalPlaneX; ///< List of x coordinates of the position in the focal plane for each second of the current time step of the simulation (centre is 0.0)
	vector<double> m_focalPlaneY; ///< List of y coordinates of the position in the focal plane for each second of the current time step of the simulation (centre is 0.0)
	double m_haloFlux; ///< halo flux intergated over 200x200 pixels, per unit incident flux
};

#endif /* _STAR_HXX_ */
