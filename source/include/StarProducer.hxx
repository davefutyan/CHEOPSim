/*
 * StarProducer.hxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to define the field of view and to add to
///        it the list of Stars
///
/// Parameter input for ModuleParams:
///		1. FOVradius (double): Field of view radius in arcseconds
///		4. minMagnitude (int): Minimum magnitude for stars in the field of view
///		5. starFilename (string): ascii file containing comma separated values
///							      for RA,dec,mag,specType for each star to be included
///
/// The field of view information is stored in an instance of SkyFieldOfView.
/// The information is passed into the constructor via an instance of
/// ModuleParams.
///
/// Stars are read in from a file names stars.txt located in Simulator/conf
/// This is an ascii file with one line for each star.  Each line consists
/// of three comma separated doubles followed by a string, specifying
/// respectively for each star: right ascension, declination, magnitude,
/// spectral type.
///
////////////////////////////////////////////////////////////////////////

#ifndef _STAR_PRODUCER_HXX_
#define _STAR_PRODUCER_HXX_

#include <iostream>
#include <fstream>
#include <string>

#include "simulator/include/Module.hxx"
#include "data/include/Star.hxx"
#include "FluxConverter.hxx"

class StarProducer : public Module {
public:

	static constexpr double kNPhoVega = 3657535819; ///< Number of photons/s incident on the telescope Aperture from Vega,
													///< obtained by running StarProducer with m_VegaSpectrum = true in
													///< StarProducer.cxx and magnitude 0.035 in mystars.txt
													///< Updated 19/02/2020, scaling previous value by (30/32)^2 to reflect new
													///< understanding of effective telescope aperture radius (32cm->30cm)

	/// @brief Enum to select density of background star field
	enum FIELD_CROWDING{low,medium,extreme};

	StarProducer() : Module("StarProducer",begin) {};
	virtual ~StarProducer() {delete m_fluxConv;};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const {};

	/** *************************************************************************
	 *  @brief Converts from hours:mins:seconds or degrees:mins:seconds to
	 *  	   arcseconds
	 *
	 *  @param [in] hoursMinsSecs  right ascension or declination in
	 *  						   hours:mins:seconds or degrees:mins:seconds
	 *  @param [in] isRA  Boolean to indicate if input is right ascension
	 *  @return value of the input converted to arcseconds
	 */
	static double convertToArcseconds(string hoursMinsSecs, bool isRA);

	/** *************************************************************************
	 *  @brief Uses linear interpolation to transform an array of values sampling
	 *  	   a function of variable x to a new array with different x sampling
	 *
	 *  @param [in] x_ref  Input x array (e.g. wavelength with binning of input file)
	 *  @param [in] y_ref  Input y array (e.g. flux with binning of input file)
	 *  @param [in] size_ref  Number of bins for input arrays
	 *  @param [in] x  Output x array (e.g. wavelength with output binning)
	 *  @param [in] y  Output y array (e.g. flux with output binning)
	 *  @param [in] size  Number of bins for output arrays
	 */
	static void interpolateFlux(double * x_ref, double * y_ref, unsigned size_ref, double * x, double * y, unsigned size);

	/** *************************************************************************
	 *  @brief Returns the difference between Gaia magnitude and V-band
	 *  	   magnitude for a star with the specified effective temperature
	 *
	 *  @param [in] effectiveTemperature  Black body effective temperature
	 */
	static double GaiaMinusVbandMagnitude(double effectiveTemperature);

	/** *************************************************************************
	 *  @brief Returns the flux spectrum of Vega in ergs/s/cm2/A, in 1141 bins,
	 *  taken from http://dls.physics.ucdavis.edu/calib/Vega.sed
	 *
	 *  @param [inout] wavelength	Units: Angstroms
	 *  @param [inout] flux   Units: ergs/s/cm2/A
	 */
	static void readVegaSpectrum(double (&wavelength)[1141], double (&flux)[1141] );

private:

	/** *************************************************************************
	 *  @brief Reads in the list of stars from the input ascii file
	 *
	 *  @param [in] fov  Pointer to the field of view
	 *  @param [in] nTimeSteps  Number of time steps in the simulation
	 *  						(for initialization of the StarData time series)
	 *  @param [in] wavelengthDependence  Pointer to instance of
	 *  								  WavelengthDependence class
	 *  @param [in] mpsVisits  Pointer to MPS_PRE_Visits instance
	 */
	void readStars(SkyFieldOfView *fov, const WavelengthDependence * wavelengthDependence, unsigned nTimeSteps, MpsPreVisits * mpsVisits) const;

	/** *************************************************************************
	 *  @brief Reads in the list of stars from the input star catalogue file
	 *
	 *  @param [in] fov  Pointer to the field of view
	 *  @param [in] nTimeSteps  Number of time steps in the simulation
	 *  						(for initialization of the StarData time series)
	 *  @param [in] wavelengthDependence  Pointer to instance of
	 *  								  WavelengthDependence class
	 *  @param [in] mpsVisits  Pointer to MPS_PRE_Visits instance
	 */
	template <class T> void readStarCatalogue(SkyFieldOfView *fov, const WavelengthDependence * wavelengthDependence, unsigned nTimeSteps, MpsPreVisits * mpsVisits) const;

	/** *************************************************************************
	 *  @brief Reads in the list of background stars from the input file
	 *  	   chosen according to the value of m_fieldCrowding
	 *
	 *  @param [in] fov  Pointer to the field of view
	 *  @param [in] nTimeSteps  Number of time steps in the simulation
	 *  						(for initialization of the StarData time series)
	 *  @param [in] wavelengthDependence  Pointer to instance of
	 *  								  WavelengthDependence class
	 */
	void readBackgroundStars(SkyFieldOfView *fov, const WavelengthDependence * wavelengthDependence, unsigned nTimeSteps) const;

	/** *************************************************************************
	 *  @brief Evaluates the flux of the Star passed as argument, given the
	 *  	   spectral type and magnitude, using the flux distribution
	 *  	   of Vega together with the V-band transmission curve to set the
	 *  	   absolute scale
	 *
	 *  @param [in] star  Pointer to the Star
	 *  @param [in] wavelengthDependence  Pointer to instance of
	 *  								  WavelengthDependence class
	 *  @param [in] target  Flag to indicate whether or not the Star is the target
	 */
	void evaluateFlux(Star * star, const WavelengthDependence * wavelengthDependence, bool target=false) const;

	/** *************************************************************************
	 *  @brief Writes the list of stars within the FOV to a file with
	 *  	   data structure EXT_PRE_StarCatalogue or EXT_DRFT_StarCatalogue
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] visitCounter  Allows to set visit counter should be 0 for T=EXT_DRFT_StarCatalogue
	 */
	template <class T> void writeStars(Data * data, unsigned visitCounter) const;

	double m_FOVradius; ///< Radius of the field of view
	double m_minMagnitude; ///< Faintest magnitude for stars to be included in the simulation
	std::string m_starFilename; ///< Name of file containing list of stars
	bool m_gaiaBand; ///< Flag to indicate whether or not to use the Gaia band rather than the V-band for stellar flux normalization
	double m_targetMagnitude; ///< Target star magnitude in either V-band or Gaia band according to m_gaiaBand
	bool m_overrideTargetMagnitude; ///< Flag to indicate if m_targetMagnitude should override the magnitude derived from the StarCatalogue file
	bool m_useTargetStarForPointing; ///< Flag to indicate whether or not to match the pointing direction to the target star
	bool m_doBackgroundStars; ///< Flag to indicate whether or not to generate background stars
	FIELD_CROWDING m_fieldCrowding; ///< Density of background star field

	bool m_VegaSpectrum; ///<Flag to indicate whether to use Vega spectrum rather than black body

	double m_wavelength[Star::kNWavelength*10]; ///< Array of wavelengths
	double m_flux_Vega[Star::kNWavelength*10]; ///< Flux of Vega for each wavelength bin
	double m_wavelength_fine[WavelengthDependence::kNWavelength]; ///< Array of wavelengths in 0.5nm steps
	double m_flux_Vega_fine[WavelengthDependence::kNWavelength]; ///< Flux of Vega for in 0.5nm steps
	double m_transmission_refBand[Star::kNWavelength*10]; ///< V-band or Gaia band transmission for each wavelength bin

	FluxConverter * m_fluxConv; ///< Pointer to FluxConverter class providing conversion from V-mag or Gaia mag to CHEOPS mag

};

inline void StarProducer::interpolateFlux(double* x_ref, double* y_ref, unsigned size_ref, double* x, double* y, unsigned size) {

	for (unsigned i=0; i<size; i++) {
		for (unsigned j=1; j<size_ref; j++) {
			if (x_ref[j]>x[i]) {
				double mu = (x[i]-x_ref[j-1])/(x_ref[j]-x_ref[j-1]);
				y[i] = y_ref[j-1]*(1.-mu)+y_ref[j]*mu;
				break;
			}
		}
	}

}

inline void StarProducer::readVegaSpectrum(double (&wavelength)[1141], double (&flux)[1141] ) {

	//Read in Vega spectrum from http://dls.physics.ucdavis.edu/calib/Vega.sed   Units: ergs/s/cm2/A
	ifstream vega_in(string(getenv("CHEOPS_SW"))+"/resources/flux_vega.txt", ios::in);
	if (!vega_in.good()) throw runtime_error("Error opening Vega flux file");
	string line;
	int i=0;
	double integratedFlux_Vega_fromFile=0.;
	while (getline(vega_in, line)) {
		if(line[0] != '#') {
			stringstream(line) >> wavelength[i] >> flux[i];
			if (i!=0) integratedFlux_Vega_fromFile += flux[i]*(wavelength[i]-wavelength[i-1]);
			//cout << "wavelength[" << i << "]=" << wavelength[i] << ", flux[" << i << "]=" << flux[i] << endl;
			i++;
		}
	}
	vega_in.close();

	//  cout << integratedFlux_Vega_fromFile << endl;
	//	ofstream file_out;
	//	file_out.open("fluxVsWavelength_Vega_fromFile.csv");
	//	for (int i=0;i<1141;i++) file_out << wavelength[i] << "," << flux[i] << endl;
	//   file_out.close();

}

#endif /* _STAR_PRODUCER_HXX_ */
