/*
 * PSFGenerator.hxx
 *
 *  Created on: 10 Feb 2014
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This module is used to generate an image of the PSF at the
///		   appropriate position on the focal plane for each star in the
///		   field of view
///
/// Parameter input for ModuleParams:
///		1. FOVradius (double): Field of view radius in arcseconds
///		2. pointingRA (double): Pointing direction Right Ascension in arcseconds
///		3. pointingDec (double): Pointing direction declination in arcseconds
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

#ifndef _PSF_GENERATOR_HXX_
#define _PSF_GENERATOR_HXX_

#include "simulator/include/Module.hxx"

class PSFGenerator: public Module {
public:

	/** *************************************************************************
	 * Telescope temperature in cm^2. Originally 767.93cm^2, from Redbook
	 * section 4.4.1. Changed to 728cm^2 following communication from A.Fortier
	 * 17/09/2014, where the 40cm^2 difference is due to spider obscuration.
	 * Changed again 3/11/2014 to pi*(32/2)^2=804.24772 since latest throughput
	 * file already includes spider and M2 obscuration
	 * Changed again 19/02/2020 to pi*(30/2)^2=706.85835 to reflect new
	 * understanding of effective telescope aperture radius (32cm->30cm)
	 */
	static constexpr double kTelescopeArea = 706.85835;
	static const int kOversampledPsfSize = 2000; ///< Array size for oversampled PSF
	static const int kPsfSize = 200; ///< Array size for non-oversampled PSF
	static const int kNpsfTemp = 4; ///< Number of temperature bins for the PSF
	static const int kNpsfWavelength = 15; ///< Number of wavelength bins for the monochromatic PSF
	static const int kPsfWavelengthBinWidth = 50; ///< Width of wavelength bins for the monochromatic PSF in nm
	static const int kPsfFirstWavelength = 375; ///< Central value of first wavelength bin for the monochromatic PSF in nm
	static const int kPsfMargin = 25; ///< PSFs are not added to the image if the central position
	  	  	  	  	  	  	  	   	  ///< is more than this number of pixels outside the edge of the frame
								   	  ///< (intended to correspond to PSF radius)

	/// @brief Algorithm used to position the PSF on the FOV
	enum POSITION_METHOD{
		Oversampling, ///< 10 times spatial oversampling
		Interpolation, ///< Bilinear spatial interpolation
		SnapToGrid ///< PSF grid directly mapped to CCD grid
	};

	/// @brief Thermal map corresponding to the PSF
	enum THERMAL_MAP{
		fixed, ///< All telescope components at the nominal temperature of -10C
		cold, ///< Thermal map corresponding to the coldest orbit
		hot1, ///< Thermal map corresponding to the coldest point of the hottest orbit
		hot2, ///< Thermal map corresponding to the hottest point of the hottest orbit
		breathing ///< Linear interpolation between the hot1 and hot2 thermal maps according to telescope temperature
	};

	PSFGenerator() : Module("PSFGenerator",timeLoop) {};
	virtual ~PSFGenerator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Reads in the PSF (can be oversampled or not, monochromatic or not)
	 *  	   for the requested thermal map, or for the two hot thermal maps
	 *  	   for the breathing case. The resulting PSF is re-normalized to 1.
	 *
	 *  @param [in] oversampled  Flag to indicate whether or not the PSF should
	 *  						 read spatially oversampled
	 *  @param [in] subarrayDimensions	Subarray dimensions
	 */
	void readPSFs(bool oversampled, Data::SubarrayDimensions subarrayDimensions);

	/** *************************************************************************
	 *  @brief Generates wavelength weights according to the black body spectrum
	 *  	   defined by the effective temperature of the target star, together
	 *  	   with the QE and throughput vs wavelength.
	 *  	   The weights are used to construct the monochromatic PSF.
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void getWavelengthWeights(Data * data);

	/** *************************************************************************
	 *  @brief Performs a bilinear spatial interpolation the PSF to an array
	 *  	   shifted w.r.t. the pixel grid and normalizes the resultant PSF to 1
	 *
	 *  @param [in] jitteredPositionX  x position of the PSF taking into account
	 *  			jitter
	 *  @param [in] jitteredPositionY  y position of the PSF taking into account
	 *  			jitter
	 *  @param [in] temperatureInterpolatedPsf  PSF after interpolation to the
	 *  			current temperature (as opposed to spatial interpolation
	 *  			performed by this method) for the breathing case.
	 *  			Only used if the thermal map is set to breathing.
	 *  @param [out] interpolatedPSF  Output PSF after spatial interpolation to
	 *  							  the jittered position
	 */
	void interpolatePSF(double jitteredPositionX, double jitteredPositionY, double ** temperatureInterpolatedPsf,
						double ** interpolatedPSF) const;

	/** *************************************************************************
	 *  @brief Called by interpolatePSF. Returns the value of a pixel in the
	 *  	   target grid (CCD) given the values of the 4 nearest pixels in the
	 *  	   source grid (PSF) and the offset between the 2 grids in x and y
	 *
	 *  @param [in] q11  Value for nearest source grid pixel below and left of
	 *  				 the target pixel position
	 *  @param [in] q12  Value for nearest source grid pixel above and left of
	 *  				 the target pixel position
	 *  @param [in] q21  Value for nearest source grid pixel below and right of
	 *  				 the target pixel position
	 *  @param [in] q22  Value for nearest source grid pixel above and right of
	 *  				 the target pixel position
	 *  @param [in] xfrac  Offset of the target grid w.r.t. the source grid in the
	 *  				   x direction, expressed as a fraction of a pixel size
	 *  @param [in] yfrac  Offset of the target grid w.r.t. the source grid in the
	 *  				   y direction, expressed as a fraction of a pixel size
	 */
	double bilinearInterpolation(double q11, double q12, double q21, double q22, double xfrac, double yfrac) const;

	string m_psfFilename; ///< Filename for fits file used to define the PSF
	bool m_oversampleJitter; ///< Flag to indicate whether or not to perform time oversampling of PSF position according to jitter and field of view rotation
	POSITION_METHOD m_targetStarPositioning; ///< Method for positioning the target star PSF onto the pixel grid: Oversampling, Interpolation, or SnapToGrid
	POSITION_METHOD m_backgroundStarPositioning; ///< Method for positioning background star PSFs onto the pixel grid: Oversampling, Interpolation, or SnapToGrid
	bool m_convertFluxToPhotons; ///< Flag to indicate whether or not to output pixel values as number of photons or ergs
	bool m_monochromaticPSF; ///< Flag to indicate whether or not the PSF is constructed from PSFs binned in wavelength
	int m_monochromaticPSFwavelength; ///< Wavelength for monochromatic PSF. A weighted combination over wavelengths is performed if value is 0.
	THERMAL_MAP m_thermalMap; ///< Thermal map corresponding to the PSF: fixed, cold, hot1, hot2 or breathing
	double m_oversampledPsf[2000][2000][2]; ///< Array to hold the oversampled PSF. 3rd dimension is thermal map used for breathing case only.
	double m_psf[200][200][2]; ///< Array to hold the non-oversampled PSF. 3rd dimension is thermal map used for breathing case only.
	double m_weight[kNpsfWavelength]; ///< Weight assigned to each PSF wavelength bin for the monochromatic case

};

#endif /* _PSF_GENERATOR_HXX_ */
