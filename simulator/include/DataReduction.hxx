/*
 * DataReduction.hxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

#ifndef _DATA_REDUCTION_HXX_
#define _DATA_REDUCTION_HXX_

#include "boost/math/tools/stats.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/poisson_distribution.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include "SCI_COR_Lightcurve.hxx"
#include "SIM_ANA_Noisecurve.hxx"

#include "Photometry.hxx"
#include "Module.hxx"

typedef boost::variate_generator<boost::mt19937,boost::poisson_distribution<int,double> > RANDOM_POISSON;
typedef boost::variate_generator<boost::mt19937,boost::normal_distribution<double> > RANDOM_GAUSS;

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module performs a basic photometric extraction from the set
 *  	   of stacked images that have been generated.
 *
 *  The module can be configured to output normalized and un-normalized
 *  incident and extracted light curves, as well as a noise curve.
 *  Aperture photometry is performed using the Photometry class. The centre
 *  of the aperture can either be the truth PSF position (default), or the
 *  barycentre determined by photometry. Options include background
 *  subtraction based on the median value within an anulus, and flat field
 *  correction.
 */

class DataReduction: public Module {
public:
	DataReduction() : Module("DataReduction",end), m_imageDirectory("null"),
			m_writeIncidentLightCurve(true),m_extractLightCurve(true), m_generateNoiseCurve(true),
			m_truthBarycentre(true),m_subtractFlatField(false){}
	virtual ~DataReduction() {delete m_photometry; delete m_poissonNoiseGenerator; delete m_gaussianNoiseGenerator;};

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const {};
	void doEnd(Data *data) const;

private:

	/// @brief Struct to hold the extracted flux time series
	struct Flux {
		/** *************************************************************************
		 *  @brief Constructor
		 *
		 *  @param [in] flux  Value of the flux from photometric extraction
		 *  @param [in] midTime  UTC mid time of the image from which the flux was extracted
		 *  @param [in] status  Flag with value 0 to 3 depending on the jitter AOCS and science validity flags
		 *  @param [in] rollAngle  Roll angle
		 *  @param [in] centroidX  X location of centroid w.r.t. left edge of exposed part of CCD
		 *  @param [in] centroidY  Y location of centroid w.r.t. bottom edge of exposed part of CCD
		 */
		Flux(double flux, UTC midTime, int status, double rollAngle, double centroidX, double centroidY) : m_flux(flux), m_midTime(midTime),
				m_status(status), m_rollAngle(rollAngle), m_centroidX(centroidX), m_centroidY(centroidY) {};
		double m_flux; ///< Value of the flux from photometric extraction
		UTC m_midTime; ///< UTC mid time of the image from which the flux was extracted
		int m_status; ///< Flag with value 0 to 3 depending on the jitter AOCS and science validity flags
		double m_rollAngle; ///< Roll angle
		double m_centroidX; ///< X location of centroid w.r.t. left edge of exposed part of CCD
		double m_centroidY; ///< Y location of centroid w.r.t. bottom edge of exposed part of CCD
	};

	/** *************************************************************************
	 *  @brief Perform photometric extraction: reads in the set of image cubes
	 *  	   and for each cube calls extractFluxesFromCube
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @return Flux vector containing the extracted flux
	 *  		for each time step of the simulation, together with
	 *  		jitter AOCS and science validity flags
	 */
	vector<Flux> extractFluxVector(Data * data) const;

	/** *************************************************************************
	 *  @brief Perform photometric extraction for images within the current
	 *  	   image cube and assigns jitter AOCS and science validity flags.
	 *  	   Calls Photometry::extractFlux for each image
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] filename  Filename for the image cube
	 *  @param [in] filename  Filename for the image cube
	 *  @param [in,out] nFiles  Counter to keep track of the number of files
	 *  						that have been read
	 *  @param [in,out] nImages  Counter to keep track of the total number of
	 *  						images that have been processed
	 *  @param [in,out] extractedFlux  Flux vector containing the extracted flux
	 *  							   for each time step of the simulation,
	 *  							   together with jitter AOCS and science
	 *  							   validity flags
	 */
	template<class T> void extractFluxesFromCube(Data* data, string filename, int & nFiles, int & nImages, vector<Flux> & extractedFlux) const;

	/** *************************************************************************
	 *  @brief For the case where data reduction is run on a pre-existing set of
	 *  	   images from an earlier simulation, use the total number of images
	 *  	   that have been processed together with header keywords from the
	 *  	   image cubes to define the time configuration, so that time
	 *  	   related header keywords can be set for the output light curves.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in,out] nImages  Total number of images that have been processed
	 */
	template <class T> void setTimeConfigurationFromImageCube(Data* data, unsigned nImages) const;

	/** *************************************************************************
	 *  @brief Write the incident light curve, including normalized curve
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void writeIncidentLightCurve(Data * data) const;

	/** *************************************************************************
	 *  @brief Write the output light curves, including normalized and O-C curves
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] extractedFlux  Flux vector containing the extracted flux
	 *  						   for each time step of the simulation,
	 *  						   together with jitter AOCS and science
	 *  						   validity flags
	 */
	void writeOutputLightCurves(Data * data, vector<Flux> & extractedFlux) const;

	/** *************************************************************************
	 *  @brief Set the header keywords for a light curve, calling either
	 *  	   setLightCurveHeaderKeywordsFromImageCube or
	 *  	   setLightCurveHeaderKeywordsWithoutImageCube
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] lightCurve  Pointer to the light curve FITS table
	 *  @param [in] lastImageMidTime  Mid time of the last image of the simulation
	 *  @param [in] photometry  Boolean to indicate whether or not LC is from photometric extraction
	 */
	void setLightCurveHeaderKeywords(Data * data, SciCorLightcurve * lightCurve, UTC lastImageMidTime, bool photometry=false) const;

	/** *************************************************************************
	 *  @brief Set the header keywords for a light curve by copying header
	 *  	   keywords from the last image cube to be processed
	 *
	 *  @param [in] lightCurve  Pointer to the light curve FITS table
	 *  @param [in] imageCube  Pointer to the last FITS image cube to be read
	 */
	template <class T> void setLightCurveHeaderKeywordsFromImageCube(SciCorLightcurve * lightCurve, T * imageCube) const;

	/** *************************************************************************
	 *  @brief For the case where only an incident light curve is generated
	 *  	   such that there are no image cubes, set the header keywords
	 *  	   explicitly
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] lightCurve  Pointer to the light curve FITS table
	 */
	void setLightCurveHeaderKeywordsWithoutImageCube(const Data * data, SciCorLightcurve * lightCurve) const;

	/** *************************************************************************
	 *  @brief Write a row of the light curve FITS table
	 *
	 *  @param [in] lightCurve  Pointer to the light curve FITS table
	 *  @param [in] flux  Value of the flux for the current time step
	 *  @param [in] midTime  UTC mid time of the current time step
	 *  @param [in] status  Status flag with value 0 to 3, depending on the
	 *  					jitter AOCS and science validity flags
	 *  @param [in] rollAngle  Roll angle
	 *  @param [in] centroidX  X location of centroid w.r.t. left edge of exposed part of CCD
	 *  @param [in] centroidY  Y location of centroid w.r.t. bottom edge of exposed part of CCD
	 *  @param [in] locationX  X intended location of target w.r.t. left edge of exposed part of CCD
	 *  @param [in] locationY  Y intended location of target w.r.t. bottom edge of exposed part of CCD
	 */
	void setLightCurveData(SciCorLightcurve * lightCurve, double flux, UTC midTime, int status, double rollAngle=-1.,
			double centroidX=0., double centroidY=0., double locationX=0., double locationY=0.) const;

	/** *************************************************************************
	 *  @brief Generate a noise curve from the data in the extracted light curve
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] flux  Vector of normalized transit subtracted flux values
	 *  				  for each time step
	 *  @param [in] lastImageMidTime  Mid time of the last image of the simulation
	 *  @param [in] filename  Name of the SIM_ANA_NoiseCurve FITS file
	 */
	void generateNoiseCurve(const Data * data, const vector<double> flux, UTC lastImageMidTime, string filename) const;

	/** *************************************************************************
	 *  @brief Set the header keywords for the noise curve by copying header
	 *  	   keywords from the last image cube to be processed
	 *
	 *  @param [in] noiseCurve  Pointer to the noise curve FITS table
	 *  @param [in] imageCube  Pointer to the last FITS image cube to be read
	 */
	template <class T> void setNoiseCurveHeaderKeywordsFromImageCube(SimAnaNoisecurve & noiseCurve, T * imageCube) const;

	/** *************************************************************************
	 *  @brief Appends the extracted flux from the first image to an ascii file.
	 *	This can be used together with the extractedFlux.sh script to generate an
	 *	array of extracted flux values for different magnitudes and spectral types.
	 *
	 *  @param [in] data  Pointer to the Data
	 *  @param [in] flux   extracted flux in ADU
	 */
	void writeExtractedFlux(Data * data, double flux) const;

	string m_imageDirectory; ///< Name of directory containing input images, for case where data reduction is performed on a pre-existing set of images
	string m_outputDirectory; ///< Name of output directory

	bool m_writeIncidentLightCurve; ///< Flag to indicate whether or not incident light curves will be output
	bool m_extractLightCurve; ///< Flag to indicate whether or not extracted light curves will be output
	bool m_generateNoiseCurve; ///< Flag to indicate whether or not a noise curve will be output
	bool m_truthBarycentre; ///< Flag to indicate whether to use the truth barycentre or the barycentre determined by photometry
	bool m_subtractFlatField; ///< Flag to indicate whether or not to perform flat field correction

	Photometry * m_photometry; ///< Pointer to an instance of the Photometry class, used to perform the photometric extraction

	RANDOM_POISSON * m_poissonNoiseGenerator; ///< Poisson random number generator to generate the ideal light curve
	RANDOM_GAUSS * m_gaussianNoiseGenerator; ///< Gaussian random number generator for read noise
};

#endif /* _DATA_REDUCTION_HXX_ */
