/*
 * FlatFieldGenerator.hxx
 *
 *  Created on: 4 Mar 2014
 *      Author: futyand
 */

#ifndef _FLAT_FIELD_GENERATOR_HXX_
#define _FLAT_FIELD_GENERATOR_HXX_

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/binomial_distribution.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/variate_generator.hpp"
#include "boost/math/tools/stats.hpp"

#include "CreateFitsFile.hxx"

#include "simulator/include/Module.hxx"

typedef boost::variate_generator<boost::mt19937,boost::uniform_int<int> > RANDOM_UNIFORM;
typedef boost::variate_generator<boost::mt19937,boost::normal_distribution<double> > RANDOM_GAUSS;

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module applies a flat field by multiplying the number of
 *  	   electrons per pixel by a flat field array.
 *
 *  The flat field can be defined either according to a Gaussian
 *  distribution, or using empirical flat field frames. In the case of an
 *  empirical flat field, the flat field array is constructed by performing
 *  a weighted sum of flat fields measured at various wavelengths, weighted
 *  according to the product of the stellar black body spectrum, the optical
 *  throughput and the quantum efficiency.
 */

class FlatFieldGenerator: public Module {
public:

	FlatFieldGenerator();
	virtual ~FlatFieldGenerator();

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/** *************************************************************************
	 *  @brief Read in the empirical flat field data
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void readFlatField(Data* data);

	/** *************************************************************************
	 *  @brief Randomly generate a Gaussian flat field
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void setGaussianFlatField(Data* data);

	/** *************************************************************************
	 *  @brief Generate the flat field image, with an option to add
	 *  a random Gaussian component as an uncertainty
	 *
	 *  @param [inout] flatFieldStats  mean and standard deviation of flat field
	 *  @param [in] withRandomNoise  Boolean to indicate whether or not to apply
	 *  							 a random Gaussian component in order to
	 *  							 define an uncertainty for the flat field
	 */
	Image * generateFlatField(boost::math::tools::stats<double> & flatFieldStats, bool withRandomNoise=false);

	/** *************************************************************************
	 *  @brief For each pixel of the image, multiply the number of electrons by
	 *  	   the flat field value for the pixel
	 *
	 *  @param [in] image  Pointer to the image
	 *  @param [in] flatField  Pointer to the combined flat field image
	 */
	void applyFlatField(Image * image, const Image * flatField) const;

	/** *************************************************************************
	 *  @brief Add each dead pixel to the truth FITS table
	 *
	 *  @param [in] image  Pointer to the image
	 *  @param [in] flatField  Pointer to the flat field image
	 */
	void storeDeadPixelTruthData(Image * image, Image * flatField) const;

	/** *************************************************************************
	 *  @brief Randomly generate dead pixel positions and write them out to
	 *  	   the REF_APP_BadPixel FITS table
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void setDeadPixels(Data * data);

	/** *************************************************************************
	 *  @brief Writes the dead pixels present in the empirical flat field
	 *  	   (sensitivity < 80% of nominal) to the REF_APP_BadPixel FITS table
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void setEmpiricalDeadPixels(Data * data);

	/** *************************************************************************
	 *  @brief Write the flat field to the SIM_TRU_FlatField data structure
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void writeTruthFlatField(Data * data) const;

	/** *************************************************************************
	 *  @brief Write the flat field to the REF_APP_WhiteFlatField data structure
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void writeWhiteFlatField(Data * data) const;

	bool m_applyFlatField; ///< Flag to indicate whether or not to apply the flat field
	double m_flatFieldTeffOffset; ///< Offset to the effective temperature in order that the flat field applied in CHEOPSim is not identical to that used for flat field correction in data reduction
	double m_flatFieldSmearSigma; ///< Standard deviation of the Gaussian distribution used to smear the flat field in order that the flat field applied in CHEOPSim is not identical to that used for flat field correction in data reduction
	double m_flatField[Image::kXDim][Image::kYDim]; ///< Flat field array
	double m_flatFieldTeff; ///< Effective temperature of target

	string m_flatFieldFilename; ///< Filename for FITS file containing the empirical flat field
	string m_flatFieldToSubtractFilename; ///< Filename for FITS file containing the empirical flat field to be used for flat field correction in photometry
	double m_flatFieldScaleFactor; ///< Scale factor to be applied to the empirical flat field

	bool m_gaussianFlatField; ///< Flag to indicate whether or not the flat field is Gaussian
	double m_flatFieldSigma; ///< standard deviation for Gaussian flat field

	bool m_writeTruthFlatField; ///< Flag to indicate whether or not to write out the truth flat field

	bool m_doDeadPixels; ///< Flag to indicate whether or not to simulate dead pixels
	unsigned m_deadPixelPositionSeed; ///< Seed for random number generation for dead pixel seed positions
	int m_nDeadPixels; ///< Number of dead pixels to be generated
	double m_deadPixelRelativeQE; ///< QE for dead pixels relative to that of normal pixels
	double m_deadPixelQE[Image::kXDim][Image::kYDim]; ///< Map of dead pixels containing their QE values relative to normal pixels

	RANDOM_GAUSS * m_gaussianNoiseGenerator; ///< Random number generator for Gaussian flat field
	RANDOM_GAUSS * m_gaussianNoiseGenerator2; ///< Random number generator to define Gaussian uncertainty for flat field correction
};

#endif /* _FLAT_FIELD_GENERATOR_HXX_ */
