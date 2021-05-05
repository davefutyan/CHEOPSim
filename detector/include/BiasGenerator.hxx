/*
 * BiasGenerator.hxx
 *
 *  Created on: 11 Mar 2014
 *      Author: futyand
 */

#ifndef _BIAS_GENERATOR_HXX_
#define _BIAS_GENERATOR_HXX_

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"
#include "boost/random/poisson_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include "simulator/include/Module.hxx"

typedef boost::variate_generator<boost::mt19937,boost::normal_distribution<double> > RANDOM_GAUSS;
typedef boost::variate_generator<boost::mt19937,boost::poisson_distribution<int,double> > RANDOM_POISSON;

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief The module applies temperature dependent gain, and generates a bias
 *  	   frame with readout noise by drawing randomly for each
 *  	   pixel from a Gaussian distribution
 */

class BiasGenerator: public Module {
public:

	BiasGenerator() : Module("BiasGenerator",timeLoop) {};
	virtual ~BiasGenerator() {};

	void initialize(const ModuleParams & params);
	void doBegin(Data *data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

	/** *************************************************************************
	 *  @brief Apply the inverse of the CCD non-linearity correction modeled using a cubic spline
	 *
	 *  @param [in,out] nElectrons  Number of electrons before [in] and after [out] inverse correction
	 */
	void applyCcdNonLinearity(double & nElectrons) const;

	/** *************************************************************************
	 *  @brief Apply the CCD non-linearity correction modeled using a cubic spline
	 *  	   Used for unit testing only
	 *
	 *  @param [in,out] nElectrons  Number of electrons before [in] and after [out] correction
	 */
	void correctCcdNonLinearity(double & nElectrons) const;

private:

	/** *************************************************************************
	 *  @brief Read in the empirical bias frame data
	 *
	 *  @param [in] serialReadRate  serial read rate in kHz
	 *  @param [in] redundantHardware  boolean to indicate whether or not redundant readout hardware is used
	 *  @param [in] ccdTemperature  temperature of the CCD
	 */
	void readBiasFrame(double serialReadRate, bool redundantHardware, double ccdTemperature);

	/** *************************************************************************
	 *  @brief Convert from number of electrons to ADU for the specified pixel,
	 *  	   by dividing by the electronic gain and
	 *  	   adding random read noise with bias offset
	 *
	 *  @param [in] image  Pointer to the image
	 *  @param [in] ix  x index of the pixel
	 *  @param [in] iy  y index of the pixel
	 *  @param [in] temperature  Current CCD temperature, used to determine the gain
	 *  @param [in,out] pixelADU_max  Parameter to keep track of the maximum ADU of the image (currently only for diagnostic printout)
	 *  @param [in] fullFrame  Flag to indicate whether or not the current image is full frame (readout noise can differ)
	 */
	void generateBiasNoise(Image * image, int ix, int iy, double temperature, double & pixelADU_max, bool fullFrame) const;

	/** *************************************************************************
	 *  @brief Read in the gain formula
	 *
	 *  @param [in] gainCorrectionFilename  Filename for the REF_APP_GainCorrection
	 *  									fits file defining the formula
	 */
	void readGainCorrection(string gainCorrectionFilename);

	/** *************************************************************************
	 *  @brief Read in the CCD non-linearity coefficients
	 *
	 *  @param [in] ccdNonLinearityFilename  Filename for the REF_APP_CCDLinearisation
	 *  									  fits file containing the coefficients
	 */
	template <class T> void readCcdNonLinearity(string ccdNonLinearityFilename);

	/// @brief Struct to store the set of values used to define each line of the gain correction formula
	struct GainFormulaLine {
		/** *************************************************************************
		 *  @brief Constructor for GainFormulaLine
		 *
		 *  @param [in] constant  Constant term for the current line of the gain correction formula
		 *  @param [in] vodExponent  Exponent of the deviation of the VOD voltage from its nominal value for the current line of the gain correction formula
		 *  @param [in] vrdExponent  Exponent of the deviation of the VRD voltage from its nominal value for the current line of the gain correction formula
		 *  @param [in] vogExponent  Exponent of the deviation of the VOG voltage from its nominal value for the current line of the gain correction formula
		 *  @param [in] vssExponent  Exponent of the deviation of the VSS voltage from its nominal value for the current line of the gain correction formula
		 *  @param [in] tempExponent  Exponent of the deviation of the CCD temperature from its nominal value for the current line of the gain correction formula
		 */
		GainFormulaLine(double constant, double vodExponent, double vrdExponent, double vogExponent, double vssExponent, double tempExponent) :
			m_constant(constant), m_vodExponent(vodExponent), m_vrdExponent(vrdExponent), m_vogExponent(vogExponent), m_vssExponent(vssExponent), m_tempExponent(tempExponent) {};
		double m_constant; ///< Constant term for the current line of the gain correction formula
		double m_vodExponent; ///< Exponent of the deviation of the VOD voltage from its nominal value for the current line of the gain correction formula
		double m_vrdExponent; ///< Exponent of the deviation of the VRD voltage from its nominal value for the current line of the gain correction formula
		double m_vogExponent; ///< Exponent of the deviation of the VOG voltage from its nominal value for the current line of the gain correction formula
		double m_vssExponent; ///< Exponent of the deviation of the VSS voltage from its nominal value for the current line of the gain correction formula
		double m_tempExponent; ///< Exponent of the deviation of the CCD temperature from its nominal value for the current line of the gain correction formula
	};

	bool m_empiricalBiasFrame; ///< Flag to indicate whether to read an empirical bias frame rather than a uniform bias offset for all pixels
	string m_biasFrameFilename; ///< Filename for FITS file containing the empirical bias frame
	double m_bias[Image::kXTotal][Image::kYTotal]; ///< Bias offset array
	double m_RON[Image::kXTotal][Image::kYTotal]; ///< Read Out Noise array
	double m_biasMean; ///< Bias offset
	double m_biasWidth; ///< Readout noise Gaussian width
	bool m_applyAnalogElectronicsStability; ///< Boolean to indicate whether or not to include analog electronics stability
	double m_nominalGain; ///< Nominal gain
	vector<GainFormulaLine> m_gainFormula; ///< Formula for gain correction as a function of bias voltages (analog electronics stability)
	bool m_doCcdNonLinearity; ///< Boolean to indicate whether or not to include gain non-linearity
	string m_ccdLinearityFilename; ///< Filename for the CCD non-linearity parameterization
	string m_ccdLinearityExtension; ///< FITS extension for the CCD non-linearity parameterization
	bool m_aduNonLinearity; ///< Boolean to indicate whether the non-linearity parameterization is in terms of ADU or electrons
	double m_ccdNonlinearityBoundaries[29]; ///< Boundaries between spline intervals used to model the CCD non-linearity
	double m_ccdNonlinearityCoeffs[4][28]; ///< Coefficients of 3rd order polynomials used to model the CCD non-linearity in each spline interval
	unsigned m_nIntervals; ///< Number of spline intervals
	double m_ccdNonlinearityExtrapSlope; ///< Slope for linear extrapolation of CCD non-linearity model for ADU>70000
	double m_ccdNonlinearityExtrapIntercept; ///< Intercept for linear extrapolation of CCD non-linearity model for ADU>70000
	RANDOM_GAUSS * m_gaussianNoiseGenerator; ///< Gaussian random number generator for read noise
	RANDOM_GAUSS * m_gaussianNoiseGenerator_fullFrame; ///< Gaussian random number generator for read noise for the full frame image
	//bool m_applyAnalogChainRandomNoise; ///< Boolean to indicate whether or not to apply analog chain random noise
	//RANDOM_POISSON * m_poissonNoiseGenerator; ///< Poisson number generator for analog chain random noise
};

#endif /* _BIAS_GENERATOR_HXX_ */
