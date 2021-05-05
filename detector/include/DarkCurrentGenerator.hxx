/*
 * DarkCurrentGenerator.hxx
 *
 *  Created on: 11 Mar 2014
 *      Author: futyand
 */

#ifndef _DARK_CURRENT_GENERATOR_HXX_
#define _DARK_CURRENT_GENERATOR_HXX_

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int.hpp"
#include "boost/random/exponential_distribution.hpp"
#include "boost/random/bernoulli_distribution.hpp"
#include "boost/random/variate_generator.hpp"

#include "simulator/include/Module.hxx"

typedef boost::variate_generator<boost::mt19937,boost::uniform_int<int> > RANDOM_UNIFORM;
typedef boost::variate_generator<boost::mt19937,boost::exponential_distribution<float> > RANDOM_EXPONENTIAL;
typedef boost::variate_generator<boost::mt19937,boost::bernoulli_distribution<> > RANDOM_BERNOULLI;

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module generates a static dark frame (shot noise is added
 *  	   later by PhotonNoiseGenerator) with mean value depending on CCD
 *  	   temperature, and generates hot, warm and telegraphic pixels
 *
 *  The expected dark current depends on the CCD temperature and takes into
 *  account the total time during which pixels can accumulate dark current,
 *  including readout time.
 */

class DarkCurrentGenerator: public Module {
public:
	DarkCurrentGenerator();
	virtual ~DarkCurrentGenerator();

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:

	/// @brief Manually defined hot pixel
	struct ManualHotPixel {
		/** *************************************************************************
		 *  @brief Constructor for manual hot pixel
		 *
		 *  @param [in] ix  x position of the pixel
		 *  @param [in] iy  y position of the pixel
		 *  @param [in] rate  Dark current rate at 233K
		 *  @param [in] isTelegraphic  Boolean to indicate if the pixel is telegraphic
		 */
		ManualHotPixel(int ix, int iy, double rate, bool isTelegraphic) :
						m_xPixel(ix), m_yPixel(iy),m_rate(rate), m_isTelegraphic(isTelegraphic) {};
		int m_xPixel; ///< x position of the pixel
		int m_yPixel; ///< y position of the pixel
		double m_rate; ///< dark current rate at 233K
		bool m_isTelegraphic; ///< boolean to indicate if the pixel is telegraphic
	};

	/** *************************************************************************
	 *  @brief Read in the empirical dark frame data
	 */
	void readDarkFrame();

	/** *************************************************************************
	 *  @brief Define the list of telegraphic pixels, randomly generating their
	 *  	   positions, initial state and list of transitions,
	 *  	   and add them to the Data
	 *
	 *  @param [in] data  Pointer to the Data
	 */
	void setTelegraphicPixels(Data * data);

	/** *************************************************************************
	 *  @brief Generate a list of times at which a transition takes place for a
	 *  	   given telegraphic pixel, by drawing randomly from an exponential
	 *  	   distribution
	 *
	 *  @param [in] duration  Simulation total duration in seconds
	 */
	vector<float> telegraphicTransitions(double duration) const;

	/** *************************************************************************
	 *  @brief Generate dark electrons in the specified pixel, taking into
	 *  account whether or not the pixel is hot, warm active or inactive
	 *  telegraphic
	 *
	 *  @param [in] image  Pointer to the Image in memory
	 *  @param [in] imageToSmearUp  Pointer to the Image extended to cover the
	 *  							full frame in vertical direction, corresponding
	 *  							to the first second of the exposure, needed to
	 *  							generate upward smear trails
	 *  @param [in] imageToSmearDown  Pointer to the Image extended to cover the
	 *  							full frame in vertical direction, corresponding
	 *  							to the last second of the exposure, needed to
	 *  							generate downward smear trails
	 *  @param [in] ix  x position of pixel
	 *  @param [in] iy  y position of pixel
	 *  @param [in] temperatureFactor  Factor by which the dark current is
	 *  							   increased at the current temperature
	 *  							   w.r.t. 233K
	 *  @param [in] telegraphicPixels  List of telegraphic pixels
	 *  @param [in] timeStep  Current time step, needed to determine whether
	 *  					  or not a telegraphic pixel is active
	 *  @param [in] repetitionPeriod  Repetition period, needed to determine whether
	 *  					  or not a telegraphic pixel is active
	 *  @param [in] doChargeTransferEoL  Flag to indicate whether or not end of
	 *  						   		 life CTI is to be simulated
	 */
	void generateDarkElectrons(Image * image, Image * imageToSmearUp, Image * imageToSmearDown,
			int ix, int iy, double temperatureFactor,
			vector<Data::TelegraphicPixel> * telegraphicPixels, int timeStep, double repetitionPeriod,
			bool doChargeTransferEoL) const;

	/** *************************************************************************
	 *  @brief Returns the mean dark current for a normal pixel at the specified temperature
	 *
	 *  @param [in] temperature  CCD temperature for the current time step of the simulation
	 */
	double darkTemperatureDependence(double temperature) const;

	bool m_empiricalDarkFrame; ///< Flag to indicate whether or not to base dark current generation on an empirical dark frame
	string m_darkFrameFilename; ///< Filename for FITS file containing the empirical dark frame
	double m_darkFrame[Image::kXTotal][Image::kYTotal]; ///< Dark frame array
	double m_darkFrameScaleFactor; ///< Scale factor to be applied to the empirical dark frame

	double m_meanDarkCurrent233K; ///< User specified mean dark current at 233K
	double m_rowDownShiftTime; ///< Time to shift a row down during readout in seconds
	double m_serialShiftTime; ///< Time to shift the serial register one pixel to the left during readout in seconds
	double m_exposureTimeWithReadout[Image::kXTotal][Image::kYTotal]; ///< Exposure duration for each pixel taking into account the delay for that pixel to be read out

	bool m_doHotPixels; ///< Flag to indicate whether or not hot pixels will be generated
	int m_nHotPixels; ///< Number of hot pixels to be generated across the full frame including margins
	double m_hotPixelRelativeDarkCurrent; ///< Factor by which dark current is increased for hot pixels
	bool m_isHotPixel[Image::kXDim+2*Image::kNDarkCols][Image::kYDim+Image::kNDarkRows]; ///< Array of flags to indicate whether or not each pixel is hot

	bool m_doWarmPixels; ///< Flag to indicate whether or not warm pixels will be generated
	int m_nWarmPixels; ///< Number of warm pixels to be generated across the full frame including margins
	double m_warmPixelRelativeDarkCurrent; ///< Factor by which dark current is increased for warm pixels
	bool m_isWarmPixel[Image::kXDim+2*Image::kNDarkCols][Image::kYDim+Image::kNDarkRows]; ///< Array of flags to indicate whether or not each pixel is warm

	bool m_doTelegraphicPixels; ///< Flag to indicate whether or not telegraphic pixels will be generated
	int m_nTelegraphicPixels; ///< Number of telegraphic pixels to be generated across the full frame including margins
	double m_telegraphicPixelRelativeDarkCurrent; ///< Factor by which dark current is increased for active telegraphic pixels
	double m_telegraphicTimeConstant; ///< Time constant for telegraphic pixel transitions

	bool m_doManualHotPixels; ///< Flag to indicate whether or not to define hot/warm/telegraphic pixels manually
	vector<ManualHotPixel> m_manualHotPixels; ///< vector of manually defined hot/warm/telegraphic pixels

	RANDOM_EXPONENTIAL * m_exponentialRandomGenerator; ///< Exponential random number generator to generate telegraphic pixel transition times
	RANDOM_BERNOULLI * m_bernoulliRandomGenerator; ///< Bernoulli random number generator to assign initial state of telegraphic pixels
    RANDOM_UNIFORM * m_uniformRandomGeneratorX; ///< Uniform number generator for hot/warm/telegraphic pixel x positions
    RANDOM_UNIFORM * m_uniformRandomGeneratorY; ///< Uniform number generator for hot/warm/telegraphic pixel y positions
};

#endif /* _DARK_CURRENT_GENERATOR_HXX_ */
