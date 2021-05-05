/*
 * TimeConfiguration.hxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Container class for time configuration information
///
/// An instance of this class is contained in the Data class
///
////////////////////////////////////////////////////////////////////////

#ifndef _TIME_CONFIGURATION_HXX_
#define _TIME_CONFIGURATION_HXX_

#include "boost/date_time/posix_time/posix_time.hpp"
#include "Utc.hxx"
#include "DeltaTime.hxx"

class TimeConfiguration {
public:

	/** *************************************************************************
	 *  @brief Constructor
	 *
	 *  @param [in] startTime  Start time of the simulation
	 *  @param [in] numberOfStackedImages  Total number of stacked images
	 *  @param [in] exposuresPerStack  Number of exposures per stacked image
	 *  @param [in] exposureTimeAsDouble  Exposure duration
	 *  @param [in] maxImagesPerCube  Maximum number of images per image cube
	 *  @param [in] doFullFrame  Boolean to indicate whether or not the sub-array images are preceded by a full frame image
	 *  @param [in] repetitionPeriod  Repetition period (optional, required only for a non-default repetition period)
	 */
	TimeConfiguration(boost::posix_time::ptime startTime, unsigned numberOfStackedImages,
					  unsigned exposuresPerStack, double exposureTimeAsDouble, unsigned maxImagesPerCube, bool doFullFrame, double repetitionPeriod=0.);
	TimeConfiguration() {};
	virtual ~TimeConfiguration() {};

	/// @brief returns the start time of the simulation. This is the start time for the sub-array
	///		   sequence, which is delayed by 30s + full frame exposure duration + 5s,
	///		   relative to the start of the visit (the delay is 30s if there is no full frame image).
	///		   To get the start time of the visit, use getVisitStartTime().
	boost::posix_time::ptime getStartTime() const {return m_startTime;}

	/// @brief returns the start time of the simulation in UTC format. This is the start time for the sub-array
	///		   sequence, which is delayed by 30s + full frame exposure duration + 5s,
	///		   relative to the start of the visit (the delay is 30s if there is no full frame image).
	///		   To get the start time of the visit, use getVisitStartTime().
	UTC getStartTimeUTC() const {return to_iso_extended_string(m_startTime);}

	/// @brief returns the start time of the visit in UTC format. The sub-array sequence is delayed by 30s + full
	///        frame exposure duration + 5s (with full frame image) or 30s (without full frame image), relative to this.
	UTC getVisitStartTimeUTC() const {return to_iso_extended_string(m_visitStartTime);}

	/// @brief Returns the interval between the start times of consecutive exposures in seconds.
	///        For exposure durations >1s, the value is always an integer in order to sync with jitter which has 1s resolution.
	///        In reality, for exposure durations >1s, the repetition period should match the exposure duration.
	///        However exposure durations longer than 1s are expected to also always be integers.
	///        For exposure durations <1s, the value is equal to the exposure duration plus 1s
	double getRepetitionPeriod() const {return m_repetitionPeriod;}

	/// @brief Returns the number of times the jitter will be sampled per exposure.
	///        For exposure durations >1s, the value is equal to the default repetition period (exposure duration rounded up to the nearest integer).
	///        For exposure durations <1s, the value is 1.
	unsigned getJitterSamplingsPerExposure() const {return static_cast<unsigned>(ceil(m_exposureTimeAsDouble));}

	/// @brief Returns the exposure duration in seconds.
	double getExposureTimeAsDouble() const {return m_exposureTimeAsDouble;}

	/// @brief returns the number of exposures per stacked image
	unsigned getExposuresPerStack() const {return m_exposuresPerStack;}

	/// @brief returns the total number of stacked images
	unsigned getNumberOfStackedImages() const {return m_numberOfStackedImages;}

	/// @brief returns the number of stacked images for which the payload is in the loop
	unsigned getNumberOfStackedImagesPITL() const {return m_numberOfStackedImages_PITL;}

	/// @brief returns the total number of time steps (unstacked images) in the simulation
	unsigned getNumberOfTimeSteps() const {return m_numberOfStackedImages*m_exposuresPerStack;}

	/// @brief Sets the number of time steps (repetition periods) between the visit start time
	///		   and the start of the first sub-array image. Should be called if StrayLightGenerator,
	///		   JitterProducer or OrbitSimulator are switched on. Otherwise the value is 0.
	void setPreSubArrayTimeSteps() {m_numberOfPreSubArrayTimeSteps = static_cast<unsigned>(ceil((getStartTimeUTC() - getVisitStartTimeUTC()).getSeconds() / m_repetitionPeriod));}

	/// @brief Returns the number of time steps (repetition periods) between the visit start time and the start of the first sub-array image
	unsigned getNumberOfPreSubArrayTimeSteps() const {return m_numberOfPreSubArrayTimeSteps;}

	/// @brief Returns the time step corresponding to the initial full frame image
	int getFullFrameTimeStep() const {return -1*static_cast<int>(round((m_repetitionPeriod+5.)/m_repetitionPeriod));}

	/// @brief returns the number of time steps (unstacked images) for which the payload is in the loop
	unsigned getNumberOfTimeStepsPITL() const {return m_numberOfStackedImages_PITL*m_exposuresPerStack;}

	/// @brief returns the total duration of the simulation in seconds
	double getDuration() const {return m_numberOfStackedImages*m_exposuresPerStack*m_repetitionPeriod;}

	/// @brief returns the number of seconds since the start of the first sub-array image for the specified time step
	double getTimeSinceStart(int iTimeStep) const {return iTimeStep*getRepetitionPeriod();}

	/// @brief returns the number of seconds since the start of the visit for the specified time step
	double getTimeSinceVisitStart(int iTimeStep) const {return getTimeSinceStart(iTimeStep) + (getStartTimeUTC()-getVisitStartTimeUTC()).getSeconds();};

	/// @brief returns the number of stacked image cubes
	unsigned getNumberOfStackedImageCubes() const;

	/// @brief returns the number of unstacked image cubes
	unsigned getNumberOfUnstackedImageCubes() const;

	/// @brief Sets the number of stacked images for which the payload is in the loop
	void setNumberOfStackedImagesPITL(unsigned numberOfStackedImages_PITL) {m_numberOfStackedImages_PITL = numberOfStackedImages_PITL;}

	 /// @brief Returns the number of imagettes per stack for imagette stacking
	unsigned getImagetteStackingNumber() const {return m_imagetteStackingNumber;}

	/// @brief Returns whether or not a full frame image will be generated before the sub-array sequence
	bool doFullFrame() const {return m_doFullFrame;}

	/// @brief Sets the imagette stacking number
	void setImagetteStackingNumber();

private:
	boost::posix_time::ptime m_visitStartTime; ///< Start time of the visit
	boost::posix_time::ptime m_startTime; ///< Start time for the sequence of sun-array images
	unsigned m_numberOfStackedImages; ///< Total number of stacked images
	unsigned m_exposuresPerStack; ///< Number of exposures per stacked image
	double m_exposureTimeAsDouble; ///< Exposure duration
	unsigned m_maxImagesPerCube; ///< Maximum number of images per image cube
	unsigned m_numberOfStackedImages_PITL; ///< Number of stacked images for which the payload is in the loop
	unsigned m_imagetteStackingNumber; ///< Number of imagettes per stack
	double m_repetitionPeriod; ///< Exposure repetition period (interval between consecutive exposure start times)
	bool m_doFullFrame; ///< Boolean to indicate whether or not the sub-array images are preceded by a full frame image
	unsigned m_numberOfPreSubArrayTimeSteps; ///< Number of time steps (repetition periods) between the visit start time and the start of the first sub-array image
};

#endif /* _TIME_CONFIGURATION_HXX_ */
