/*
 * SatelliteData.hxx
 *
 *  Created on: Dec 20, 2013
 *      Author: futyand
 */

#ifndef _SATELLITE_DATA_HXX_
#define _SATELLITE_DATA_HXX_

#include <vector>
#include <stdexcept>

#include "data/include/SkyPosition.hxx"

using namespace std;

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Container for information relating to the satellite conditions
///		   at a given time (temperature, roll angle, jitter, stray light,
///        particle flux)
///
/// A vector of SatelliteData
/// (each vector element corresponds to a different time step) can be
/// accessed from the Data class via Data::getSatelliteData()
///
////////////////////////////////////////////////////////////////////////

class SatelliteData {
public:

	/// @brief readout hardware: main or redundant
	enum READOUT_HARDWARE{main,redundant};

	/// @brief Bias voltage type
	enum VOLTAGE_TYPE{
		VOD, ///< Output drain level voltage
		VRD, ///< Reset drain level voltage
		VOG, ///< Output gate level voltage
		VSS, ///< Substrate level voltage
		TEMP ///< CCD temperature
	};

	static constexpr double kDefaultCcdTemperature = 233.15; ///< Default CCD temperature
	static constexpr double kDefaultFeeBiasTemperature = 263.15; ///< Default front end electronics bias temperature
	static constexpr double kDefaultFeeAdcTemperature = 265.15; ///< Default front end electronics ADC temperature
	static constexpr double kDefaultTelescopeTemperature = 263.15; ///< Default telescope temperature
	static constexpr double kDefaultDpuTemperature = 294.15; ///< Default DPU temperature
	static constexpr double kDefaultDpuVoltage = -5.121; ///< Default DPU voltage

	static constexpr double kCcdTemperatureSigma[2] = {1.11,1.16}; ///< RMS of CCD temperature fluctuations for main and redundant readout channels (mK)
	static constexpr double kFeeBiasTemperatureSigma[2] = {0.51,0.52}; ///< RMS of FEE Bias temperature fluctuations for main and redundant readout channels (mK)
	static constexpr double kFeeAdcTemperatureSigma[2] = {0.7,1.63}; ///< RMS of FEE ADC temperature fluctuations for main and redundant readout channels (mK)

	//NOTE: for the VOD, VOG and VRD voltages, the following sigmas and drifts apply to voltage-VSS rather than the voltage itself
	static constexpr double kVoltageSigma[2][4] = {{21.82,16.17,11.61,6.61},{20.21,15.75,11.41,6.11}}; ///< Voltage fluctuation RMS, for each voltage type, for main and redundant readout channels (muV)
	static constexpr double kVoltageDrift[2][4] = {{28.42,-2.15,-3.13,11.1},{-2.77,-4.9,-10.13,9.37}}; ///< Voltage drift, for each voltage type, for main and redundant readout channels (muV/day)

	/// @brief Absolute Pointing Error due to jitter
	struct APE {
		/** *************************************************************************
		 *  @brief Constructor for APE
		 *
		 *  @param [in] X  X APE (roll offset in arcseconds)
		 *  @param [in] Y  Y APE (pitch offset in arcseconds)
		 *  @param [in] Z  Z APE (yaw offset in arcseconds)
		 */
		APE(double X, double Y, double Z) : m_X(X), m_Y(Y), m_Z(Z) {};
		APE() : m_X(0.), m_Y(0.), m_Z(0.) {}; //Default constructor
		double getRollOffset() {return m_X;} ///< X APE (roll offset in arcseconds)
		double getPitchOffset() {return m_Y;} ///< Y APE (roll offset in arcseconds)
		double getYawOffset() {return m_Z;} ///< Z APE (roll offset in arcseconds)
		double m_X; ///< X APE (roll offset in arcseconds)
		double m_Y; ///< Y APE (roll offset in arcseconds)
		double m_Z; ///< Z APE (roll offset in arcseconds)
	};

	SatelliteData();
	virtual ~SatelliteData() {};

	/** *************************************************************************
	 *  @brief Add roll angle value to the list of jittered roll angles
	 *  	   corresponding to the current exposure
	 *
	 *  @param [in] rollAngle  jittered roll angle in degrees
	 */
	void addRollAngle(double rollAngle) {m_rollAngles.push_back(rollAngle);}

	/** *************************************************************************
	 *  @brief Return the list of jittered roll angles corresponding to the
	 *  	   current exposure (degrees)
	 */
	vector<double> getRollAngles() const {return m_rollAngles;}

	/** *************************************************************************
	 *  @brief Append the list of jittered pointing directions corresponding to
	 *  	   the current exposure
	 *
	 *  @param [in] pointingDirection  jittered pointing direction
	 */
	void addPointingDirection(SkyPosition pointingDirection) {m_pointingDirections.push_back(pointingDirection);}

	/** *************************************************************************
	 *  @brief Return the list of jittered pointing directions corresponding to
	 *  	   the current exposure
	 */
	vector<SkyPosition> getPointingDirections() const {return m_pointingDirections;}

	/** *************************************************************************
	 *  @brief Add values of jitter offsets together with corresponding
	 *  	   AOCS and science validity flags the lists of these quantities
	 *  	   corresponding to the current exposure
	 *
	 *  @param [in] ape  APE offsets (roll, pitch and yaw) due to jitter
	 *  @param [in] validAocs  AOCS validity flag from the input jitter file
	 *  @param [in] validScience  Science validity flag from the input jitter file
	 */
	void addJitterAPE(APE ape, bool validAocs, bool validScience);

	/** *************************************************************************
	 *  @brief Return the list of APE offsets due to jitter corresponding to the
	 *  	   current exposure
	 */
	vector<APE> getJitterAPEs() const {return m_jitterAPEs;}

	/** *************************************************************************
	 *  @brief Return the list of AOCS validity flags corresponding to the
	 *  	   current exposure
	 */
	vector<bool> getJitterValidAocsFlags() const {return m_validAocsFlags;}

	/** *************************************************************************
	 *  @brief Return the list of science validity flags corresponding to the
	 *  	   current exposure
	 */
	vector<bool> getJitterValidScienceFlags() const {return m_validScienceFlags;}

	/** *************************************************************************
	 *  @brief Returns true if the AOCS flag is valid for every second of the
	 *  	   current exposure
	 */
	bool validAocs() const;

	/** *************************************************************************
	 *  @brief Returns true if the science flag is valid for every second of the
	 *  	   current exposure
	 */
	bool validScience() const;

	/** *************************************************************************
	 *  @brief Returns true if both the AOCS and science flags are valid for
	 *  	   every second of the current exposure
	 */
	bool validAocsScience() const {return (validAocs() && validScience());};

	/** *************************************************************************
	 *  @brief Returns the value of the CCD temperature for the current time step
	 */
	double getCcdTemperature() const {return m_ccdTemperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the CCD temperature for the current time step
	 */
	void setCcdTemperature(double temperature) {m_ccdTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the front end electronics bias temperature for the
	 *  	   current time step
	 */
	void setFeeBiasTemperature(double temperature) {m_feeBiasTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Returns the value of the CCD temperature for the current time step
	 *  	   after applying Gaussian fluctuation
	 */
	double getFluctuatedCcdTemperature() {return m_fluctuatedCcdTemperature;}

	/** *************************************************************************
	 *  @brief Returns the value of the front end electronics bias temperature for the
	 *  	   current time step after applying Gaussian fluctuation
	 */
	double getFluctuatedFeeBiasTemperature() {return m_fluctuatedFeeBiasTemperature;}

	/** *************************************************************************
	 *  @brief Returns the value of the front end electronics ADC temperature for the
	 *  	   current time step after applying Gaussian fluctuation
	 */
	double getFluctuatedFeeAdcTemperature() {return m_fluctuatedFeeAdcTemperature;}

	/** *************************************************************************
	 *  @brief Returns the value of the specified bias voltage for the
	 *  	   current time step after applying Gaussian fluctuation
	 *
	 *  @param [in] type  VOD, VRD, VOG or VSS
	 */
	double getFluctuatedVoltage(VOLTAGE_TYPE type) {return m_fluctuatedVoltage[type];}

	/** *************************************************************************
	 *  @brief Sets the value of the Gaussian fluctuated CCD temperature
	 *  	   for the current time step
	 */
	void setFluctuatedCcdTemperature(double temperature) {m_fluctuatedCcdTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the Gaussian fluctuated FEE bias temperature
	 *  	   for the current time step
	 */
	void setFluctuatedFeeBiasTemperature(double temperature) {m_fluctuatedFeeBiasTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the Gaussian fluctuated FEE ADC temperature
	 *  	   for the current time step
	 */
	void setFluctuatedFeeAdcTemperature(double temperature) {m_fluctuatedFeeAdcTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the specified Gaussian fluctuated bias voltage
	 *  	   for the current time step
	 *
	 *  @param [in] type  VOD, VRD, VOG or VSS
	 *  @param [in] voltage  voltage value to be set
	 */
	void setFluctuatedVoltage(VOLTAGE_TYPE type, double voltage) {m_fluctuatedVoltage[type] = voltage;}

	/** *************************************************************************
	 *  @brief Returns the value of the telescope temperature for the current
	 *  	   time step
	 */
	double getTelescopeTemperature() const {return m_telescopeTemperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the telescope temperature for the current
	 *         time step
	 */
	void setTelescopeTemperature(double temperature) {m_telescopeTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Returns the value of the DPU temperature for the current time step
	 */
	double getDpuTemperature() const {return m_dpuTemperature;}

	/** *************************************************************************
	 *  @brief Sets the value of the DPU temperature for the current time step
	 */
	void setDpuTemperature(double temperature) {m_dpuTemperature = temperature;}

	/** *************************************************************************
	 *  @brief Returns the value of the angle between the pointing direction and
	 *  	   the moon for the current time step
	 */
	double getMoonAngle() const {return m_moonAngle;}

	/** *************************************************************************
	 *  @brief Sets the value of the angle between the pointing direction and
	 *  	   the moon for the current time step
	 */
	void setMoonAngle(double angle) {m_moonAngle = angle;}

	/** *************************************************************************
	 *  @brief Returns the value of the angle between the pointing direction and
	 *  	   the sun for the current time step
	 */
	double getSunAngle() const {return m_sunAngle;}

	/** *************************************************************************
	 *  @brief Sets the value of the angle between the pointing direction and
	 *  	   the sun for the current time step
	 */
	void setSunAngle(double angle) {m_sunAngle = angle;}

	/** *************************************************************************
	 *  @brief Returns the value of the angle between the pointing direction and
	 *  	   the Earth limb for the current time step
	 */
	double getEarthLimbAngle() const {return m_earthLimbAngle;}

	/** *************************************************************************
	 *  @brief Sets the value of the angle between the pointing direction and
	 *  	   the Earth limb for the current time step
	 */
	void setEarthLimbAngle(double angle) {m_earthLimbAngle = angle;}

	/** *************************************************************************
	 *  @brief Returns the value of the spacecraft latitude for the current time step
	 */
	double getLatitude() const {return m_latitude;}

	/** *************************************************************************
	 *  @brief Sets the value of the spacecraft latitude for the current time step
	 */
	void setLatitude(double latitude) {m_latitude = latitude;}

	/** *************************************************************************
	 *  @brief Returns the value of the spacecraft longitude for the current time step
	 */
	double getLongitude() const {return m_longitude;}

	/** *************************************************************************
	 *  @brief Sets the value of the spacecraft longitude for the current time step
	 */
	void setLongitude(double longitude) {m_longitude = longitude;}

	/** *************************************************************************
	 *  @brief Returns whether or not the satellite is in Earth occultation
	 *  	   for the current time step
	 */
	bool getEarthOccultationFlag() const {return m_earthOccultationFlag;}

	/** *************************************************************************
	 *  @brief Sets the value of the flag to indicate whether or not the
	 *  	   satellite is in Earth occultation for the current time step
	 */
	void setEarthOccultationFlag(bool flag) {m_earthOccultationFlag = flag;}

	/** *************************************************************************
	 *  @brief Returns the stray light flux in photons/pix/s for the
	 *  	   current time step
	 */
	double getStrayLightFlux() const {return m_strayLightFlux;}

	/** *************************************************************************
	 *  @brief Sets the stray light flux in photons/pix/s for the current time step
	 */
	void setStrayLightFlux(double flux) {m_strayLightFlux = flux;}

	/** *************************************************************************
	 *  @brief Returns whether or not the stray light is above threshold for the
	 *  	   current time step
	 */
	bool getStrayLightFlag() const {return m_strayLightFlag;}

	/** *************************************************************************
	 *  @brief Sets the value of the flag to indicate whether or not the
	 *  	   stray light is above threshold for the current time step
	 */
	void setStrayLightFlag(bool flag) {m_strayLightFlag = flag;}

	/** *************************************************************************
	 *  @brief Returns whether or not the satellite is over the SAA for the
	 *  	   current time step
	 */
	bool getSAAFlag() const {return m_SAAFlag;}

	/** *************************************************************************
	 *  @brief Sets the value of the flag to indicate whether or not the
	 *  	   satellite is within the SAA for the current time step
	 */
	void setSAAFlag(bool flag) {m_SAAFlag = flag;}

	/** *************************************************************************
	 *  @brief Returns whether or not the
	 *  	   stacked image corresponding to the timestep is to be discarded
	 *  	   due to SAA, Earth occultation, stray light, option to omit first N images,
	 *  	   or option to only write out every Nth image
	 */
	bool getDiscardFlag() const {return m_discardFlag;}

	/** *************************************************************************
	 *  @brief Sets the value of the flag to indicate whether or not the
	 *  	   stacked image corresponding to the timestep is to be discarded
	 *  	   due to SAA, Earth occultation, stray light, option to omit first N images,
	 *  	   or option to only write out every Nth image
	 */
	void setDiscardFlag(bool flag) {m_discardFlag = flag;}

private:
	double m_ccdTemperature; ///< CCD temperature for current time step
	double m_feeBiasTemperature; ///< Front end electronics bias temperature for current time step
	double m_telescopeTemperature; ///< Telescope temperature for current time step
	double m_fluctuatedCcdTemperature; ///< CCD temperature for current time step after applying Gaussian fluctuation
	double m_fluctuatedFeeBiasTemperature; ///< Front end electronics bias temperature for current time step after applying Gaussian fluctuation
	double m_fluctuatedFeeAdcTemperature; ///< Front end electronics ADC temperature for current time step after applying Gaussian fluctuation
	double m_fluctuatedVoltage[4]; ///< Bias voltage value for each voltage type, for current time step after applying Gaussian fluctuation
	double m_dpuTemperature; ///< DPU temperature used to assign ADC_TEMP1
	double m_moonAngle; ///< Angle between pointing direction and moon for current time step
	double m_sunAngle; ///< Angle between pointing direction and sun for current time step
	double m_earthLimbAngle; ///< Angle between the pointing direction and the Earth limb for the current time step
	double m_latitude; ///< Latitude of the spacecraft
	double m_longitude; ///< Latitude of the spacecraft
	bool m_earthOccultationFlag; ///< Flag to indicate whether or not the satellite is in Earth occultation for the current time step
	bool m_SAAFlag; ///< Flag to indicate whether or not the satellite is within the SAA for the current time step
	double m_strayLightFlux; ///< Stray light flux in photons/pixel/second for the current time step
	bool m_strayLightFlag; ///< Flag to indicate whether or not the stray light is above threshold for the current time step
	bool m_discardFlag; ///< Flag to indicate whether or not the stacked image corresponding to the timestep is to be discarded due to SAA, Earth occultation, stray light, option to omit first N images, or option to only write out every Nth image
	vector<APE> m_jitterAPEs; ///< List of APE offsets due to jitter corresponding to the current exposure
	vector<bool> m_validAocsFlags; ///< List of AOCS validity flags corresponding to the current exposure
	vector<bool> m_validScienceFlags; ///< List of science validity flags corresponding to the current exposure
	vector<double> m_rollAngles; ///< List of jittered roll angles corresponding to the current exposure
	vector<SkyPosition> m_pointingDirections; ///< List of jittered pointing directions corresponding to the current exposure
};

#endif /* _SATELLITE_DATA_HXX_ */
