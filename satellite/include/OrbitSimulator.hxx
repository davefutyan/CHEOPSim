/*
 * OrbitSimulator.hxx
 *
 *  Created on: 14 Feb 2014
 *      Author: futyand
 */

/** *************************************************************************
 *  @author David Futyan UGE
 *
 *  @brief This module simulates the variation of relevant quantities which
 *  	   depend on the position of the satellite within the its orbit.
 *
 *  The module models the rotation of the field of view, the variation of
 *  telescope temperature (affecting PSF breathing) and CCD temperature
 *  (affecting dark current and gain). The module also outputs the position
 *  and velocity of the satellite in its orbit.
 */

#ifndef _ORBIT_SIMULATOR_HXX_
#define _ORBIT_SIMULATOR_HXX_

#include <Eigen/Dense>

#include "SCI_RAW_Attitude.hxx"
#include "MPS_PRE_VisitConstraints.hxx"
#include "EXT_APP_SAAMap.hxx"

#include "OrbitData.hxx"

#include "simulator/include/Module.hxx"

class OrbitSimulator: public Module {
public:

	/// @brief CHEOPS component
	enum COMPONENT{
		telescope, ///< telescope
		ccd, ///< CCD
		fee, ///< Front end electronics
		dpu	 ///< DPU
	};

	static const unsigned kCheopsOrbitPeriod = 5910;  ///< CHEOPS orbit period in seconds
	static const unsigned kCheopsOrbitTimeStep = 60;  ///< Time step in seconds for calculating CHEOPS orbit
	static constexpr double kDpuTemperatureAmplitude = 1.; ///< Amplitude of sinusoidal variation for DPU temperature

	OrbitSimulator() : Module("OrbitSimulator",timeLoop), m_SAAMap(nullptr), m_attitude(nullptr), m_MPSVisitConstraints(nullptr) {};
	virtual ~OrbitSimulator();

	void initialize(const ModuleParams & params);
	void doBegin(Data * data, bool fullFrame=false);
	void process(Data * data, int timeStep, bool fullFrame=false) const {};

	/** *************************************************************************
	 *  @brief Reads the input temperature time series
	 *
	 *  @param [in] component  CHEOPS component: telescope or CCD
	 */
	void readTemperatures(COMPONENT component);

	/** *************************************************************************
	 *  @brief Returns the current temperature based on the file input time series
	 *
	 *  @param [in] currentTime  Number of seconds since the start of the simulation
	 *  @param [in] component  CHEOPS component: telescope or CCD
	 */
	double getTemperatureFromFile(double currentTime, COMPONENT component) const;

	/** *************************************************************************
	 *  @brief Returns the current temperature according to sinusoidal variation
	 *
	 *  @param [in] currentTime  Number of seconds since the start of the simulation
	 *  @param [in] component  CHEOPS component: telescope or CCD
	 */
	double getSinusoidTemperature(double currentTime, COMPONENT component) const;

	/** *************************************************************************
	 *  @brief Reads the input orbit data time series. The data is stored in
	 *  	   memory, and is also written out at 5 minute intervals to a
	 *  	   AUX_RES_Orbit output data structure.
	 *
	 *  @param [in] data  Pointer to the data
	 */
	void readOrbit(Data * data);

	/** *************************************************************************
	 *  @brief Reads the SAA map and stores the information in memory.
	 */
	void readSAAMap();

	/** *************************************************************************
	 *  @brief Sets the spacecraft orbit position for each minute of the
	 *  	   simulation by linear interpolation from the time resolution of
	 *  	   the input file
	 *
	 *  @param [in] orbitData_ref  Vector of input orbit data
	 */
	void setOrbitPositions(const vector<OrbitData> orbitData_ref);

	/** *************************************************************************
	 *  @brief Sets the spacecraft orbit velocity for each time step based on the
	 *  	   delta between the positions at the adjacent time steps
	 */
	void setOrbitVelocities();

	/** *************************************************************************
	 *  @brief Sets the angle between the pointing direction and the moon
	 *  	   for each time step
	 *
	 *  @param [in] startTime  UTC start time of the simulation
	 *  @param [in] endTime  UTC end time of the simulation
	 */
	void setMoonAngles(UTC startTime, UTC endTime);

	/** *************************************************************************
	 *  @brief Sets the angle between the pointing direction and the moon
	 *  	   for each time step
	 *
	 *  @param [in] startTime  UTC start time of the simulation
	 *  @param [in] endTime  UTC end time of the simulation
	 */
	void setSunAngles(UTC startTime, UTC endTime);

	/** *************************************************************************
	 *  @brief Calculates the jittered roll angle and pointing direction for each time step
	 *
	 *  @param [in] data  Pointer to the Data
	 *	@param [in] timeStep  Index of the current time step
	 */
	void setRollAnglesAndPointingDirections(Data* data, int timeStep) const;

	/** *************************************************************************
	 *  @brief Returns the rotation matrix to transform from the inertial frame
	 *  	   to the satellite frame at the specified time
	 *
	 *  @param [in] currentTime  UTC time for which the matrix is to be calculated
	 */
	Eigen::Matrix3d inertialToSatMatrix(UTC currentTime) const;

	/** *************************************************************************
	 *  @brief Calculates the jittered roll angle
	 *
	 *  @param [in] M_inertialToSat  Rotation matrix to transform from the
	 *  							 inertial frame to the satellite frame
	 *  @param [in] rollOffset  Roll angle offset due to jitter (X APE) (arcseconds)
	 */
	double calculateRollAngle(Eigen::Matrix3d M_inertialToSat, double rollOffset) const;

	/** *************************************************************************
	 *  @brief Calculates the jittered pointing direction
	 *
	 *  @param [in] M_inertialToSat  Rotation matrix to transform from the
	 *  							 inertial frame to the satellite frame
	 *  @param [in] ape  Jitter Absolute Pointing Errors (arcseconds)
	 */
	SkyPosition calculatePointingDirection(Eigen::Matrix3d M_inertialToSat, SatelliteData::APE ape) const;

	/** *************************************************************************
	 *  @brief Prints the angle between the pointing direction and the orbital
	 *         plane to the log file
	 *
	 *  @param [in] pointingDirection  Pointing direction
	 */
	void printAngleToOrbitalPlane(SkyPosition pointingDirection);

	/** *************************************************************************
	 *  @brief Returns the larger of two double precision values
	 *
	 *  @param [in] a  first value
	 *  @param [in] b  second value
	 */
	double max(double a, double b) const {return (a > b) ? a : b;}

private:
	bool m_ccdTemperatureVariation; ///< Flag to indicate whether or not to apply CCD temperature variation
	double m_ccdMeanTemperature; ///< CCD mean temperature in Kelvin
	double m_ccdTemperatureAmplitude; ///< Amplitude of sinusoidal variation of CCD temperature in Kelvin
	double m_ccdTemperaturePeriod; ///< Period of sinusoidal variation of CCD temperature in seconds
	bool m_feeTemperatureVariation; ///< Flag to indicate whether or not to apply FEE temperature variation
	double m_feeMeanTemperature; ///< FEE mean temperature in Kelvin
	double m_feeTemperatureAmplitude; ///< Amplitude of sinusoidal variation of FEE temperature in Kelvin
	double m_feeTemperaturePeriod; ///< Period of sinusoidal variation of FEE temperature in seconds
	bool m_telescopeTemperatureVariation; ///< Flag to indicate whether or not to apply telescope temperature variation
	double m_telescopeMeanTemperature; ///< telescope mean temperature in Kelvin
	double m_telescopeTemperatureAmplitude; ///< Amplitude of sinusoidal variation of telescope temperature in Kelvin
	double m_telescopeTemperaturePeriod; ///< Period of sinusoidal variation of telescope temperature in seconds
	bool m_telescopeTemperatureFromFile; ///< Flag to indicate whether or not to read the telescope temperature time series from a file
	string m_telescopeTemperatureFilename; ///< Filename for telescope temperature vs time
	bool m_rotateFOV; ///< Flag to indicate whether or not to apply FOV rotation
	string m_orbitFilename; ///< Filename for orbit position time series
	double m_minAngleToOrbitalPlane; ///< Minimum angle (degrees) between pointing direction and orbital plane for roll angle calculation
	int m_attitudeCadence; ///< Cadence in seconds for writing out spacecraft attitude data to SCI_RAW_Attitude data structure
	int m_orbitCadence; ///< Cadence in minutes for writing out spacecraft orbit data to AUX_RES_Orbit data structure
	string m_SAAMapFilename; ///< Filename for FITS file containing the SAA map
	bool m_saaFlagFromVisitConstraints; ///< Flag to indicate whether to read the SAA from MPS_PRE_VisitConstraints rather than calculating it from the latitude/longitude in the orbit file

	vector<double> m_telescopeTemperature; ///< Telescope temperature for each entry in the time series read from the input file
	vector<double> m_telescopeTime; ///< Time in seconds corresponding to each entry in the telescope temperature time series read from the input file
	vector<OrbitData> m_orbitData; ///< Vector of orbit data with one element per minute
	Eigen::Vector3d m_lineOfSight; ///< Line of sight vector for telescope pointing direction in inertial frame coordinates

	int m_SAAMap_lat[90]; ///< Latitude values for SAA map
	int m_SAAMap_long[121]; ///< Longitude values for SAA map
	bool m_SAAMap_value[90][121]; ///< SAA flag values for each lat/long bin

	ExtAppSaamap * m_SAAMap; ///< Input FITS table containing the SAA map
	SciRawAttitude * m_attitude; ///< Output FITS table containing spacecraft attitude
	MpsPreVisitconstraints * m_MPSVisitConstraints; ///< Output FITS table containing mission planning visit constraints

};

#endif /* _ORBIT_SIMULATOR_HXX_ */
