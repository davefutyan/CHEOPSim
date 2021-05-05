/*
 * Simulator.hxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief This class is used to configure and run CHEOPSim.
///
/// The CHEOPSim configuration is passed to the Simulator constructor
/// through an instance of the ProgramParams class, in turn constructed using
/// the configuration xml file.  The following parameters must be provided
/// in the configuration file:
///		1. startTime(double): Julian start date of the simulation
///		2. duration (int): Duration of the simulation in seconds.
///		   Must be an integer multiple of exposureTime*exposuresPerStack
///		3. exposuresPerStack (int) Number of exposures to be stacked
///		4. exposureTime (int): Exposure time in seconds
///		5. modulesToRun (string): Comma separated list of the names
///								  of modules to be run
///		6. ModuleParams parameter set for each module.  The \<name\>
///		   parameter must match:
///				(a) the name in moduesToRun
///				(b) the name mapping to module to the corresponding class
///					instantiation in ModuleRegistry::addModule(string name).
///		   Refer to the documentation of the relevant module(s) for module
///		   specific parameters.
///
/// Once the Simulator constructor has been called, CHEOPSim is run by calling the
/// method Simulator::process().  This executes the method Module::process()
/// for each of the sequence of modules defined by the user in the modulesToRun
/// parameter in the configuration xml file.
///
/// Modules with Module::Schedule=begin are run first,
/// followed by a loop over time steps in which the sequence of modules with
/// Module::Schedule=timeLoop is run for each time step.  Modules with
/// Module::Schedule=end are run after the time loop.  Aside from these
/// conditions, modules are run in the order defined in modulesToRun in the
/// configuration xml file.
///
/// Instantiation of the classes deriving from Module is performed in
/// ModuleRegistry::addModule(string name).  New modules classes thus need
/// to be instantiated there - see ModuleRegistry documentation for details.
///
////////////////////////////////////////////////////////////////////////

#ifndef _SIMULATOR_HXX_
#define _SIMULATOR_HXX_

#include <pqxx/pqxx>

#include "Module.hxx"

class Simulator {
public:

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Constructor takes the configuration defined in an instance of
	///        ProgramParams as input
	///
	/// @param [in] params:  Pointer to an instance of ProgramParams containing the
	///						 CHEOPSim input parameters read in from the CHEOPSim
	///						 configuration xml file
	Simulator(const ParamsPtr params);
	virtual ~Simulator();

	/// @brief Runs the Simulator by initializing and running all Modules
	void process();

	/// @brief Returns the list of modules in the ModuleRegistry
	std::vector<Module*> getModules() {return m_modules; }

private:

	/// @brief Initializes the member instances of ProgramParams, TimeConfiguration and Data
	void initialize();
	/// @brief Initializes the connection to the job submission database
	void connectToDatabase(string databaseServer);
	/// @brief Initializes all modules
	void initializeModules();
	/// @brief Generates bias voltages and temperatures with drift and random fluctuations
	///	       at 1.2s cadence, to be used for image metadata and HK data
	void generateVoltagesAndTemperatures();
	/// @brief Runs the method Module::process() for all modules with Module::Schedule=begin
	void processBegin();
	/// @brief Loops over time steps and for each time step runs the method Module::process()
	///        for all modules with Module::Schedule=timeLoop
	void processTimeLoop();
	/// @brief Runs the method Module::process() for all modules with Module::Schedule=end
	void processEnd();

	const ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	vector<Module*> m_modules; ///< List of modules to be run
	Data * m_data; ///< Pointer to the Data
	pqxx::connection * m_database; ///< Connection to the job submission database
	bool m_modulesInitialized; ///< Flag to indicate whether or not modules have been initialized
	bool m_doFullFrame; ///< Flag to indicate of a full frame image should be generated
	string m_jobid; ///< Job ID for the simulation, extracted from the current directory name
};

#endif /* _SIMULATOR_HXX_ */
