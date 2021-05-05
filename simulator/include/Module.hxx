/*
 * Module.hxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Base class for all CHEOPSim modules.
///
/// Classes derived from Module should be instantiated in
/// ModuleRegistry::addModule(string name) where the argument is the name
/// of the module as it appears in the configuration xml file.
///
////////////////////////////////////////////////////////////////////////

#ifndef _MODULE_HXX_
#define _MODULE_HXX_

#include "ProgramParams.hxx"
#include "data/include/Data.hxx"

class Module {
public:

	/// @brief enum to define the scheduling of the module in Simulator::process()
	enum SCHEDULE{begin,    ///< module should run before the time loop
		          timeLoop, ///< module should run for each time step (default)
		          end       ///< module should run after the time loop
	};

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Constructor initializes module name and scheduling
	///
	/// @param [in] name:  String specifying the name of the module as it appears in the
	///                    configuration xml file
	/// @param [in] schedule:  enum (Module::SCHEDULE) with one of the following values:
	///                 - begin if the module should run before the time loop
	///                 - timeLoop if the module should run for each time step (default)
	///                 - end if the module should run after the time loop
	Module(std::string name,SCHEDULE schedule=timeLoop) : m_name(name),m_schedule(schedule) {};

	virtual ~Module() {};

	////////////////////////////////////////////////////////////////////////////////
	/// @brief Initialize the module by reading in parameters from the configuration file
    ///
	/// @param [in] params:  Pointer to an instance of ProgramParams containing the
	///                      CHEOPSim input parameters read in from the CHEOPSim
	///                      configuration xml file
	virtual void initialize(const ModuleParams & params) = 0;

	////////////////////////////////////////////////////////////////////////////////
	/// @brief Further module initialization with access to data initialized by other modules.
	///        Called after preceding modules have been initialized
	///
	/// @param [in] data:  Pointer to the data
	/// @param [in] fullFrame:  Boolean to indicate whether or not a full frame image is to be generated
	virtual void doBegin(Data * data, bool fullFrame=false) {};

	////////////////////////////////////////////////////////////////////////////////
	/// @brief Processing of the data.  Called for each time step of the simulation if
	///        Module::Schedule is timeLoop (default), otherwise called once.
	///
	/// @param [in] data:  Pointer to the data
	/// @param [in] timeStep:  Index of the current time step (only for Module::Schedule==timeloop)
	/// @param [in] fullFrame:  Boolean to indicate whether or not a full frame image is to be generated
	virtual void process(Data * data, int timeStep, bool fullFrame=false) const = 0;

	////////////////////////////////////////////////////////////////////////////////
	/// @brief For processing of data after the time loop (data reduction)
	///
	/// @param [in] data:  Pointer to the data
	virtual void doEnd(Data * data) const {};

	/// @brief Returns the name of the module
	std::string getName() {return m_name;}

	/// @brief Returns the scheduling of the module in the processing sequence
	SCHEDULE schedule() {return m_schedule;}

protected:
	std::string m_name; ///< Name of the module
	SCHEDULE m_schedule; ///< Scheduling of the module in the processing sequence
};

#endif /* _MODULE_HXX_ */
