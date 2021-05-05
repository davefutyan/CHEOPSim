/*
 * ExampleModule.hxx
 *
 *  Created on: Dec 19, 2013
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Template for classes derived from Module
///
/// This class can be used as a template for creating a new CHEOPSim module.
/// All CHEOPSim modules inherit from the base class Module.  The class
/// should be instantiated in ModuleRegistry::addModule(string name)
/// where the argument is the name of the module as it appears in the
/// configuration xml file.
///
/// Scheduling, initialization and processing of the module are
/// handled by the Simulator through Simulator::process(),
/// so there is no need for the user to directly call
/// ExampleModule::initialize() and ExampleModule::process(), although
/// implementation of these methods must be provided.
///
////////////////////////////////////////////////////////////////////////

#ifndef _EXAMPLE_MODULE_HXX_
#define _EXAMPLE_MODULE_HXX_

#include "Module.hxx"

class ExampleModule: public Module {
public:

	/////////////////////////////////////////////////////////////////////////////////
	/// The constructor of the Module base class takes two arguments:
	///  1: string specifying the name of the module as it appears in the configuration xml file
	///  2: enum (Module::SCHEDULE) with value - begin if the module should run before the time loop
	///                                        - timeLoop if the module should run for each time step (default)
	///                                        - end if the module should run after the time loop
	ExampleModule() : Module("ExampleModule",end) {};
	virtual ~ExampleModule() {};

	/// @brief Use this method to initialize the module by reading in parameters from the configuration file
	void initialize(const ModuleParams & params);

	/// @brief Use this method to implement processing of the data.  The second argument (timeStep)
	///        is only used if Module::Schedule==timeLoop (default).  In this case, the method is called
	///        for each time step of the simulation.
	void process(Data * data, int timeStep, bool fullFrame=false) const;

private:
	//Add data members here as required.  Typically there will be a data member
	//corresponding to each parameter for the module in the configuration xml file.
	//Naming convention:  m_parameterName
	double m_parameter1; ///<< First parameter
	int m_parameter2; ///<< Second parameter
};

#endif /* _EXAMPLE_MODULE_HXX_ */
