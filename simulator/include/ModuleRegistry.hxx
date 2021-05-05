/*
 * ModuleRegistry.hxx
 *
 *  Created on: Jan 14, 2014
 *      Author: futyand
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Class used to instantiate modules
///
/// The list of modules to be run is read from the parameter xml file via the global
/// parameter modulesToRun, which consists of a comma separated list of modules.
/// The list is parsed within the instanciateModules method and the addModule
/// method is called for each module in the list.
///
/// Note that ModuleRegistry::addModule is the only place where classes derived
/// from module are referred to explicitly, rather than being accessed via the Module base class.
/// Hence once a new class deriving from Module is defined, it is sufficient to
/// perform the following three steps in order to include the new module in the simulator:
///		1. Instantiate the class in ModuleRegistry::addModule as follows:
///		   else if (name == "moduleName") module = new className(\[optional constructor parameters]);
///		2. Add the parameters of the module to the parameter xml file, where the \<name\> attribute
///		   must match moduleName in step 1.
///		3. Add moduleName to the list of modules in the global parameter modulesToRun in the xml file.
///
////////////////////////////////////////////////////////////////////////

#ifndef _MODULE_REGISTRY_HXX_
#define _MODULE_REGISTRY_HXX_

#include "Module.hxx"

class ModuleRegistry {
public:
	ModuleRegistry() {};
	virtual ~ModuleRegistry() {};

	/////////////////////////////////////////////////////////////////////////////////
	/// @brief Instantiation of classes deriving from Module is performed exclusively in this method
	///
	/// @param [in] params:  Pointer to an instance of ProgramParams containing the
	///						 CHEOPSim input parameters read in from the CHEOPSim
	///						 configuration xml file
	void instanciateModules(const ParamsPtr params);

	/// @brief Returns the list of instantiated modules for inclusion in Simulator
	std::vector<Module*> getModules() {return m_modules;}

private:
	/// @brief Creates an instance of module corresponding to the specified name and adds it to the list
	void addModule(string name);

	std::vector<Module*> m_modules; ///< List of modules
};

#endif /* _MODULE_REGISTRY_HXX_ */
