/*
 * ModuleRegistry.cxx
 *
 *  Created on: Jan 14, 2014
 *      Author: futyand
 */

#include "boost/tokenizer.hpp"
#include "boost/token_iterator.hpp"
#include "boost/token_functions.hpp"

#include "source/include/StarProducer.hxx"
#include "source/include/TransitFluxModulator.hxx"
#include "source/include/StellarNoiseFluxModulator.hxx"
#include "source/include/StellarVariationFluxModulator.hxx"
#include "source/include/UserFluxModifier.hxx"
#include "source/include/ZodiacalLightGenerator.hxx"
#include "source/include/StrayLightGenerator.hxx"
#include "satellite/include/JitterProducer.hxx"
#include "satellite/include/OrbitSimulator.hxx"
#include "telescope/include/GlobalThroughputGenerator.hxx"
#include "telescope/include/FocalPlaneGenerator.hxx"
#include "telescope/include/PSFGenerator.hxx"
#include "telescope/include/HaloGenerator.hxx"
#include "detector/include/FrameTransferSmearer.hxx"
#include "detector/include/PhotonNoiseGenerator.hxx"
#include "detector/include/FlatFieldGenerator.hxx"
#include "detector/include/DarkCurrentGenerator.hxx"
#include "detector/include/CosmicRayGenerator.hxx"
#include "detector/include/FullWellSimulator.hxx"
#include "detector/include/ChargeTransferSimulator.hxx"
#include "detector/include/BiasGenerator.hxx"
#include "IdealLightCurveGenerator.hxx"
#include "ImageWriter.hxx"
#include "HKWriter.hxx"
#include "DataReduction.hxx"
#include "ExampleModule.hxx"
#include "ModuleRegistry.hxx"

using namespace std;

////////////////////////////////////////////////////////////////////////
/// Classes deriving from Module must be instantiated in this method as follows:
///		   else if (name == "moduleName") module = new className([optional constructor parameters]);
/// where moduleName matches the name attribute of the module in the parameter xml file.
void ModuleRegistry::addModule(string name) {

	Module * module;

	unsigned starIndex = 0;
	if (name.find("_star") != string::npos) {
		starIndex = boost::lexical_cast<unsigned>(name.substr(name.find("_star")+5));
		name = name.substr(0,name.find("_star"));
	}

	if (name == "StarProducer") module = new StarProducer();
	else if (name == "TransitFluxModulator") module = new TransitFluxModulator(starIndex);
	else if (name == "StellarNoiseFluxModulator") module = new StellarNoiseFluxModulator(starIndex);
	else if (name == "StellarVariationFluxModulator") module = new StellarVariationFluxModulator(starIndex);
	else if (name == "UserFluxModifier") module = new UserFluxModifier(starIndex);
	else if (name == "ZodiacalLightGenerator") module = new ZodiacalLightGenerator();
	else if (name == "StrayLightGenerator") module = new StrayLightGenerator();
	else if (name == "JitterProducer") module = new JitterProducer();
	else if (name == "OrbitSimulator") module = new OrbitSimulator();
	else if (name == "GlobalThroughputGenerator") module = new GlobalThroughputGenerator();
	else if (name == "FocalPlaneGenerator") module = new FocalPlaneGenerator();
	else if (name == "PSFGenerator") module = new PSFGenerator();
	else if (name == "HaloGenerator") module = new HaloGenerator();
	else if (name == "FrameTransferSmearer") module = new FrameTransferSmearer();
	else if (name == "PhotonNoiseGenerator") module = new PhotonNoiseGenerator();
	else if (name == "FlatFieldGenerator") module = new FlatFieldGenerator();
	else if (name == "DarkCurrentGenerator") module = new DarkCurrentGenerator();
	else if (name == "CosmicRayGenerator") module = new CosmicRayGenerator();
	else if (name == "FullWellSimulator") module = new FullWellSimulator();
	else if (name == "ChargeTransferSimulator") module = new ChargeTransferSimulator();
	else if (name == "BiasGenerator") module = new BiasGenerator();
	else if (name == "IdealLightCurveGenerator") module = new IdealLightCurveGenerator();
	else if (name == "ImageWriter") module = new ImageWriter();
	else if (name == "HKWriter") module = new HKWriter();
	else if (name == "DataReduction") module = new DataReduction();
	else if (name == "ExampleModule") module = new ExampleModule();
	else throw runtime_error("Module name not recognized: "+name);

	m_modules.push_back(module);

}

void ModuleRegistry::instanciateModules(const ParamsPtr params) {

	// Get the list of modules to be run from the parameter xml file
	string modulesToRun = params->GetAsString("modulesToRun");
	cout << "The following modules will be run: " << modulesToRun << endl;

	// Parse the list of modules, and call addModule for each module in the list
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep(",");
	tokenizer tokens(modulesToRun, sep);
	for (tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter) {
		addModule(string(*tok_iter));
	}

}
