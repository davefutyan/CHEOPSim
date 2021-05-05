/*
 * ExampleModule.cxx
 *
 *  Created on: Dec 19, 2013
 *      Author: futyand
 */

#include "ExampleModule.hxx"

using namespace std;

//Use the following method to initialize the module by reading in parameters from the configuration file
void ExampleModule::initialize(const ModuleParams& params) {

	//Use GetAsDouble, GetAsInt, GetAsBoolean or GetAsString according to the
	//parameter type specified in the configuration xml file
	m_parameter1 = params.GetAsDouble("parameter1");
	m_parameter2 = params.GetAsInt("parameter2");

	cout << "initializing ExampleModule.  Input parameters are:" << endl;
	cout << "   parameter1: " << m_parameter1 << endl;
	cout << "   parameter2: " << m_parameter2 << endl;
}

//Use the following method to implement processing which is to be performed at each time step of the simulation
void ExampleModule::process(Data* data, int timeStep, bool fullFrame) const {

	//Print target star flux as a function of time
    unsigned currentTime = data->getTimeConfiguration().getTimeSinceStart(timeStep);
    const Star * targetStar = data->getFieldOfView()->getStars()[0];
    double flux = targetStar->getMeanFlux() * targetStar->getTimeSeries()[timeStep]->getTransitFluxFactor();
    cout << "Processing ExampleModule (time loop): Target star flux after " << currentTime << " seconds = " << flux << endl;

}
