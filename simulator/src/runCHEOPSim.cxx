/*
 * testCHEOPSim2.cxx
 *
 *  Created on: Dec 17, 2013
 *      Author: futyand
 */

#include <iostream>

#include "Simulator.hxx"

using namespace std;

int main(int argc, char * argv [])
{
	try {

		//Create an instance of ProgramParams which reads the user configuration from conf/testCHEOPSim.xml
		ParamsPtr params = CheopsInit(argc,argv);

		//Create an instance of the Simulator and run it
		Simulator *simulator = new Simulator(params);
		simulator->process();
		delete simulator;

	} catch(std::exception& e) {
		throw runtime_error(e.what());
	}

    return  0;
}

