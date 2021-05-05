/*
 * runReferenceFileConversion.cxx
 *
 *  Created on: Jun 10, 2020
 *      Author: David Futyan UGE
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Program to convert reference file data provided in ascii
///		   format to the CHEOPS SOC data structure formats defined in
///		   common_sw/fits_data_model
///
/// Generates executable runReferenceFileConversion
///	Configured using reference_file_generation/conf/runReferenceFileConversion.xml
///
/// Note that The ProcNum and ArchRev keywords are not set by this program
/// and must be assigned by manually editing the output fits file and
/// updating the version number in the output filename accordingly.
///
////////////////////////////////////////////////////////////////////////

#include <iostream>

#include "ProgramParams.hxx"

#include "ReferenceFileConverter.hxx"

using namespace std;

int main(int argc, char * argv [])
{
	try {

		//Create an instance of ProgramParams which reads the user configuration from conf/runInputConverter.xml
		ParamsPtr params = CheopsInit(argc,argv);

		//Create an instance of the converter class and run it
		ReferenceFileConverter *converter = new ReferenceFileConverter(params);
		converter->process();
		delete converter;

	} catch(std::exception& e) {
		throw runtime_error(e.what());
	}

    return  0;
}
