/*
 * runFluxConversion.cxx
 *
 *  Created on: Jun 1, 2020
 *      Author: David Futyan UGE
 */

////////////////////////////////////////////////////////////////////////
/// @author David Futyan UGE
///
/// @brief Program to generate a REF_APP_FluxConversion reference file
///
///	Implements the formulae described here:
/// https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion
///	using parameters calculated using
///	https://svn.isdc.unige.ch/svn-cheops/06_cheopsim/software/CHEOPSim/trunk/resources/FluxConversion.html
///
/// Generates executable runFluxConversion
///	Configured using reference_file_generation/conf/runFluxConversion.xml
///
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <vector>

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/gregorian_types.hpp"
#include "boost/algorithm/string.hpp"

#include "CreateFitsFile.hxx"
#include "REF_APP_FluxConversion.hxx"
#include "REF_APP_SEDTeff.hxx"

#include "source/include/StarProducer.hxx"

using namespace std;

int main(int argc, char * argv [])
{
	try {

		//Create an instance of ProgramParams which reads the user configuration from conf/runFluxConversion.xml
		ParamsPtr params = CheopsInit(argc,argv);

		//Read in the configuration parameters
		double blackBodyOffset = params->GetAsDouble("blackBodyOffset");
		double fitIntercept = params->GetAsDouble("fitIntercept");
		double fitSlope = params->GetAsDouble("fitSlope");
        string throughputFilename = params->GetAsString("throughputFilename");
        string qeFilename = params->GetAsString("QEFilename");

		//Initialize the wavelength dependence class
		WavelengthDependence * wavelengthDependence = new WavelengthDependence(true,throughputFilename,true,qeFilename,
																			   params->GetAsString("SEDFilename"),"");

		//Calculate and fill the contents of REF_APP_FluxConversion given the specified throughput and QE files

		//Open the output FITS file
		UTC startValidity = UTC(to_iso_extended_string(boost::posix_time::time_from_string(params->GetAsString("validityStart"))));
		UTC endValidity = UTC(to_iso_extended_string(boost::posix_time::ptime(boost::gregorian::date(2040,1,1))));
		RefAppFluxconversion * fluxConversion = new RefAppFluxconversion(buildFileName(".", startValidity,
				VisitId(), PassId(), std::string(), RefAppFluxconversion::getExtName()),"CREATE");

		//Calculate value of F0_X, the flux of Vega in the CHEOPS passband
		double integratedThroughputQE_CHEOPS_Vega = wavelengthDependence->getIntegratedThroughputQE(0.,SatelliteData::kDefaultCcdTemperature);
		double F0_X = StarProducer::kNPhoVega * integratedThroughputQE_CHEOPS_Vega;

		//Multiply F0_X by scale factor to account for empirical flux deficit at the effective temperature of Vega
		//(see https://redmine.isdc.unige.ch/projects/cheops/wiki/Flux_Conversion for details)
		double fitScaleFactor = fitIntercept + fitSlope * 9600;
		F0_X *= fitScaleFactor;

		//Set header keywords
		fluxConversion->setKeyF0X(F0_X);
		fluxConversion->setKeyX0(0.);
		fluxConversion->setKeyVStrtU(startValidity);
		fluxConversion->setKeyVStopU(endValidity);

		//Assign ArchRev and ProcNum values for the QE and throughput input reference files as header keywords
        vector<string> substrings;
        boost::split(substrings,throughputFilename,boost::is_any_of("_"));
        fluxConversion->setKeyThrArev(int((substrings[substrings.size()-1])[2] - '0'));
		fluxConversion->setKeyThrPnum(int((substrings[substrings.size()-1])[4] - '0'));
        boost::split(substrings,qeFilename,boost::is_any_of("_"));
		fluxConversion->setKeyQeArev(int((substrings[substrings.size()-1])[2] - '0'));
		fluxConversion->setKeyQePnum(int((substrings[substrings.size()-1])[4] - '0'));

		//Initialize the FluxConverter class, passing it the fit parameters to correct for the
		//empirically observed/predicted flux ratio as a function of effective temperature
		//The default parameter values in runFluxConversion.xml were calculated using
		//https://svn.isdc.unige.ch/svn-cheops/06_cheopsim/software/CHEOPSim/trunk/resources/FluxConversion.html
		FluxConverter * fluxConv = new FluxConverter(blackBodyOffset,fitIntercept,fitSlope);

		//Loop over effective temperatures
		for (int effectiveTemperature=2300; effectiveTemperature <= 40000; effectiveTemperature+=100) {

			//Calculate the difference between CHEOPS magnitude and Gaia magnitude
			WavelengthDependence::REFERENCE_BAND refBand = params->GetAsBool("gaiaBand") ? WavelengthDependence::GaiaBand : WavelengthDependence::Vband;
			double deltaMag = fluxConv->cheopsMinusRefBandMagnitude(refBand,effectiveTemperature,wavelengthDependence);

			//Calculate the electrons to photon ratio
			double integratedThroughput_CHEOPS = wavelengthDependence->getIntegratedThroughput(double(effectiveTemperature));
			double integratedThroughputQE_CHEOPS = wavelengthDependence->getIntegratedThroughputQE(double(effectiveTemperature),SatelliteData::kDefaultCcdTemperature);
			double electronsPerPhoton = integratedThroughputQE_CHEOPS/integratedThroughput_CHEOPS;

			fluxConversion->setCellTEff(effectiveTemperature);
			fluxConversion->setCellCheopsmagMinusGmag(deltaMag);
			fluxConversion->setCellElectronsPerPhoton(electronsPerPhoton);
			fluxConversion->WriteRow();

		}

		delete fluxConv;
		delete fluxConversion;

	} catch(std::exception& e) {
		throw runtime_error(e.what());
	}

    return  0;
}

