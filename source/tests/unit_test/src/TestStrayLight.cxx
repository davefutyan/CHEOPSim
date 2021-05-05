#define BOOST_TEST_MODULE TestCHEOPSim
#include "boost/test/unit_test.hpp"
#include "boost/test/unit_test_suite.hpp"
#include "boost/test/unit_test_log.hpp"
#include "boost/test/test_tools.hpp"
#include "boost/test/detail/global_typedef.hpp"
#include "boost/filesystem.hpp"

#include <cmath>
#include <pqxx/pqxx>
#include <stdlib.h>

#include "LeapSeconds.hxx"
#include "ObtUtcCorrelation.hxx"
#include "BarycentricOffset.hxx"
#include "REF_APP_FluxConversion.hxx"

#include "source/include/StarProducer.hxx"
#include "source/include/TransitFluxModulator.hxx"
#include "source/include/StellarNoiseFluxModulator.hxx"
#include "source/include/StellarVariationFluxModulator.hxx"
#include "source/include/ZodiacalLightGenerator.hxx"
#include "source/include/StrayLightGenerator.hxx"
#include "telescope/include/FocalPlaneGenerator.hxx"
#include "telescope/include/PSFGenerator.hxx"

#define private public
#include "ProgramParams.hxx"

using namespace boost;
using namespace boost::unit_test;

/// @brief unit test fixture
struct StrayLightFixture {
	StrayLightFixture() {
		m_params = CheopsInit(framework::master_test_suite().argc,
						      framework::master_test_suite().argv);

		// Get the time configuration
		boost::posix_time::ptime startTime = boost::posix_time::time_from_string(m_params->GetAsString("startTime"));
		unsigned numberOfStackedImages = m_params->GetAsInt("numberOfStackedImages");
		unsigned exposuresPerStack = m_params->GetAsInt("exposuresPerStack");
		double exposureTime = m_params->GetAsDouble("exposureTime");
		unsigned maxImagesPerCube = m_params->GetAsInt("maxImagesPerCube");

		// Initialize the time configuration
		TimeConfiguration timeConf(startTime,numberOfStackedImages,exposuresPerStack,exposureTime,maxImagesPerCube,m_params->GetAsBool("doFullFrame"));
		// Initialize the data
		m_data = new Data(timeConf,m_params->GetAsString("modulesToRun"));

	    //Create the output directory
		boost::filesystem::create_directory(m_data->getOutputDirectory());

		//Initialize the field of view by defining the pointing direction
		string pointingRA_str = m_params->GetAsString("pointingRA");
		string pointingDec_str = m_params->GetAsString("pointingDec");
		double pointingRA = StarProducer::convertToArcseconds(pointingRA_str,true);
		double pointingDec = StarProducer::convertToArcseconds(pointingDec_str,false);
		SkyPosition pointingDirection(pointingRA,pointingDec);
		SkyFieldOfView *fov = new SkyFieldOfView(pointingDirection);
		m_data->setFieldOfView(fov);

		//Initialize the wavelength dependence
		bool applyThroughput = m_params->GetAsBool("applyThroughput");
		string throughputFilename = m_params->GetAsString("throughputFilename");
		bool applyQE = m_params->GetAsBool("applyQE");
		string QEFilename = m_params->GetAsString("QEFilename");
		string SEDFilename = m_params->GetAsString("SEDFilename");
		string targetSEDFilename = m_params->GetAsString("targetSEDFilename");
		m_data->setWavelengthDependence(new WavelengthDependence(applyThroughput,throughputFilename,applyQE,QEFilename,SEDFilename,targetSEDFilename));

		//Set the reference files needed to calculate the transit time in BJD
		list<string> leapSecondFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU1972-01-01T00-00-00_SOC_APP_LeapSeconds_V0002.fits"};
		LeapSeconds::setLeapSecondsFileNames(leapSecondFile);
		list<string> correlationFile {string(getenv("CHEOPS_SW"))+"/resources/CH_TU2015-01-01T00-00-00_AUX_RES_ObtUtcCorrelation_V000.fits"};
		OBTUTCCorrelation::SetCorrelationFileNames(correlationFile);
		BarycentricOffset::setEphemerisFileName(string(getenv("CHEOPS_TESTDATA"))+"/resources/CH_TU1949-12-14T00-00-00_EXT_APP_DE1_V0100.fits");

	}
	~StrayLightFixture() {
		delete m_data;
	}
	ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	Data * m_data; ///< Pointer to the data
};

BOOST_FIXTURE_TEST_SUITE( TestCHEOPSim, StrayLightFixture )

BOOST_AUTO_TEST_CASE( testStrayLight )
{
	Module * strayLightGenerator = new StrayLightGenerator();
	strayLightGenerator->initialize(m_params->Module("StrayLightGenerator"));
	strayLightGenerator->doBegin(m_data);

	Module * focalPlaneGenerator = new FocalPlaneGenerator();
	focalPlaneGenerator->initialize(m_params->Module("FocalPlaneGenerator"));
	focalPlaneGenerator->doBegin(m_data);

	unsigned nTimesteps = m_data->getTimeConfiguration().getNumberOfTimeSteps();
	for (unsigned iTimeStep=0; iTimeStep<nTimesteps;  iTimeStep++) {
		focalPlaneGenerator->process(m_data,iTimeStep);
		strayLightGenerator->process(m_data,iTimeStep);
	}

	BOOST_CHECK_CLOSE( (m_data->getImages()[8])->getPixelValue(512,512), 0.001686053, 1E-3);

	delete focalPlaneGenerator;
	delete strayLightGenerator;
}

BOOST_AUTO_TEST_SUITE_END()
